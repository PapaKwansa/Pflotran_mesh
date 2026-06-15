#!/usr/bin/env python3
"""Build a PFLOTRAN surrogate-training dataset from a single continuous run.

This version matches the current workflow:
- one coupled PFLOTRAN deck per sample (injection 0-19 h, shut-in 19-96 h)
- no checkpoint/restart handoff
- one flow HDF5 output and one geomechanics HDF5 output per sample
- Latin hypercube sampling over the four layer permeabilities

The script expects a template deck similar to `layers4_geomech.in` and patches
its MATERIAL_PROPERTY PERM_ISO values for each sample.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import h5py
import numpy as np
from scipy.stats import qmc


MATERIALS = ["overburden", "bartlesville_sand", "basal_layer", "underburden"]

# Wellbore indices previously identified for the North Avant model.
WELLBORE_H5_INDICES = np.array([
    354057, 354058, 354059, 354060, 354061, 354062, 354063, 354064,
    375524, 375525, 375526, 375527, 375528, 375529, 375530, 375531,
    437633, 437634, 437635, 437636, 437637, 437638, 437639, 437640,
    471238, 471239, 471240, 471241, 471242, 471243, 471244, 471245,
], dtype=int)

STRAIN_COMPONENTS = [
    "strain_xx", "strain_yy", "strain_zz",
    "strain_xy", "strain_yz", "strain_zx",
]

PRESSURE_DATASET_CANDIDATES = [
    "LIQUID_PRESSURE",
    "Liquid Pressure [Pa]",
    "Liquid Pressure",
]

DEFAULT_LOG10_BOUNDS: Dict[str, Tuple[float, float]] = {
    "overburden": (-18.0, -16.0),
    "bartlesville_sand": (-14.0, -12.0),
    "basal_layer": (-19.0, -17.0),
    "underburden": (-18.0, -16.0),
}

DEFAULT_BASE_PERM: Dict[str, float] = {
    "overburden": 5.0e-18,
    "bartlesville_sand": 1.0e-13,
    "basal_layer": 2.0e-18,
    "underburden": 5.0e-18,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Generate a PFLOTRAN surrogate dataset from a single continuous "
            "injection->shut-in run."
        )
    )
    p.add_argument("--model-dir", type=str, default=".",
                   help="Directory containing the PFLOTRAN input files.")
    p.add_argument("--out-dir", type=str, default="./surrogate_dataset",
                   help="Output dataset directory.")
    p.add_argument("--n-samples", type=int, default=20,
                   help="Number of Latin hypercube samples.")
    p.add_argument("--seed", type=int, default=1234,
                   help="Random seed for the LHS sampler.")
    p.add_argument("--nprocs", type=int, default=int(os.environ.get("SLURM_NTASKS", "64")),
                   help="MPI ranks to use per PFLOTRAN run.")
    p.add_argument("--pflotran-bin", type=str,
                   default=os.environ.get("PFLOTRAN_BIN", "pflotran"),
                   help="Path to the PFLOTRAN executable.")
    p.add_argument("--mpiexec", type=str, default="mpiexec",
                   help="MPI launcher command.")
    p.add_argument("--deck-template", type=str, default="layers4_geomech.in",
                   help="Template deck filename to patch for each sample.")
    p.add_argument("--copy-static", action="store_true",
                   help="Copy static files instead of symlinking them.")
    p.add_argument("--keep-runs", action="store_true",
                   help="Keep sample run directories even after successful extraction.")
    return p.parse_args()


def safe_unlink(path: Path) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def link_or_copy(src: Path, dst: Path, copy_mode: bool) -> None:
    safe_unlink(dst)
    if copy_mode:
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
    else:
        os.symlink(src.resolve(), dst)


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def replace_perm_iso_in_block(text: str, material_name: str, perm_value: float) -> str:
    """Replace PERM_ISO in the MATERIAL_PROPERTY block for one material."""
    pattern = rf"(MATERIAL_PROPERTY\s+{re.escape(material_name)}.*?PERM_ISO\s+)([0-9eEdD+\-.]+)(\s*/)"
    match = re.search(pattern, text, flags=re.DOTALL)
    if not match:
        raise RuntimeError(f"Could not find PERM_ISO for material '{material_name}'.")

    def repl(m: re.Match) -> str:
        return f"{m.group(1)}{perm_value:.6e}{m.group(3)}"

    new_text, n = re.subn(pattern, repl, text, count=1, flags=re.DOTALL)
    if n != 1:
        raise RuntimeError(f"Expected exactly one replacement for material '{material_name}', got {n}.")
    return new_text


def generate_lhs_log10_samples(
    n_samples: int,
    bounds_log10: Dict[str, Tuple[float, float]],
    seed: int,
) -> Tuple[np.ndarray, List[str]]:
    names = MATERIALS[:]
    lower = np.array([bounds_log10[m][0] for m in names], dtype=float)
    upper = np.array([bounds_log10[m][1] for m in names], dtype=float)
    sampler = qmc.LatinHypercube(d=len(names), seed=seed)
    unit = sampler.random(n=n_samples)
    scaled = qmc.scale(unit, lower, upper)
    return scaled, names


def parse_time_from_group_name(group_name: str) -> float | None:
    """Parse strings like '0 Time 1.000000e+00 h'."""
    m = re.search(r"Time\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s*h", group_name)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def find_time_groups(h5obj: h5py.File, dataset_candidates: Sequence[str]) -> List[Tuple[float, str]]:
    """Return [(time_hours, group_path), ...] for time groups containing candidate datasets."""
    groups: List[Tuple[float, str]] = []
    norm_cands = [normalize_name(c) for c in dataset_candidates]

    def visitor(name: str, obj) -> None:
        if not isinstance(obj, h5py.Group):
            return
        if "Time" not in name:
            return

        has_candidate = False
        for leaf_name, _ in obj.items():
            if normalize_name(leaf_name) in norm_cands:
                has_candidate = True
                break
        if not has_candidate:
            # search recursively one level deeper for names that contain candidates
            def leaf_visitor(subname: str, subobj) -> None:
                nonlocal has_candidate
                if has_candidate:
                    return
                if isinstance(subobj, h5py.Dataset):
                    leaf = normalize_name(Path(subname).name)
                    if any(c == leaf or c in leaf or leaf in c for c in norm_cands):
                        has_candidate = True
            obj.visititems(leaf_visitor)

        if not has_candidate:
            return

        t = parse_time_from_group_name(Path(name).name)
        if t is not None:
            groups.append((t, name))

    h5obj.visititems(visitor)

    dedup: Dict[str, float] = {}
    for t, p in groups:
        dedup[p] = t
    return sorted([(t, p) for p, t in dedup.items()], key=lambda x: x[0])


def find_dataset_in_group(group: h5py.Group, candidates: Sequence[str]) -> np.ndarray:
    """Return the first matching dataset inside a group (recursive search)."""
    norm_cands = [normalize_name(c) for c in candidates]
    found: np.ndarray | None = None

    def visitor(name: str, obj) -> None:
        nonlocal found
        if found is not None:
            return
        if isinstance(obj, h5py.Dataset):
            leaf = normalize_name(Path(name).name)
            if any(c == leaf or c in leaf or leaf in c for c in norm_cands):
                found = np.asarray(obj, dtype=float)

    group.visititems(visitor)
    if found is None:
        raise KeyError(f"None of the candidate datasets were found: {candidates}")
    return found


def compute_well_stats_at_time(arr: np.ndarray, well_idx: np.ndarray) -> Dict[str, float]:
    values = arr[well_idx]
    return {
        "median": float(np.nanmedian(values)),
        "min": float(np.nanmin(values)),
        "max": float(np.nanmax(values)),
    }


def extract_pressure_series(h5_path: Path, well_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as f:
        groups = find_time_groups(f, PRESSURE_DATASET_CANDIDATES)
        if not groups:
            raise RuntimeError(f"No pressure time groups found in {h5_path}")

        times: List[float] = []
        med: List[float] = []
        pmin: List[float] = []
        pmax: List[float] = []

        for t, group_path in groups:
            grp = f[group_path]
            pressure = find_dataset_in_group(grp, PRESSURE_DATASET_CANDIDATES)
            if well_idx.max() >= len(pressure):
                raise IndexError(f"Well index out of bounds for pressure array in {h5_path} at time {t} h.")
            stats = compute_well_stats_at_time(pressure, well_idx)
            times.append(t)
            med.append(stats["median"])
            pmin.append(stats["min"])
            pmax.append(stats["max"])

    return (
        np.asarray(times, dtype=float),
        np.asarray(med, dtype=float),
        np.asarray(pmin, dtype=float),
        np.asarray(pmax, dtype=float),
    )


def extract_geomech_series(h5_path: Path, well_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as f:
        groups = find_time_groups(f, STRAIN_COMPONENTS)
        if not groups:
            raise RuntimeError(f"No geomechanics time groups found in {h5_path}")

        times: List[float] = []
        med_strains: List[np.ndarray] = []
        vol_strain: List[float] = []

        for t, group_path in groups:
            grp = f[group_path]
            comp_vals = []
            for comp in STRAIN_COMPONENTS:
                arr = find_dataset_in_group(grp, [comp])
                if well_idx.max() >= len(arr):
                    raise IndexError(f"Well index out of bounds for strain array in {h5_path} at time {t} h.")
                comp_vals.append(np.nanmedian(arr[well_idx]))

            comp_vals_arr = np.asarray(comp_vals, dtype=float)
            times.append(t)
            med_strains.append(comp_vals_arr)
            vol_strain.append(float(comp_vals_arr[0] + comp_vals_arr[1] + comp_vals_arr[2]))

    return (
        np.asarray(times, dtype=float),
        np.asarray(med_strains, dtype=float),
        np.asarray(vol_strain, dtype=float),
    )


def prepare_sample_run_dir(
    model_dir: Path,
    run_root: Path,
    sample_id: int,
    k_map: Dict[str, float],
    deck_template_name: str,
    copy_static: bool,
) -> Path:
    sample_dir = run_root / f"sample_{sample_id:04d}"
    sample_dir.mkdir(parents=True, exist_ok=True)

    static_files = [
        "layers4.uge",
        "layers4.ugi",
        "layers4.mapping",
        "layers4_material_ids.h5",
        "wellbore.vset",
        "top.ex", "bottom.ex", "north.ex", "south.ex", "east.ex", "west.ex",
        "overburden.vset",
        "bartlesville_sand.vset",
        "basal_layer.vset",
        "underburden.vset",
    ]

    for fname in static_files:
        src = model_dir / fname
        if not src.exists():
            raise FileNotFoundError(f"Missing required input file: {src}")
        link_or_copy(src, sample_dir / fname, copy_static)

    deck_src = model_dir / deck_template_name
    if not deck_src.exists():
        raise FileNotFoundError(f"Missing deck template: {deck_src}")

    deck_text = deck_src.read_text(encoding="utf-8")
    for material, perm in k_map.items():
        deck_text = replace_perm_iso_in_block(deck_text, material, perm)

    # Write the active PFLOTRAN deck as pflotran.in to avoid symlink confusion.
    write_text(sample_dir / "pflotran.in", deck_text)
    return sample_dir


def run_pflotran(run_dir: Path, pflotran_bin: str, mpiexec: str, nprocs: int) -> None:
    cmd = [mpiexec, "-n", str(nprocs), pflotran_bin]
    subprocess.run(cmd, cwd=run_dir, check=True)


def read_sample_outputs(sample_dir: Path, well_idx: np.ndarray) -> Dict[str, np.ndarray]:
    flow_h5 = sample_dir / "pflotran.h5"
    geomech_h5 = sample_dir / "pflotran-geomech.h5"
    if not flow_h5.exists():
        raise FileNotFoundError(f"Flow output missing: {flow_h5}")
    if not geomech_h5.exists():
        raise FileNotFoundError(f"Geomechanics output missing: {geomech_h5}")

    t_p, p_med, p_min, p_max = extract_pressure_series(flow_h5, well_idx)
    t_s, s_med, ev = extract_geomech_series(geomech_h5, well_idx)

    return {
        "pressure_times": t_p,
        "pressure_median": p_med,
        "pressure_min": p_min,
        "pressure_max": p_max,
        "strain_times": t_s,
        "strain_median": s_med,
        "volumetric_strain": ev,
    }


def main() -> int:
    args = parse_args()

    model_dir = Path(args.model_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    run_root = out_dir / "runs"
    run_root.mkdir(parents=True, exist_ok=True)

    lhs_log10, names = generate_lhs_log10_samples(args.n_samples, DEFAULT_LOG10_BOUNDS, args.seed)

    k_log10_all: List[np.ndarray] = []
    k_all: List[np.ndarray] = []
    pressure_times_ref: np.ndarray | None = None
    strain_times_ref: np.ndarray | None = None
    pressure_median_all: List[np.ndarray] = []
    pressure_min_all: List[np.ndarray] = []
    pressure_max_all: List[np.ndarray] = []
    strain_median_all: List[np.ndarray] = []
    volumetric_strain_all: List[np.ndarray] = []

    manifest_rows: List[Dict[str, object]] = []
    failures: List[Tuple[int, str]] = []

    for i in range(args.n_samples):
        sample_id = i + 1
        sample_log10 = lhs_log10[i]
        k_map = {mat: float(10.0 ** sample_log10[j]) for j, mat in enumerate(names)}

        sample_dir = prepare_sample_run_dir(
            model_dir=model_dir,
            run_root=run_root,
            sample_id=sample_id,
            k_map=k_map,
            deck_template_name=args.deck_template,
            copy_static=args.copy_static,
        )

        try:
            run_pflotran(sample_dir, args.pflotran_bin, args.mpiexec, args.nprocs)
            obs = read_sample_outputs(sample_dir, WELLBORE_H5_INDICES)

            if pressure_times_ref is None:
                pressure_times_ref = obs["pressure_times"]
            elif not np.allclose(pressure_times_ref, obs["pressure_times"]):
                raise RuntimeError("Pressure time grid changed across samples. Check deck output times.")

            if strain_times_ref is None:
                strain_times_ref = obs["strain_times"]
            elif not np.allclose(strain_times_ref, obs["strain_times"]):
                raise RuntimeError("Geomechanics time grid changed across samples. Check deck output times.")

            k_log10_all.append(sample_log10)
            k_all.append(np.array([k_map[m] for m in MATERIALS], dtype=float))
            pressure_median_all.append(obs["pressure_median"])
            pressure_min_all.append(obs["pressure_min"])
            pressure_max_all.append(obs["pressure_max"])
            strain_median_all.append(obs["strain_median"])
            volumetric_strain_all.append(obs["volumetric_strain"])

            manifest_rows.append({
                "sample_id": sample_id,
                "status": "ok",
                "overburden_k": k_map["overburden"],
                "bartlesville_sand_k": k_map["bartlesville_sand"],
                "basal_layer_k": k_map["basal_layer"],
                "underburden_k": k_map["underburden"],
                "run_dir": str(sample_dir),
            })
            print(f"[OK] sample {sample_id:04d}")

            if not args.keep_runs:
                # Keep only the extracted results in the master dataset, not the run folder.
                # Comment out the next line if you want to inspect each run directory later.
                shutil.rmtree(sample_dir, ignore_errors=True)

        except Exception as e:
            failures.append((sample_id, str(e)))
            manifest_rows.append({
                "sample_id": sample_id,
                "status": f"failed: {e}",
                "overburden_k": k_map.get("overburden", np.nan),
                "bartlesville_sand_k": k_map.get("bartlesville_sand", np.nan),
                "basal_layer_k": k_map.get("basal_layer", np.nan),
                "underburden_k": k_map.get("underburden", np.nan),
                "run_dir": str(sample_dir),
            })
            print(f"[FAIL] sample {sample_id:04d}: {e}", file=sys.stderr)

    if not k_log10_all:
        raise RuntimeError("No successful samples were generated.")

    k_log10_arr = np.asarray(k_log10_all, dtype=float)
    k_arr = np.asarray(k_all, dtype=float)
    pressure_median_arr = np.asarray(pressure_median_all, dtype=float)
    pressure_min_arr = np.asarray(pressure_min_all, dtype=float)
    pressure_max_arr = np.asarray(pressure_max_all, dtype=float)
    strain_median_arr = np.asarray(strain_median_all, dtype=float)
    volumetric_strain_arr = np.asarray(volumetric_strain_all, dtype=float)

    np.savez_compressed(
        out_dir / "dataset_master.npz",
        material_names=np.array(MATERIALS, dtype="U"),
        pressure_times=pressure_times_ref,
        strain_times=strain_times_ref,
        k_log10=k_log10_arr,
        k_values=k_arr,
        pressure_median=pressure_median_arr,
        pressure_min=pressure_min_arr,
        pressure_max=pressure_max_arr,
        strain_median=strain_median_arr,
        volumetric_strain=volumetric_strain_arr,
    )

    manifest_path = out_dir / "sample_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sample_id", "status",
            "overburden_k", "bartlesville_sand_k", "basal_layer_k", "underburden_k",
            "run_dir",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    meta = {
        "workflow": "single_continuous_run_injection_to_shutin",
        "deck_template": args.deck_template,
        "n_requested": args.n_samples,
        "n_successful": int(len(k_log10_arr)),
        "n_failed": int(len(failures)),
        "failures": [{"sample_id": sid, "error": err} for sid, err in failures],
        "materials": MATERIALS,
        "pressure_times": pressure_times_ref.tolist() if pressure_times_ref is not None else None,
        "strain_times": strain_times_ref.tolist() if strain_times_ref is not None else None,
        "wellbore_h5_indices": WELLBORE_H5_INDICES.tolist(),
        "base_permeability_values": DEFAULT_BASE_PERM,
        "log10_bounds": DEFAULT_LOG10_BOUNDS,
        "notes": [
            "One coupled PFLOTRAN run per sample.",
            "Injection is active from 0 to 19 h, then set to zero through 96 h.",
            "Pressure and geomechanics outputs are extracted from the single run outputs.",
        ],
    }
    (out_dir / "dataset_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print("\nDone.")
    print(f"Successful samples: {len(k_log10_arr)} / {args.n_samples}")
    print(f"Dataset: {out_dir / 'dataset_master.npz'}")
    print(f"Manifest: {manifest_path}")
    print(f"Metadata: {out_dir / 'dataset_metadata.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
