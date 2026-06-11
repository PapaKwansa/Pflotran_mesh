#!/usr/bin/env python3
"""
Build a PFLOTRAN surrogate-training dataset using Latin hypercube sampling.

What it does:
1) Samples log10(permeability) values for the four geologic units.
2) Copies the injection and recovery decks into a sample-specific run folder.
3) Replaces the PERM_ISO values in both decks for that sample.
4) Runs PFLOTRAN for injection, then recovery.
5) Extracts wellbore pressure and geomechanics observables from HDF5 outputs.
6) Saves a compressed master dataset plus a CSV manifest.

Assumptions:
- You already have working decks:
    layers4_injection_geomech.in
    layers4_recovery_geomech.in
- The decks write:
    flow output  -> pflotran.h5
    geomech output -> pflotran-geomech.h5
- The geomechanics output times in the injection/recovery decks are already
  set the way you want for training data.
- The model directory contains the static mesh / region files:
    layers4.uge, layers4.mapping, layers4_material_ids.h5,
    top.ex, bottom.ex, north.ex, south.ex, east.ex, west.ex,
    overburden.vset, bartlesville_sand.vset, basal_layer.vset,
    underburden.vset, wellbore.vset
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
from typing import Dict, List, Tuple

import h5py
import numpy as np
from scipy.stats import qmc


MATERIALS = ["overburden", "bartlesville_sand", "basal_layer", "underburden"]

# Your wellbore HDF5 indices, as provided earlier.
WELLBORE_H5_INDICES = np.array([
    354057, 354058, 354059, 354060, 354061, 354062, 354063, 354064,
    375524, 375525, 375526, 375527, 375528, 375529, 375530, 375531,
    437633, 437634, 437635, 437636, 437637, 437638, 437639, 437640,
    471238, 471239, 471240, 471241, 471242, 471243, 471244, 471245
], dtype=int)

STRAIN_COMPONENTS = [
    "strain_xx", "strain_yy", "strain_zz",
    "strain_xy", "strain_yz", "strain_zx"
]

# Candidate dataset names that PFLOTRAN may use.
PRESSURE_DATASET_CANDIDATES = [
    "Liquid Pressure [Pa]",
    "LIQUID_PRESSURE",
    "Liquid Pressure",
]

# Default permeability bounds in log10(m^2).
# These are intentionally modest starting ranges: +/- 1 order of magnitude
# around your current values. Tighten or widen as needed.
DEFAULT_LOG10_BOUNDS = {
    "overburden": (-18.0, -16.0),
    "bartlesville_sand": (-14.0, -12.0),
    "basal_layer": (-19.0, -17.0),
    "underburden": (-18.0, -16.0),
}

DEFAULT_BASE_PERM = {
    "overburden": 5.0e-18,
    "bartlesville_sand": 1.0e-13,
    "basal_layer": 2.0e-18,
    "underburden": 5.0e-18,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate a PFLOTRAN surrogate dataset with LHS.")
    p.add_argument("--model-dir", type=str, default=".", help="Directory containing the PFLOTRAN input files.")
    p.add_argument("--out-dir", type=str, default="./surrogate_dataset", help="Output dataset directory.")
    p.add_argument("--n-samples", type=int, default=20, help="Number of LHS samples.")
    p.add_argument("--seed", type=int, default=1234, help="Random seed.")
    p.add_argument("--nprocs", type=int, default=int(os.environ.get("SLURM_NTASKS", "64")),
                   help="MPI ranks to use per PFLOTRAN run.")
    p.add_argument("--pflotran-bin", type=str,
                   default=os.environ.get("PFLOTRAN_BIN", "pflotran"),
                   help="Path to PFLOTRAN executable.")
    p.add_argument("--mpiexec", type=str, default="mpiexec", help="MPI launcher command.")
    p.add_argument("--inj-template", type=str, default="layers4_injection_geomech.in",
                   help="Injection deck template filename.")
    p.add_argument("--rec-template", type=str, default="layers4_recovery_geomech.in",
                   help="Recovery deck template filename.")
    p.add_argument("--copy-static", action="store_true",
                   help="Copy static files instead of symlinking them.")
    return p.parse_args()


def safe_unlink(path: Path) -> None:
    if path.is_symlink() or path.exists():
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()


def make_link_or_copy(src: Path, dst: Path, copy_mode: bool) -> None:
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


def replace_perm_iso_in_block(text: str, material_name: str, perm_value: float) -> str:
    """
    Replace the PERM_ISO value inside a MATERIAL_PROPERTY <material_name> block.
    """
    pattern = rf"(MATERIAL_PROPERTY\s+{re.escape(material_name)}.*?PERM_ISO\s+)([0-9eE+\-\.]+)(\s*/)"
    m = re.search(pattern, text, flags=re.DOTALL)
    if not m:
        raise RuntimeError(f"Could not find PERM_ISO for material '{material_name}'.")

    def repl(match: re.Match) -> str:
        return f"{match.group(1)}{perm_value:.6e}{match.group(3)}"

    new_text, n = re.subn(pattern, repl, text, count=1, flags=re.DOTALL)
    if n != 1:
        raise RuntimeError(f"Expected exactly one replacement for material '{material_name}', got {n}.")
    return new_text


def patch_deck_permeabilities(template_path: Path, out_path: Path, k_map: Dict[str, float]) -> None:
    text = template_path.read_text(encoding="utf-8")
    for material in MATERIALS:
        text = replace_perm_iso_in_block(text, material, k_map[material])
    write_text(out_path, text)


def generate_lhs_log10_samples(
    n_samples: int,
    bounds_log10: Dict[str, Tuple[float, float]],
    seed: int
) -> Tuple[np.ndarray, List[str]]:
    names = MATERIALS[:]
    lower = np.array([bounds_log10[m][0] for m in names], dtype=float)
    upper = np.array([bounds_log10[m][1] for m in names], dtype=float)

    sampler = qmc.LatinHypercube(d=len(names), seed=seed)
    unit = sampler.random(n=n_samples)
    scaled = qmc.scale(unit, lower, upper)
    return scaled, names


def find_time_groups(h5obj: h5py.File, dataset_candidates: List[str]) -> List[Tuple[float, str]]:
    """
    Find groups that look like time snapshots and contain one of the candidate datasets.
    Returns [(time_hours, group_path), ...].
    """
    groups: List[Tuple[float, str]] = []

    def visitor(name: str, obj) -> None:
        if not isinstance(obj, h5py.Group):
            return
        if "Time" not in name:
            return

        # Look for any candidate dataset directly inside this group.
        has_candidate = any(candidate in obj for candidate in dataset_candidates)
        if not has_candidate:
            return

        group_name = Path(name).name
        t = parse_time_from_group_name(group_name)
        if t is not None:
            groups.append((t, name))

    h5obj.visititems(visitor)

    # Deduplicate by group path and sort by time.
    unique = {}
    for t, path in groups:
        unique[path] = t
    out = sorted([(t, p) for p, t in unique.items()], key=lambda x: x[0])
    return out


def parse_time_from_group_name(group_name: str) -> float | None:
    """
    Parse strings like:
      '0 Time 1.000000e+00 h'
      '1 Time 19.0 h'
    """
    m = re.search(r"Time\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s*h", group_name)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def find_dataset_in_group(group: h5py.Group, candidates: List[str]) -> np.ndarray:
    """
    Return the first matching dataset inside a group (searching recursively).
    """
    found = None

    def visitor(name: str, obj) -> None:
        nonlocal found
        if found is not None:
            return
        if isinstance(obj, h5py.Dataset):
            leaf = Path(name).name
            if leaf in candidates:
                found = np.asarray(obj, dtype=float)

    group.visititems(visitor)

    if found is None:
        raise KeyError(f"None of the candidate datasets were found: {candidates}")
    return found


def compute_well_stats_at_time(
    arr: np.ndarray,
    well_idx: np.ndarray
) -> Dict[str, float]:
    values = arr[well_idx]
    return {
        "median": float(np.nanmedian(values)),
        "min": float(np.nanmin(values)),
        "max": float(np.nanmax(values)),
    }


def extract_pressure_series(h5_path: Path, well_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract wellbore pressure statistics vs time from a PFLOTRAN flow HDF5 file.
    """
    with h5py.File(h5_path, "r") as f:
        groups = find_time_groups(f, PRESSURE_DATASET_CANDIDATES)
        if not groups:
            raise RuntimeError(f"No pressure time groups found in {h5_path}")

        times = []
        med = []
        pmin = []
        pmax = []

        for t, group_path in groups:
            grp = f[group_path]
            pressure = find_dataset_in_group(grp, PRESSURE_DATASET_CANDIDATES)
            if well_idx.max() >= len(pressure):
                raise IndexError(
                    f"Well index out of bounds for pressure array in {h5_path} at time {t} h."
                )

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
    """
    Extract median strain tensor components and volumetric strain at the wellbore
    from a PFLOTRAN geomechanics HDF5 file.
    """
    with h5py.File(h5_path, "r") as f:
        groups = find_time_groups(f, STRAIN_COMPONENTS)
        if not groups:
            raise RuntimeError(f"No geomechanics time groups found in {h5_path}")

        times = []
        med_strains = []
        vol_strain = []

        for t, group_path in groups:
            grp = f[group_path]
            comp_vals = []
            for comp in STRAIN_COMPONENTS:
                if comp not in grp:
                    raise KeyError(f"Component '{comp}' not found in {h5_path} group '{group_path}'")
                arr = np.asarray(grp[comp], dtype=float)
                if well_idx.max() >= len(arr):
                    raise IndexError(
                        f"Well index out of bounds for strain array in {h5_path} at time {t} h."
                    )
                comp_vals.append(np.nanmedian(arr[well_idx]))

            comp_vals = np.asarray(comp_vals, dtype=float)
            times.append(t)
            med_strains.append(comp_vals)
            vol_strain.append(comp_vals[0] + comp_vals[1] + comp_vals[2])

    return (
        np.asarray(times, dtype=float),
        np.asarray(med_strains, dtype=float),   # shape: (ntime, 6)
        np.asarray(vol_strain, dtype=float),     # shape: (ntime,)
    )


def merge_time_series(
    times_a: np.ndarray,
    data_a: np.ndarray,
    times_b: np.ndarray,
    data_b: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Merge two time series arrays and sort by time.
    If there are duplicate times, keep the first occurrence.
    """
    times = np.concatenate([times_a, times_b], axis=0)
    data = np.concatenate([data_a, data_b], axis=0)

    order = np.argsort(times)
    times = times[order]
    data = data[order]

    uniq_times, uniq_idx = np.unique(times, return_index=True)
    return uniq_times, data[uniq_idx]


def find_restart_checkpoint(run_dir: Path) -> Path:
    """
    Find the 19 h restart checkpoint produced by the injection run and return it.
    We expect an HDF5 checkpoint, but a .chk fallback is included.
    """
    h5_candidates = sorted(run_dir.glob("*19*hour*.h5")) + sorted(run_dir.glob("*19*.h5"))
    chk_candidates = sorted(run_dir.glob("*19*hour*.chk")) + sorted(run_dir.glob("*19*.chk"))

    if h5_candidates:
        return h5_candidates[0]
    if chk_candidates:
        return chk_candidates[0]

    raise FileNotFoundError(f"No 19-hour checkpoint found in {run_dir}")


def run_pflotran(run_dir: Path, pflotran_bin: str, mpiexec: str, nprocs: int) -> None:
    cmd = [mpiexec, "-n", str(nprocs), pflotran_bin]
    subprocess.run(cmd, cwd=run_dir, check=True)


def prepare_sample_run_dir(
    model_dir: Path,
    run_root: Path,
    sample_id: int,
    k_map: Dict[str, float],
    inj_template_name: str,
    rec_template_name: str,
    copy_static: bool
) -> Path:
    """
    Create a sample-specific run directory and populate it with symlinks/copies
    and deck files patched with the sample permeability values.
    """
    sample_dir = run_root / f"sample_{sample_id:04d}"
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Static files needed by your model.
    static_files = [
        "layers4.uge",
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
        make_link_or_copy(src, sample_dir / fname, copy_static)

    # Patch and write the two decks.
    inj_src = model_dir / inj_template_name
    rec_src = model_dir / rec_template_name
    if not inj_src.exists():
        raise FileNotFoundError(f"Missing injection template: {inj_src}")
    if not rec_src.exists():
        raise FileNotFoundError(f"Missing recovery template: {rec_src}")

    inj_text = inj_src.read_text(encoding="utf-8")
    rec_text = rec_src.read_text(encoding="utf-8")

    for material, perm in k_map.items():
        inj_text = replace_perm_iso_in_block(inj_text, material, perm)
        rec_text = replace_perm_iso_in_block(rec_text, material, perm)

    write_text(sample_dir / inj_template_name, inj_text)
    write_text(sample_dir / rec_template_name, rec_text)

    return sample_dir


def read_model_outputs_for_sample(
    sample_dir: Path,
    well_idx: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Read injection + recovery outputs from the sample directory and merge them.
    """
    inj_flow = sample_dir / "inj_pflotran.h5"
    inj_geo = sample_dir / "inj_pflotran-geomech.h5"
    rec_flow = sample_dir / "rec_pflotran.h5"
    rec_geo = sample_dir / "rec_pflotran-geomech.h5"

    # Read injection outputs
    t_p_inj, p_med_inj, p_min_inj, p_max_inj = extract_pressure_series(inj_flow, well_idx)
    t_s_inj, s_med_inj, ev_inj = extract_geomech_series(inj_geo, well_idx)

    # Read recovery outputs
    t_p_rec, p_med_rec, p_min_rec, p_max_rec = extract_pressure_series(rec_flow, well_idx)
    t_s_rec, s_med_rec, ev_rec = extract_geomech_series(rec_geo, well_idx)

    # Merge pressure and geomechanics time series.
    t_p, p_med = merge_time_series(t_p_inj, p_med_inj, t_p_rec, p_med_rec)
    _, p_min = merge_time_series(t_p_inj, p_min_inj, t_p_rec, p_min_rec)
    _, p_max = merge_time_series(t_p_inj, p_max_inj, t_p_rec, p_max_rec)

    t_s, s_med = merge_time_series(t_s_inj, s_med_inj, t_s_rec, s_med_rec)
    _, ev = merge_time_series(t_s_inj, ev_inj, t_s_rec, ev_rec)

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

    # Bounds in log10 space.
    bounds = DEFAULT_LOG10_BOUNDS

    # Generate LHS samples in log10(k).
    lhs_log10, names = generate_lhs_log10_samples(args.n_samples, bounds, args.seed)

    # Collect dataset in memory.
    k_log10_all = []
    k_all = []
    pressure_times_ref = None
    strain_times_ref = None
    pressure_median_all = []
    pressure_min_all = []
    pressure_max_all = []
    strain_median_all = []
    volumetric_strain_all = []

    manifest_rows = []
    failures = []

    for i in range(args.n_samples):
        sample_id = i + 1
        sample_dir = prepare_sample_run_dir(
            model_dir=model_dir,
            run_root=run_root,
            sample_id=sample_id,
            k_map={},  # filled below
            inj_template_name=args.inj_template,
            rec_template_name=args.rec_template,
            copy_static=args.copy_static,
        )

        # Convert sample values to permeability map.
        sample_log10 = lhs_log10[i]
        k_map = {mat: float(10.0 ** sample_log10[j]) for j, mat in enumerate(names)}

        # Re-patch decks with actual k values.
        for deck_name in [args.inj_template, args.rec_template]:
            deck_path = sample_dir / deck_name
            text = deck_path.read_text(encoding="utf-8")
            for mat in MATERIALS:
                text = replace_perm_iso_in_block(text, mat, k_map[mat])
            write_text(deck_path, text)

        # Run injection.
        try:
            inj_deck = sample_dir / args.inj_template
            rec_deck = sample_dir / args.rec_template

            # Make PFLOTRAN read the injection deck.
            make_link_or_copy(inj_deck, sample_dir / "pflotran.in", copy_mode=False)
            run_pflotran(sample_dir, args.pflotran_bin, args.mpiexec, args.nprocs)

            # Keep a copy of the injection outputs before they get overwritten.
            inj_flow = sample_dir / "pflotran.h5"
            inj_geo = sample_dir / "pflotran-geomech.h5"
            if not inj_flow.exists():
                raise FileNotFoundError(f"Injection flow output missing: {inj_flow}")
            if not inj_geo.exists():
                raise FileNotFoundError(f"Injection geomech output missing: {inj_geo}")

            safe_unlink(sample_dir / "inj_pflotran.h5")
            safe_unlink(sample_dir / "inj_pflotran-geomech.h5")
            os.symlink(inj_flow.resolve(), sample_dir / "inj_pflotran.h5")
            os.symlink(inj_geo.resolve(), sample_dir / "inj_pflotran-geomech.h5")

            # Find the restart checkpoint and make the recovery deck use it.
            restart_file = find_restart_checkpoint(sample_dir)
            restart_link = sample_dir / "layers4-19hour-restart.h5"
            safe_unlink(restart_link)
            os.symlink(restart_file.resolve(), restart_link)

            # Run recovery.
            make_link_or_copy(rec_deck, sample_dir / "pflotran.in", copy_mode=False)
            run_pflotran(sample_dir, args.pflotran_bin, args.mpiexec, args.nprocs)

            rec_flow = sample_dir / "pflotran.h5"
            rec_geo = sample_dir / "pflotran-geomech.h5"
            if not rec_flow.exists():
                raise FileNotFoundError(f"Recovery flow output missing: {rec_flow}")
            if not rec_geo.exists():
                raise FileNotFoundError(f"Recovery geomech output missing: {rec_geo}")

            safe_unlink(sample_dir / "rec_pflotran.h5")
            safe_unlink(sample_dir / "rec_pflotran-geomech.h5")
            os.symlink(rec_flow.resolve(), sample_dir / "rec_pflotran.h5")
            os.symlink(rec_geo.resolve(), sample_dir / "rec_pflotran-geomech.h5")

            # Extract observables.
            obs = read_model_outputs_for_sample(sample_dir, WELLBORE_H5_INDICES)

            # Keep a reference time grid across samples.
            if pressure_times_ref is None:
                pressure_times_ref = obs["pressure_times"]
            elif not np.allclose(pressure_times_ref, obs["pressure_times"]):
                raise RuntimeError("Pressure time grid changed across samples. Check output times.")

            if strain_times_ref is None:
                strain_times_ref = obs["strain_times"]
            elif not np.allclose(strain_times_ref, obs["strain_times"]):
                raise RuntimeError("Geomechanics time grid changed across samples. Check output times.")

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

    # Convert to arrays.
    k_log10_all = np.asarray(k_log10_all, dtype=float)                 # (n_ok, 4)
    k_all = np.asarray(k_all, dtype=float)                             # (n_ok, 4)
    pressure_median_all = np.asarray(pressure_median_all, dtype=float) # (n_ok, ntp)
    pressure_min_all = np.asarray(pressure_min_all, dtype=float)       # (n_ok, ntp)
    pressure_max_all = np.asarray(pressure_max_all, dtype=float)       # (n_ok, ntp)
    strain_median_all = np.asarray(strain_median_all, dtype=float)     # (n_ok, nts, 6)
    volumetric_strain_all = np.asarray(volumetric_strain_all, dtype=float)  # (n_ok, nts)

    # Save master dataset.
    np.savez_compressed(
        out_dir / "dataset_master.npz",
        material_names=np.array(MATERIALS, dtype="U"),
        pressure_times=pressure_times_ref,
        strain_times=strain_times_ref,
        k_log10=k_log10_all,
        k_values=k_all,
        pressure_median=pressure_median_all,
        pressure_min=pressure_min_all,
        pressure_max=pressure_max_all,
        strain_median=strain_median_all,
        volumetric_strain=volumetric_strain_all,
    )

    # Save a CSV manifest for bookkeeping.
    manifest_path = out_dir / "sample_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sample_id", "status",
            "overburden_k", "bartlesville_sand_k", "basal_layer_k", "underburden_k",
            "run_dir",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    # Save metadata.
    meta = {
        "n_requested": args.n_samples,
        "n_successful": int(len(k_log10_all)),
        "n_failed": int(len(failures)),
        "failures": [{"sample_id": sid, "error": err} for sid, err in failures],
        "materials": MATERIALS,
        "pressure_times": pressure_times_ref.tolist() if pressure_times_ref is not None else None,
        "strain_times": strain_times_ref.tolist() if strain_times_ref is not None else None,
        "wellbore_h5_indices": WELLBORE_H5_INDICES.tolist(),
        "base_permeability_values": DEFAULT_BASE_PERM,
        "log10_bounds": DEFAULT_LOG10_BOUNDS,
        "notes": [
            "Pressure and geomechanics outputs are stored separately by PFLOTRAN.",
            "The dataset contains median wellbore pressure and median wellbore strain tensors.",
        ],
    }
    (out_dir / "dataset_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print("\nDone.")
    print(f"Successful samples: {len(k_log10_all)} / {args.n_samples}")
    print(f"Dataset: {out_dir / 'dataset_master.npz'}")
    print(f"Manifest: {manifest_path}")
    print(f"Metadata: {out_dir / 'dataset_metadata.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())