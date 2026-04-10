import numpy as np
import argparse

def read_tetgen_node_file(file_path):
    """
    Reads a TetGen.node file and returns a numpy array of shape (num_points, dim) containing the coordinates of the points.

    :param file_path: The path to the.node file.
    :return: A numpy array of shape (num_points, dim) containing the coordinates of the points.
    """
    
    with open(file_path, 'r') as f:
        lines = f.readlines()

    num_points = int(lines[0].strip().split()[0])
    dim = int(lines[0].strip().split()[1])

    points = np.zeros((num_points, dim))
    matids = np.zeros((num_points), dtype=int)
    for i, line in enumerate(lines[1:]):
        data = line.strip().split()
        # print(i)
        points[i] = np.array(data[1:4]).astype(np.float64)
        if len(data) > 4:
            matids[i] = np.array(data[5]).astype(np.int32)
        

    return points, matids


def read_tetgen_ele_file(file_path):
    """
    Reads a TetGen.ele file and returns a numpy array of shape (num_elements, num_vertices) containing the indices of the vertices of each element.

    :param file_path: The path to the.ele file.
    :return: A numpy array of shape (num_elements, num_vertices) containing the indices of the vertices of each element.
    """
    
    with open(file_path, 'r') as f:
        lines = f.readlines()

    num_elements = int(lines[0].strip().split()[0])
    num_vertices = int(lines[0].strip().split()[1])

    elements = np.zeros((num_elements, num_vertices), dtype=int)

    for i, line in enumerate(lines[1:-1]):
        data = line.strip().split()
        elements[i] = np.array(data[1:5]).astype(np.int32)

    return elements

def write_to_avs(nodes, elements, matids, file_path):
    """
    Writes a numpy array of nodes and elements to AVS format.

    :param nodes: A numpy array of shape (num_nodes, 3) containing the coordinates of the nodes.
    :param elements: A numpy array of shape (num_elements, 4) containing the indices of the vertices of each element.
    :param file_path: The path to the output file.
    """
    
    with open(file_path, 'w') as f:
        f.write(f'{nodes.shape[0]} {elements.shape[0]} 1 0 0\n')
        for i in range(nodes.shape[0]):
            f.write(f'{i+1} {nodes[i,0]} {nodes[i,1]} {nodes[i,2]}\n')

        for i in range(elements.shape[0]):
            f.write(f'{i+1} 1 tet {elements[i,0]} {elements[i,1]} {elements[i,3]} {elements[i,2]}\n')

        f.write(f'001 1\n')
        f.write(f'materialid, integer\n')            
        for i in range(matids.shape[0]):
            f.write(f'{i+1} {matids[i]}\n')
            
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Convert TetGen files to AVS format.')
    parser.add_argument('node_file', type=str, help='The path to the TetGen.node file.')
    parser.add_argument('ele_file', type=str, help='The path to the TetGen.ele file.')
    parser.add_argument('output_file', type=str, help='The path to the output AVS file.')

    args = parser.parse_args()

    nodes, matids = read_tetgen_node_file(args.node_file)
    elements = read_tetgen_ele_file(args.ele_file)

    write_to_avs(nodes, elements, matids, args.output_file)