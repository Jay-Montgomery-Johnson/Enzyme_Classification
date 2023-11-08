import networkx as nx
from protein_properties import aa_index
from GraphUtils import load_PDB_dir, determine_SSE, remove_SSE_with_aa_length_less_than_or_equal_to_n, \
                    create_structural_edges, create_nodes_and_adjacent_edges, get_sequence, create_graph
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import matplotlib.pyplot as plt
import json

with open('dataset/pdb_files_in_database.txt') as f:
    data = f.read()
pdb_ids = json.loads(data)

count = 0
for EC in range(1, 7):
    for pdb_id in pdb_ids[str(EC)]:
        try:
            structure, residues = load_PDB_dir(pdb_id, f'dataset/ec_{EC}/pdb/')
            SSE = determine_SSE(pdb_id, f'dataset/ec_{EC}/pdb/')
            SSE = remove_SSE_with_aa_length_less_than_or_equal_to_n(SSE, 2, structure, residues)
            structural_edges, structural_distances = create_structural_edges(structure, SSE)
            nodes, sequential_edges, sequential_edge_lengths, sse_lengths = create_nodes_and_adjacent_edges(SSE)
            
            property_calc = aa_index()
            properties = []
            for i in SSE:
                properties.append(property_calc.get_properties(get_sequence(i, structure, residues)))
            
            graph = create_graph(nodes, sequential_edges, structural_edges, structural_distances, 
                sequential_edge_lengths, sse_lengths, properties)
            
            print(f'EC_{pdb_id}:{count}')
            nx.write_gexf(graph, f'dataset/ec_{EC}/graph/' + pdb_id + '.gexf')
            count += 1
        except:
            print('error')
            continue
