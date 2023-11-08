import networkx as nx
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import Bio.PDB as pdb
from protein_properties import aa_index
import numpy as np
from Bio.PDB.Selection import unfold_entities

def load_PDB_dir(pdb_id, path):
    parser = PDBParser(PERMISSIVE=0)
    structure = parser.get_structure(pdb_id, path + pdb_id + '.pdb')
    residues = unfold_entities(structure[0], 'R')
    #remove waters
    for n, residue in enumerate(residues):
        if residue.get_resname() == 'HOH':
            residues.pop(n)
    return structure, residues

def determine_SSE(pdb_id, path):
    #2. determine SSE
    dssp_tuple = dssp_dict_from_pdb_file( path + pdb_id + '.pdb')
    #3. Parser SSE so they are in a list of tuples
    SSE_list = []
    for key, value in dssp_tuple[0].items():
        SSE_list.append((key, value[1]))
    SSE = []
    e = ''
    prev_res = ''
    for i, n in enumerate(SSE_list):
        if n[1] != e:
            if len(SSE) > 0:
                SSE[-1][2] = prev_res
            SSE.append([n[1], n[0], ''])
        e = n[1]
        prev_res = n[0]
    SSE[-1][2] = prev_res
    return SSE

def remove_SSE_with_aa_length_less_than_or_equal_to_n(SSE, n, structure, residues):
    new_SSE = []
    for i in SSE:
        aa_length = get_aa_len(i, structure, residues)
        if aa_length > n:
            new_SSE.append(i)
    return new_SSE

def residue_midpoint(residue1, residue2):
	coord1 = residue1['CA'].get_coord()
	coord2 = residue2['CA'].get_coord()
	
	midpoint = (coord1 + coord2)/2
	
	return midpoint

def create_structural_edges(structure, SSE):
    SSE_midpoints = []
    for i in SSE:
        try:
            residue1 = structure[0][i[1][0]][i[1][1]]
            residue2 = structure[0][i[2][0]][i[2][1]]
            midpoint = residue_midpoint(residue1, residue2)
            SSE_midpoints.append(midpoint)
        except:
            print('NO STRUCTURAL EDGES')
    
    structural_edges = []
    edges_distance = []
    for n1, mid1 in enumerate(SSE_midpoints):
        s_edges = []
        edge_d = []
        for n2, mid2 in enumerate(SSE_midpoints):
            if n2 < n1 and n2 != n1+1 and n2 != n1-1:
                distance = np.linalg.norm(mid2-mid1)
                s_edges.append([n1, n2])
                edge_d.append(distance)
        
        smallest_3_indeces = sorted(range(len(edge_d)), key = lambda k: edge_d[k])[:3]
        for i in smallest_3_indeces:
            structural_edges.append(s_edges[i])
            edges_distance.append(edge_d[i])
    return structural_edges, edges_distance

def get_sequence2(SSE, structure):
    #here SSE refers to an element of the SSE list
    sequence = []
    start_res = SSE[1]
    end_res = SSE[2]
    for i in range(start_res[1][1], end_res[1][1]):
        #current_res[1][1] = i 
        residue = structure[0][start_res[0]][(' ', i, ' ')]
        res_name = residue.get_resname()
        sequence.append(res_name)
    
    return sequence

def get_sequence(SSE, structure, residues):
    sequence = []
    start_res_id = SSE[1]
    end_res_id = SSE[2]
    #print(start_res_id,structure[0][start_res_id[0]][start_res_id[1]])
    start_res = structure[0][start_res_id[0]][start_res_id[1]]
    end_res = structure[0][end_res_id[0]][end_res_id[1]]
    stop = residues.index(end_res)
    start = residues.index(start_res)
    for i in range(start, stop+1):
        residue = residues[i].get_resname()
        sequence.append(residue)
    return sequence

def get_aa_len(SSE, structure, residues):
    sequence = []
    start_res_id = SSE[1]
    end_res_id = SSE[2]
    start_res = structure[0][start_res_id[0]][start_res_id[1]]
    end_res = structure[0][end_res_id[0]][end_res_id[1]]
    return residues.index(end_res)-residues.index(start_res)
    #return residues.index(end_res)-residues.index(start_res)

def create_nodes_and_adjacent_edges(SSE):
	nodes = []
	edges = []
	edge_lengths = []
	sse_lengths = []
	count = 0
	for i in SSE:
		tup = (count, {'SSE': i[0], 'start_res': i[1], 'end_res': i[2]})
		sse_lengths.append(tup[1]['end_res'][1][1] - tup[1]['start_res'][1][1])
		nodes.append(tup)
		if count != (len(SSE)-1):
			if tup[1]['end_res'][0] == SSE[count+1][1][0]:
				edges.append((count, count+1))
		count += 1
	edges.pop()
	for i in edges:
		sse1 = SSE[i[0]]
		sse2 = SSE[i[1]]
		residue1 = sse1[2]
		residue2 = sse2[1]
	
		if residue1[0] == residue2[0]:
			length = residue2[1][1] - residue1[1][1]
			edge_lengths.append(length)
		
	return nodes, edges, edge_lengths, sse_lengths

def create_graph(nodes, edges, structural_edges, structural_distances, sequential_lengths, sse_length, properties):
	#add nodes and edges to graph
	graph = nx.Graph()
	#graph.add_nodes_from(nodes)
	#graph.add_edges_from(edges, Type = 'sequential')
	#graph.add_edges_from(contact_edges, Type = 'structural')
	for n, node in enumerate(nodes):
		graph.add_node(node[0], structure=node[1]['SSE'], length=sse_length[n],
					hydrophobicity=properties[n][0], polarisability=properties[n][1], 
					vdw_volume=properties[n][2], polarity=properties[n][3])
		
	for n, edge in enumerate(edges):
		graph.add_edge(*edge, Type = 'sequential', length=sequential_lengths[n])
	
	for n, edge in enumerate(structural_edges):
		graph.add_edge(*edge, Type = 'structural', length=structural_distances[n])
	return graph
