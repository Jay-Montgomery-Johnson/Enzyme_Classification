'''
This file turns a protin structure in the form of a PDB file into a 
NetworkX graph object. The nodes are secondary structre elements (SSE) 
and the edges are bonds.

SSE are defined by DSSP
'''
import networkx as nx
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import Bio.PDB.Residue as residue
import Bio.PDB as pdb
import matplotlib.pyplot as plt
import Bio.PDB.Selection as Selection
import numpy as np
from protein_properties import aa_index


def get_SSE(chain, residue, SSE): #this is not a good solution I should have used a 
    # given the list of SSE and a particular residue return the index of 
    #which SSE you are in
    res_id = residue.get_id()
    current_SSE = None
    for i, n  in enumerate(SSE):
        #print(chain.get_id(), n[1][0], res_id[1], n[1][1][1], n[1][1][1])
        #print(res_id[1], n[1][1][1], n[2][1][1])
        if chain.get_id() == n[1][0] and res_id[1] >= n[1][1][1] and res_id[1] <= n[2][1][1]:
            current_SSE = i
    if current_SSE == None:
        return 'error'
    else:
        return current_SSE

def filter_contacts(contacts, chain, res_id, SSE, current_SSE):
	#remove contacts from a contact list if they are in the same SSE
	#or if the residue id a water...
	new_contacts = []
	for i in contacts:
		contact_SSE = get_SSE(chain, i, SSE)
		if current_SSE != contact_SSE and i.get_resname() != 'HOH' and res_id[1] != i.get_id()[1]+1 and res_id[1] != i.get_id()[1]-1:
			#remove contacts which are from adjacent residues
			new_contacts.append(i)
	return new_contacts

def add_SSE_to_dict(contact_dict, new_contacts_SSE, current_SSE):
	if len(new_contacts_SSE) != 0:
		if current_SSE in contact_dict:
			contact_dict[current_SSE] += new_contacts_SSE
		else:
			contact_dict[current_SSE] = new_contacts_SSE
	return contact_dict

def distance_between_residues(residue1, residue2):
	#determine the distance between the Ca atoms of two residues
	print(residue1['CA'] - residue2['CA'])

def residue_midpoint(residue1, residue2):
	coord1 = residue1['CA'].get_coord()
	coord2 = residue2['CA'].get_coord()
	
	midpoint = (coord1 + coord2)/2
	
	return midpoint

def generate_contact_map(structure, SSE):
	#Determines a list of contacts for a given structure
	#iterate through each SSE
	
	atoms = Selection.unfold_entities(structure, 'A')
	neighbours = pdb.NeighborSearch(atoms)
	
	contact_dict = {}
	for chain in structure[0]:
		for residue in structure[0][chain.get_id()]:
			res_id = residue.get_id()
			#enumerate secondary structure elements
			#determine which SSE the current residue is in and add to overall dict
			if residue.get_resname() != 'HOH':
				current_SSE = get_SSE(chain, residue, SSE)
				contacts = neighbours.search(structure[0][chain.get_id()][residue.get_id()[0], residue.get_id()[1], residue.get_id()[2]]['CA'].coord, 15, level='R')    
					
				#remove contacts which are in the same SSE also remove contacts with water
				new_contacts = filter_contacts(contacts, chain, res_id, SSE, current_SSE)
				new_contacts_SSE = []
				for i in new_contacts:
					i_code, i_model, i_chain, i_residue = i.get_full_id()
					b = get_SSE(structure[i_model][i_chain], i, SSE)
					if b != 'error': #the errors are not helpful
						new_contacts_SSE.append(b)
				contact_dict = add_SSE_to_dict(contact_dict, new_contacts_SSE, current_SSE)
	
	#for i in contact_dict:
	for key in contact_dict:
		contact_dict[key] = set(contact_dict[key])
	return contact_dict

def load_PDB(pdb_id):
	parser = PDBParser(PERMISSIVE=0)
	structure = parser.get_structure(pdb_id, pdb_id + '.pdb')
	return structure
	
def load_PDB_dir(pdb_id, path):
	parser = PDBParser(PERMISSIVE=0)
	structure = parser.get_structure(pdb_id, path + pdb_id + '.pdb')
	return structure

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
	
def create_contact_edges(structure, SSE):
	#get contacts
	contact_dict = generate_contact_map(structure, SSE)
	contact_edges = []
	for k in contact_dict:
		for i in contact_dict[k]:
			e = [k, i]
			contact_edges.append(e)
	return contact_edges

def remove_SSE_with_aa_length_less_than_or_equal_to_n(SSE, n):
	new_SSE = []
	for i in SSE:
		aa_length = i[2][1][1]-i[1][1][1]
		if aa_length > n:
			new_SSE.append(i)
	return new_SSE

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

def create_structural_edges(structure, SSE):
	SSE_midpoints = []
	for i in SSE:
		residue1 = structure[0][i[1][0]][i[1][1]]
		residue2 = structure[0][i[2][0]][i[2][1]]
		midpoint = residue_midpoint(residue1, residue2)
		SSE_midpoints.append(midpoint)
	
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

def get_sequence(SSE, structure):
	#here SSE refers to an element of the SSE list
	sequence = []
	start_res = SSE[1]
	end_res = SSE[2]
	current_res = start_res
	for i in range(start_res[1][1], end_res[1][1]):
		#current_res[1][1] = i 
		residue = structure[0][current_res[0]][i]
		res_name = residue.get_resname()
		sequence.append(res_name)
	
	return sequence

