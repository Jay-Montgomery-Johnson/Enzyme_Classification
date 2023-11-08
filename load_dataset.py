import networkx as nx
import os, glob
import numpy as np
import torch
from torch_geometric.data import Data
from torch_geometric.utils.convert import from_networkx
from torch_geometric.loader import DataLoader

def load_data():
    structure_element = {'-':0., 'H':1., 'E':2., 'T':3., 'S':4., 'G':5., 'I':6., 'B':7.}
    y_elem = {1: torch.tensor([1,0,0,0,0,0]), 2: torch.tensor([0,1,0,0,0,0]), 
            3:torch.tensor([0,0,1,0,0,0]) , 4:torch.tensor([0,0,0,1,0,0]), 
            5:torch.tensor([0,0,0,0,1,0]), 6:torch.tensor([0,0,0,0,0,1])}
    count = 0
    graphs = []
    for ec in range(1,7):
        path = f'dataset/ec_{ec}/graph/'
        print(ec)
        ec_graph = []
        for filename in glob.glob(os.path.join(path,'*.gexf')):
            try:
                graph = nx.read_gexf(filename)
                #change the string elements to integers
                for node in graph.nodes(data=True):
                    node[1]['structure'] = structure_element[node[1]['structure']]
                    del node[1]['label']
                for edge in graph.edges(data=True):
                    if edge[2]['Type'] == 'structural':
                        edge[2]['Type'] = 1.
                    else:
                        edge[2]['Type'] = 0.
                    del edge[2]['id']
                graph = nx.convert_node_labels_to_integers(graph)
                ec_graph.append(graph)
            except:
                print('error')
        graphs.append(ec_graph)
    
    tensors = []
    for num, ec in enumerate(graphs):
        for g in ec:
            n = list(g.nodes(data=True))
            e = list(g.edges(data=True))
            ng = nx.Graph()
            ng.add_nodes_from(n)
            ng.add_edges_from(e)
            datum = from_networkx(ng)
            datum.y = y_elem[num + 1]
            tensors.append(datum)
            
    #return a dataloader
    for i, datum in enumerate(tensors):
        x = torch.stack((datum.structure, datum.length, datum.hydrophobicity, 
                     datum.polarisability, datum.vdw_volume, datum.polarity),dim=1)
        #print(x)
        means = x.mean(dim=0, keepdim=True)
        stds = x.std(dim=1, keepdim=True)
        x = (x - means) / stds
        #print(x)
        datum.x = x
        
    np.random.shuffle(tensors)
    train = tensors[:450]
    test = tensors[450:]
    train_loader = DataLoader(train, batch_size=64, shuffle=True)
    test_loader = DataLoader(test, batch_size=64, shuffle=False)
    return train_loader, test_loader
        

