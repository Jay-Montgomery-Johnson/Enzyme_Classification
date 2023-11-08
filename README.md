# Enzyme_Classification
Part of my ongoing explorations of graph machine learning with relation to proteins. In particular I am classifying the protein structures of enzymes into one of 5 enzyme classes which define their function.

The scripts can be split into two main funtions. Preparing the dataset and training the model. I have constructed the dataset based on my attempt to replicate the result found in this paper: https://academic.oup.com/bioinformatics/article/21/suppl_1/i47/202991
which uses graph kenel methods to classify enzymes. 

This project takes protein structure files and finds secondary structure elements (SSE) such as alpha helices and beta pleated sheets. These SSEs are made into nodes on the graph and the bonds between nodes are edges. 

# Preparing the dataset
## download_pdb.py 
The PDB is the centalrepository where protein crystal structures are stored. I used the API request constructor to write an api request to search for only proteins I requied. I used the requests library in python to download a dataset of protein files.

## fix_xml_flaw.py
Some of the XML files were formatted incorreclt preventing python from reading them. So I adjusted the XML files with this script.

## Graph_Utils.py
This file contains a collections of functions which are used to process a pdb protein structure file into a networkX network.

## protein_properties.py
Amino acids are the building blocks of proteins and each have different chemical and physical properties. The amino acid index is a dataset of published results on amino acid physical properties. This file contains a class which assigns properties to polymer chains of amino acids. 

## download_dataset.py
This where pdb files are transformed into networkX networks and the nodes and edges are assigned properies with using protein_properties.py

# The Model
load_data.py, model.py and train.py are an implementation of a graph convolutional neural network. The model is taken from the pyTorch Geometric website. 
This model reached an accuracy of 40% on the test dataset so it is sucessfully classifying the enzymes but it is far from the 75% accuracy rates I have seen in the literature. Next I want to implement a graph kernal to use a SVM and after that I will tinker with graphs as it is likely that engineering some features can improve my results. 
