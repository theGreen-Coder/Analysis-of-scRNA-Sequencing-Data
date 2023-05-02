import networkx as nx
import omnipath as op
import pandas as pd
import json
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

PATH = "./geneNetworks/bHLH/"

# Get list of all HSA21 genes
HSA21genesDataframe = pd.read_csv("../一Network Visualization/data/HSA21_genes_biomaRt_conversion.csv")
HSA21genes = [x for x in HSA21genesDataframe["hgnc_symbol"] if str(x) != 'nan']

networks = [f for f in listdir(PATH) if isfile(join(PATH, f))]
newNetworks = []

for network in networks:
    if "Protein.csv" in network:
        newNetworks.append(network)
networks = newNetworks

networksDataframe = []

print(newNetworks)
for network in networks:
    print(f"Reading File: {network}")
    networksDataframe.append(pd.read_csv(PATH+network))

allNetworks = pd.concat(networksDataframe)
allNetworks = allNetworks.drop_duplicates(subset=['source', 'target'])
print(allNetworks)
allNetworks = allNetworks[~(allNetworks["target"].isin(HSA21genes))]
print(allNetworks)
allNetworks.to_csv(PATH+f"joinedNetworks.csv", index=False)

####################################################################
# Label Cytoscape
####################################################################
import sys
import os
sys.path.insert(1, '../一Network Visualization/')
os.chdir('../一Network Visualization/')

import labelCytoscape

genesToLabelled = allNetworks["source"].to_list() + allNetworks["target"].to_list()
genesToLabelled = list(dict.fromkeys(genesToLabelled))

sys.path.insert(1, '../一OmniPath/')
os.chdir('../一OmniPath/')

labelCytoscape.createLabellingTable(genesToLabelled, PATH+f"joinedNetworksLabels.csv")