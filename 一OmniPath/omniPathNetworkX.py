import networkx as nx
import requests
import omnipath as op
import pandas as pd
import json
import matplotlib.pyplot as plt

regulatorGenes = ["SOX4", "DLX5", "DLX4", "DLX6", "DLX2", "RBFOX2", "CELF4", "CELF5", "EBF1", "SIX3", "PBX3", "HMGB2"]
regProts = ["Q06945", "P56178", "Q92988", "P56179", "Q07687", "O43251", "Q9BZC1", "Q8N6W0", "Q9UH73", "O95343", "P40426", "P26583"]
bHLHGenes = ["ID1", "ID2", "ID3", "ID4", "HES1", "HES2", "HES3", "HES5", "HES6", "HEY1", "HEY2", "TCF3", "TCF4", "TCF12", "ATOH7", "ATOH1", 
             "NEUROD1", "NEUROD2", "NEUROD4", "NEUROD6", "NEUROG1", "NEUROG2", "NEUROG3", "OLIG1", "OLIG2", "OLIG3", "BHLHE22", "BHLHE23",
             "ASCL1", "ASCL2", "ASCL3", "ASCL4", "NHLH1", "NHLH2"]

FILTER_BY_SOURCES = False
FILTER_NUM_SOURCES = 2
FILTER_BY_SHORTEST_PATH_MAX = True
FILTER_INHIBITION_STIMULATION = False
FILTER_SHORTEST_PATH_MAX = 3
DEV_GEN = "NHLH2"
OUTPUT_NAME = DEV_GEN + f"_filterSOURCE{str(FILTER_BY_SOURCES)+str(FILTER_NUM_SOURCES)}" + f"_filterSHORTPATH{str(FILTER_BY_SHORTEST_PATH_MAX)+str(FILTER_SHORTEST_PATH_MAX)}"

####################################################################
# Transform Protein Dataframe to Gene Dataframe
####################################################################
def transformToGeneDF(proteinDF):
    sourceNames = list(proteinDF["source"])
    sourceNewNames = []

    targetNames = list(proteinDF["target"])
    targetNewNames = []

    for name in sourceNames:
        try:
            sourceNewNames.append(dictProteinToGene[str(name)])
        except:
            sourceNewNames.append(name)
    
    for name in targetNames:
        try:
            targetNewNames.append(dictProteinToGene[str(name)])
        except:
            targetNewNames.append(name)

    proteinDF["source"] = sourceNewNames
    proteinDF["target"] = targetNewNames
    return proteinDF

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'    

print(bcolors.FAIL + "Importing Data..." + bcolors.ENDC)
# Get list of all HSA21 genes
HSA21genesDataframe = pd.read_csv("../一Network Visualization/data/HSA21_genes_biomaRt_conversion.csv")
HSA21genes = [x for x in HSA21genesDataframe["hgnc_symbol"] if str(x) != 'nan']

# Get Dict from Protein to Genes
proteinToGene = pd.read_csv('./proteinToGene.tsv', sep='\t')
proteinToGene = proteinToGene.set_index("To")
dictProteinToGene = proteinToGene.to_dict()["From"]

# Get Dict from Protein to Genes
geneToProtein = pd.read_csv('./proteinToGene.tsv', sep='\t')
dictGeneToProtein = geneToProtein.groupby('From').apply(lambda dfg: dfg.drop('From', axis=1).to_dict(orient='list')).to_dict()

# Get All Interactions From OmniPath
output = op.interactions.AllInteractions()
allInteractions = output.get()
allInteractions = transformToGeneDF(allInteractions)

####################################################################
# Find All Shortest Paths
####################################################################
proteinNames = HSA21genes

if FILTER_BY_SOURCES:
    print(bcolors.OKBLUE + f"Flitering by N_SOURCES >= {FILTER_NUM_SOURCES}..." + bcolors.ENDC)
    allInteractions = allInteractions[allInteractions["n_sources"] >= FILTER_NUM_SOURCES]

if FILTER_INHIBITION_STIMULATION:
    print(bcolors.OKBLUE + f"Flitering non-INH/STIM undefined interactions..." + bcolors.ENDC)
    allInteractions = allInteractions[~((allInteractions["is_inhibition"] == False) & (allInteractions["is_stimulation"] == False))]

print(bcolors.FAIL + "Finding All Shortest Paths..." + bcolors.ENDC)

G = nx.from_pandas_edgelist(allInteractions, source='source', target='target', create_using=nx.DiGraph())

allShortestPathGenes = []
count = 0
for protein in proteinNames:
    try:
        tempShortPath = []
        if FILTER_BY_SHORTEST_PATH_MAX:
            tempShortPath = []
            tempShortPaths = [p for p in nx.all_shortest_paths(G, source=protein, target=DEV_GEN)]
            print(bcolors.OKCYAN + f"Nº of Paths: {str(len(tempShortPaths))}" + bcolors.ENDC)
            for path in tempShortPaths:
                if len(path) <= FILTER_SHORTEST_PATH_MAX:
                    tempShortPath += path
            tempShortPath = list(dict.fromkeys(tempShortPath)) 
        else:
            tempShortPath = list(dict.fromkeys(sum([p for p in nx.all_shortest_paths(G, source=protein, target=DEV_GEN)], []))) 
        allShortestPathGenes += tempShortPath
        print(bcolors.BOLD + f"Gene Path Found: {str(tempShortPath)}" + bcolors.ENDC)

    except Exception as e:
        print(e)
        count += 1

allShortestPathGenes = list(dict.fromkeys(allShortestPathGenes))
print(bcolors.FAIL + f"Total Found Genes: {len(allShortestPathGenes)}" + bcolors.ENDC)
print(bcolors.FAIL + f"Total Error Genes: {count}" + bcolors.ENDC)

####################################################################
# Cytoscape List
####################################################################
print(bcolors.FAIL + "Cytoscape Interactions Result..." + bcolors.ENDC)
filterAllInteractions = allInteractions[
    (allInteractions["source"].isin(allShortestPathGenes)) &
    (allInteractions["target"].isin(allShortestPathGenes))]
print(filterAllInteractions)

filterAllInteractions.to_csv(f"./geneNetworks/{OUTPUT_NAME}_Protein.csv", index=False)

####################################################################
# Label Cytoscape
####################################################################
import sys
import os
sys.path.insert(1, '../一Network Visualization/')
os.chdir('../一Network Visualization/')

import labelCytoscape

genesToLabelled = filterAllInteractions["source"].to_list() + filterAllInteractions["target"].to_list()
genesToLabelled = list(dict.fromkeys(genesToLabelled))

sys.path.insert(1, '../一OmniPath/')
os.chdir('../一OmniPath/')

labelCytoscape.createLabellingTable(genesToLabelled, f"./geneNetworks/{OUTPUT_NAME}_ProteinLabels.csv")