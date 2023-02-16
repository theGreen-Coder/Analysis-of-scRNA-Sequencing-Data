import networkx as nx
import omnipath as op
import pandas as pd
import json
import matplotlib.pyplot as plt

regGenes = ["SOX4", "DLX5", "DLX4", "DLX6", "DLX2", "RBFOX2", "CELF4", "CELF5", "EBF1", "SIX3", "PBX3", "HMGB2"]
regProts = ["Q06945", "P56178", "Q92988", "P56179", "Q07687", "O43251", "Q9BZC1", "Q8N6W0", "Q9UH73", "O95343", "P40426", "P26583"]

FILTER_BY_SOURCES = False
FILTER_NUM_SOURCES = 1
FILTER_BY_SHORTEST_PATH_MAX = False
FILTER_INHIBITION_STIMULATION = False
FILTER_SHORTEST_PATH_MAX = 5
DEV_GEN = regGenes[0]
DEV_GEN_NAME = regProts[0]
OUTPUT_NAME = DEV_GEN_NAME + f"_filterSOURCE{str(FILTER_BY_SOURCES)+str(FILTER_NUM_SOURCES)}" + f"_filterSHORTPATH{str(FILTER_BY_SHORTEST_PATH_MAX)+str(FILTER_SHORTEST_PATH_MAX)}"

# Get list of all HSA21 genes
HSA21genesDataframe = pd.read_csv("../一Network Visualization/data/HSA21_genes_biomaRt_conversion.csv")
HSA21genes = [x for x in HSA21genesDataframe["hgnc_symbol"] if str(x) != 'nan']

# Get All Interactions From OmniPath
output = op.interactions.AllInteractions()
allInteractions = output.get()

# Get Dict from Protein to Genes
proteinToGene = pd.read_csv('./proteinToGene.tsv', sep='\t')
proteinToGene = proteinToGene.set_index("To")
dictProteinToGene = proteinToGene.to_dict()["From"]

# Get Dict from Protein to Genes
geneToProtein = pd.read_csv('./proteinToGene.tsv', sep='\t')
dictGeneToProtein = geneToProtein.groupby('From').apply(lambda dfg: dfg.drop('From', axis=1).to_dict(orient='list')).to_dict()

####################################################################
# Convert HSA21 Genes to Protein Accession Numbers
####################################################################

# Convert geneNames to Protein Accession Numbers
proteinNames = []
count = 0

for gene in HSA21genes:
    try:
        proteinNames += dictGeneToProtein[gene]["To"]
    except:
        count +=1

print(f"Total genes skipped: {count}")

proteinNames = list(dict.fromkeys(proteinNames))
print(proteinNames)

####################################################################
# Find All Shortest Paths
####################################################################

if FILTER_BY_SOURCES:
    allInteractions = allInteractions[allInteractions["n_sources"] >= FILTER_NUM_SOURCES]

if FILTER_INHIBITION_STIMULATION:
    allInteractions = allInteractions[~((allInteractions["is_inhibition"] == False) & (allInteractions["is_stimulation"] == False))]

G = nx.from_pandas_edgelist(allInteractions, source='source', target='target', create_using=nx.DiGraph())

allShortestPathGenes = []
count = 0
##### ------ Implement All Protein Accession numbers from Gene ------
for protein in proteinNames:
    try:
        tempShortPath = []
        if FILTER_BY_SHORTEST_PATH_MAX:
            tempShortPath = []
            tempShortPaths = [p for p in nx.all_shortest_paths(G, source=protein, target=DEV_GEN)]
            for path in tempShortPaths:
                if len(path) <= FILTER_SHORTEST_PATH_MAX:
                    tempShortPath += path
            tempShortPath = list(dict.fromkeys(tempShortPath)) 
        else:
            tempShortPath = list(dict.fromkeys(sum([p for p in nx.all_shortest_paths(G, source=protein, target=DEV_GEN)], []))) 
        allShortestPathGenes += tempShortPath
    except:
        count += 1

allShortestPathGenes = list(dict.fromkeys(allShortestPathGenes))
print("ALL SHORTEST PAHTS GENES:")
print(f"Total Skipped Genes: {count}")
print(allShortestPathGenes)

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

####################################################################
# Cytoscape List
####################################################################
filterAllInteractions = allInteractions[
    (allInteractions["source"].isin(allShortestPathGenes)) &
    (allInteractions["target"].isin(allShortestPathGenes))]
print(filterAllInteractions)
filterAllInteractions = transformToGeneDF(filterAllInteractions)
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