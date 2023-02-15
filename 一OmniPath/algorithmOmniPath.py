import omnipath as op
import pandas as pd
import json

# Get All Interactions From OmniPath
output = op.interactions.AllInteractions()
allInteractions = output.get()

# Get Dict from Protein to Genes
proteinToGene = pd.read_csv('./proteinToGene.tsv', sep='\t')
proteinToGene = proteinToGene.set_index("To")
dictProteinToGene = proteinToGene.to_dict()["From"]

####################################################################
# Get Dict of Source and their Targets
####################################################################

totalDict = allInteractions.groupby('target').apply(lambda dfg: dfg.drop('target', axis=1).to_dict(orient='list')).to_dict()
jsonDict = json.dumps(totalDict, sort_keys=True, indent=4)

####################################################################
# Select Genes To Show (Algorithm)
####################################################################
RUNS = 1
selectedGenes = ["Q8WXX7"]
try:
    labelFile = dictProteinToGene["Q8WXX7"]
except:
    labelFile = ""

for run in range(RUNS):
    tempGenes = []
    for gene in selectedGenes:
        try:
            tempGenes += totalDict[gene]["source"]
            # print(f"Getting targets of {gene}")
        except Exception as e:
            RUNS = RUNS
            # print(f"Skipped {gene}")
        # print(tempGenes)
        # print("")
    selectedGenes += tempGenes
    # selectedGenes = list(dict.fromkeys(selectedGenes))
    print(len(selectedGenes))

print("\n Final Selected Genes:")
print(selectedGenes)

print(len(selectedGenes))
removeDuplicates = list(dict.fromkeys(selectedGenes))
print(len(removeDuplicates))

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
fallInteractions = allInteractions[allInteractions["target"].isin(removeDuplicates)]
print(fallInteractions.head(20))
fallInteractions = transformToGeneDF(fallInteractions)
print(fallInteractions.head(20))

fallInteractions.to_csv(f"./{labelFile}Proteins.csv", index=False)

####################################################################
# Label Cytoscape
####################################################################
import sys
import os
sys.path.insert(1, '../一Network Visualization/')
os.chdir('../一Network Visualization/')

import labelCytoscape

genesToLabelled = fallInteractions["source"].to_list() + fallInteractions["target"].to_list()
genesToLabelled = list(dict.fromkeys(genesToLabelled))

sys.path.insert(1, '../一OmniPath/')
os.chdir('../一OmniPath/')

labelCytoscape.createLabellingTable(genesToLabelled, f"./{labelFile}ProteinLabels.csv")