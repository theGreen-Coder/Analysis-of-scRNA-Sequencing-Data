import pandas as pd
import json
from tqdm.auto import tqdm

####################################################################
# Import network
####################################################################

networkDS = pd.read_csv('./data/adjDatasetFullGenes.tsv', sep='\t')

# Filter network at 50 importance
networkDS = networkDS[networkDS["importance"]>32]

####################################################################
# Import lists of genes
####################################################################

# Import HSA21 genes
HSA21genesDataframe = pd.read_csv("./data/HSA21_genes_biomaRt_conversion.csv")
HSA21genes = [x for x in HSA21genesDataframe["hgnc_symbol"] if str(x) != 'nan']

# Import more HSA21 genes
HSA21genesDataframeFull = pd.read_csv("./data/fullHSA21genes.csv")
HSA21genesFull = [x for x in HSA21genesDataframeFull["hgnc_symbol"] if str(x) != 'nan']

count = 0
for geneA in HSA21genes:
    if geneA in HSA21genesFull:
        count += 1
    else:
        HSA21genesFull.append(geneA)
HSA21genes = HSA21genesFull

# Import Transcription Factors
TFsList = open("./data/TFs.txt", "r").read().split("\n")

# Create list for TFs that also are HSA21 genes
HSA21_TFs = []
for gene in HSA21genes:
    if gene in TFsList:
        HSA21_TFs.append(gene)

HSA21Regulators = [
    "ZNF294", "LTN1", "RNF17", 
    "ZNF295", "ZBTB21", "KIAA1227",
    "Pred65", "ZNF355P", "PRED65", 
    "ZNF298", "PRDM15", 
    "APECED", 
    "KIAA0136", "MORC3", "ZCWCC3", "NXP2",
    "GCFC", "PAXBP1", "GCFC1",
    "SON", "NREBP", "BASS1",
    "PKNOX1", "PREP1", 
    "HSF2BP", "MEILB2", "POF19",
    "NRIP1", "RIP140", "NRIP1"
]

HSA21_TFs += HSA21Regulators

# Reselect which genes are on the TF list
newHSA21TFs = []
for gene in HSA21_TFs:
    if gene in TFsList:
        newHSA21TFs.append(gene)
HSA21_TFs = newHSA21TFs

# Import Neuro Development Disease genes
neuroDDdf = pd.read_csv('./data/NDDgenes.csv')
neuroDDdf = neuroDDdf[neuroDDdf["High Confidence NDD genes"] == True]
neuroDD = neuroDDdf["Symbol"].tolist()

# Import NeuroDevelopmental Genes
neuroDevGenesDF = pd.read_csv('./data/NervousSystemDevelopmentGO.tsv', sep='\t')
neuroDevGenes = neuroDevGenesDF["Gene"].tolist()


regulatedGenes = pd.read_csv("../一Expression Heatmaps/geneExpressionData/Log2FoldAllGenes.csv")
regulatedGenes.rename(columns={"Unnamed: 0": "Genes"}, inplace=True)
regulatedGenes = regulatedGenes.set_index('Genes')
upDownRegulated = regulatedGenes['Mean'].to_dict()

####################################################################
# Map list of genes on networkDS
####################################################################

# Type of lists
# HSA21genes TFsList HSA21_TFs neuroDD neuroDevGenes
listDict = {
    "GeneralTF": TFsList,
    "NeuroDev": neuroDevGenes,
    "NeuroDisease": neuroDD,
    "HSA21gene": HSA21genes,
    "HSA21TF": HSA21_TFs,
}

# Give type to TFs
listOfTFs = networkDS["TF"].tolist()
TF_Type = []
geneExpressionTF = []

for gene in listOfTFs:
    value = "None"
    for key in listDict:
        if gene in listDict[key]:
            value = key
    TF_Type.append(value)
    try:
        geneExpressionTF.append(upDownRegulated[gene]*-1)
    except:
        geneExpressionTF.append("Not Known")

networkDS["TF_Type"] =  TF_Type
networkDS["TF_Expression"] =  geneExpressionTF

# Give type to Targets
listOfTargets = networkDS["target"].tolist()
Target_Type = []
geneExpressionTarget = []

for gene in listOfTargets:
    value = "None"
    for key in listDict:
        if gene in listDict[key]:
            value = key
    Target_Type.append(value)
    try:
        geneExpressionTarget.append(upDownRegulated[gene]*-1)
    except:
        geneExpressionTarget.append("Not Known")

networkDS["Target_Type"] =  Target_Type
networkDS["Target_Expression"] =  geneExpressionTarget

print(networkDS)

####################################################################
# Get Dict of TFs and their Targets
####################################################################

dictionaryNetworkDataframe = networkDS.set_index('TF')
dictNet = dictionaryNetworkDataframe['target'].to_dict()

totalDict = networkDS.groupby('TF').apply(lambda dfg: dfg.drop('TF', axis=1).to_dict(orient='list')).to_dict()

jsonDict = json.dumps(totalDict, sort_keys=True, indent=4)

# print(totalDict["SOD1"])
# print(totalDict["SOD1"]["target"])
# print(totalDict["BACH1"]["target"])
# print("")

####################################################################
# Select Genes To Show (Algorithm)
####################################################################
RUNS = 3
selectedGenes = HSA21_TFs

for run in tqdm(range(RUNS)):
    tempGenes = []
    for gene in selectedGenes:
        try:
            tempGenes += totalDict[gene]["target"]
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

####################################################################
# Cytoscape List
####################################################################
fNetworkDS = networkDS[networkDS["TF"].isin(selectedGenes)]
print(fNetworkDS)

fNetworkDS.to_csv("./algorithmGenesi50.csv", index=False)

####################################################################
# Cytoscape Labelling Table
####################################################################
import labelCytoscape

genesToLabelled = fNetworkDS["TF"].to_list() + fNetworkDS["target"].to_list()
genesToLabelled = list(dict.fromkeys(genesToLabelled))

labelCytoscape.createLabellingTable(genesToLabelled, "./algorithmGenesLabelsi50.csv")