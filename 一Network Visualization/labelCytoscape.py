import pandas as pd

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
    "ZNF294", "LTN1", "RNF160", 
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

def createLabellingTable(geneNamesList, fileName):
    labellingTable = pd.DataFrame(columns=['name', 'shared name','type','avgExpression'])
    labellingTable["name"] = geneNamesList
    labellingTable["shared name"] = geneNamesList
    # Give type to TFs
    listOfGenes = labellingTable["name"].tolist()
    Gene_Type = []
    geneExpression = []

    for gene in listOfGenes:
        value = "None"
        for key in listDict:
            if gene in listDict[key]:
                value = key
        Gene_Type.append(value)
        try:
            geneExpression.append(upDownRegulated[gene]*-1)
        except:
            geneExpression.append(0)
    labellingTable["type"] =  Gene_Type
    labellingTable["avgExpression"] =  geneExpression
    
    labellingTable.to_csv(fileName, index=False)