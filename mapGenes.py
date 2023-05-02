import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

numberGenesChromosomes = {
    "1": 2061,
    "2": 1299,
    "3": 1081,
    "4": 757,
    "5": 882,
    "6": 1051,
    "7": 1010,
    "8": 701,
    "9": 778,
    "10": 730,
    "11": 1317,
    "12": 1037,
    "13": 322,
    "14": 821,
    "15": 617,
    "16": 863,
    "17": 1186,
    "18": 266,
    "19": 1476,
    "20": 546,
    "21": 221,
    "22": 495,
    "X": 859,
    "Y": 63,
}
geneLocations = {
    "1": 0, "2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0,"10": 0,"11": 0,"12": 0,
    "13": 0,"14": 0,"15": 0,"16": 0,"17": 0,"18": 0,"19": 0,"20": 0,"21": 0,"22": 0,"X": 0,"Y": 0
}
significantGenes = list(pd.read_csv("./一Expression Heatmaps/significantGenes.csv")["x"])
print(f"Len Significant Genes: {len(significantGenes)}")
listOfSignificantGenes = significantGenes # ["AUTS2", "OLIG2", "DSCAM"]
count = 0

for gene in listOfSignificantGenes:
    try:
        chromosomeNumber = str(data.genes_by_name(gene)[0].contig)
        if chromosomeNumber in geneLocations.keys():
            geneLocations[chromosomeNumber] += 1
    except:
        print("Unknown Gene!")
        count += 1

print(f"Total skipped genes: {count}")
print(geneLocations)
plt.figure(figsize = (5,4))
fig = sns.barplot(x=list(geneLocations.keys()),y=list(geneLocations.values()), color="black")
for item in fig.get_xticklabels():
    item.set_rotation(90)
plt.savefig('./Figures/Figure 2/Total Significant Genes.pdf', bbox_inches='tight')

relativeLocation = geneLocations

for chromo in relativeLocation:
    try:
        percentatge = (relativeLocation[chromo]/numberGenesChromosomes[chromo])*100
        relativeLocation[chromo] = percentatge
    except:
        error = 0

print(relativeLocation)
plt.figure(figsize = (5,4))
fig = sns.barplot(x=list(relativeLocation.keys()),y=list(relativeLocation.values()), color="black")
for item in fig.get_xticklabels():
    item.set_rotation(90)
plt.savefig('./Figures/Figure 2/Percentatge Significant Genes.pdf', bbox_inches='tight')