import os
import pandas as pd
import modules.GO.geneOntology as geneOntology
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap
from matplotlib.backends.backend_pdf import PdfPages

inputDirectory = "./outputPDFs/DiffExpOut/Significant Genes/"

# For Loop

fileNamesDirectory = os.listdir(inputDirectory)
fileNames = []

for element in fileNamesDirectory:
    print(element)
    if ".csv" in element and "results" in element:
        fileNames.append(element)

print(fileNames)

for fileName in fileNames:
    fileName = fileName.split(".")[0]
    DESeq2Results = pd.read_csv(inputDirectory+fileName+".csv")
    DESeq2Results = DESeq2Results.sort_values(['log2FoldChange'], ascending = [False])

    positiveGenes = DESeq2Results[DESeq2Results['log2FoldChange'] > 0]
    negativeGenes = DESeq2Results[DESeq2Results['log2FoldChange'] < 0]
    positiveGenesList = list(positiveGenes["Unnamed: 0"])
    negativeGenesList = list(negativeGenes["Unnamed: 0"])

    if(len(positiveGenes) > 0):
        results = geneOntology.geneOntologyAnalysis(positiveGenesList)

        print(results)

        resultsDF = results
        resultsDF['per'] = resultsDF.n_genes/resultsDF.n_go
        resultsDF = resultsDF[0:10]

        print(resultsDF)
        if(len(resultsDF) > 0):
            fig, ax = plt.subplots(figsize = (0.5, 2.75))
            cmap = mpl.cm.bwr_r
            norm = mpl.colors.Normalize(vmin = resultsDF.p_corr.min(), vmax = resultsDF.p_corr.max())
            mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)
            cbl = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
            plt.figure(figsize = (10,10))
            ax = sns.barplot(data = resultsDF, x = 'per', y = 'term', palette = "mako")
            ax.set_yticklabels([textwrap.fill(e, 22) for e in resultsDF['term']])

            plt.savefig("./outputPDFs/DiffExpOut/GeneOntology/"+fileName+"Positive.pdf", bbox_inches='tight')

    if(len(negativeGenes) > 0):
        results = geneOntology.geneOntologyAnalysis(negativeGenesList)

        print(results)

        resultsDF = results
        resultsDF['per'] = resultsDF.n_genes/resultsDF.n_go
        resultsDF = resultsDF[0:10]

        print(resultsDF)
        if(len(resultsDF) > 0):
            fig, ax = plt.subplots(figsize = (0.5, 2.75))
            cmap = mpl.cm.bwr_r
            norm = mpl.colors.Normalize(vmin = resultsDF.p_corr.min(), vmax = resultsDF.p_corr.max())
            mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)
            cbl = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
            plt.figure(figsize = (10,10))
            ax = sns.barplot(data = resultsDF, x = 'per', y = 'term', palette = "mako")
            ax.set_yticklabels([textwrap.fill(e, 22) for e in resultsDF['term']])

            plt.savefig("./outputPDFs/DiffExpOut/GeneOntology/"+fileName+"Negative.pdf", bbox_inches='tight')
