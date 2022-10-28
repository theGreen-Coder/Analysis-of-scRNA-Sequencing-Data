from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import os
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


####################################################################################################
# Global Settings 
####################################################################################################

np.set_printoptions(threshold=sys.maxsize)
pd.options.display.max_columns = sys.maxsize
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

resultsFile = './output/1QualityData.h5ad'
outputDirectory = "./outputPDFs/"

# Arguments Settings

sysArgs = sys.argv

showPlots = False
if "-showPlots" in sysArgs:
    showPlots = True
else:
    showPlots = False

arg = ""
for word in sysArgs:
    if word != "1qualityControl.py" and word != "-showPlots":
        arg = arg + word + " "
        if "-" not in word:
            OUTPUT_DIR = outputDirectory+word+"/"
            CHECK_FOLDER = os.path.isdir(OUTPUT_DIR)
            if not CHECK_FOLDER:
                os.makedirs(OUTPUT_DIR)
            outputDirectory = OUTPUT_DIR
            

####################################################################################################
# Import Data 
####################################################################################################

if "-loadGenomics" in sysArgs:
    adata = loadScanpy.readData()
else:
    adata = sc.read("./dataSaveOriginal/rawDataset5000.h5ad")
adata.var_names_make_unique()

geneMarkers = pd.read_csv('./dataSaveOriginal/cellTypeMarkers.csv')
geneNames = geneMarkers["gene"].tolist()


####################################################################################################
# QC Calculations 
####################################################################################################

# Percentatge of Mitochondrial Genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Get Top 50 Genes
percentageList = []
print("Calculating Percentatge of Top 50 Genes")
for row in adata.X:
    firstRow = np.sort(row.toarray()[0])
    sum50 = firstRow[-50:]
    percentatge = (sum(sum50)/sum(firstRow))*100
    percentageList.append(percentatge)
adata.obs["percentageTop50"] = percentageList

####################################################################################################
# QC Filtering 
####################################################################################################

if "-filter" in sysArgs:
    print("Filtering Damaged Cells")
    adata = adata[adata.obs.n_genes_by_counts > 1000, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
else:
    print("Not Filtering Damaged Cells")

####################################################################################################
# QC Plots 
####################################################################################################
print("Saving QC Plots")
with PdfPages(outputDirectory+'Quality Control Plots '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    violinPlots = sc.pl.violin(adata,['nCount_RNA', 'nFeature_RNA', 'pct_counts_mt', "percentageTop50"], groupby="sample", jitter=0.4, multi_panel=True, show=showPlots)
    scatterPlots = sc.pl.scatter(adata, x='pct_counts_mt', y='percentageTop50', show=showPlots)
    scatterPlots = sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=showPlots)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

adata.write(resultsFile)
