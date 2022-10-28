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

resultsFile = './dataSave/savedData.h5ad'
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
    if word != "2clustering.py.py" and word != "-showPlots":
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
    adata = sc.read("./output/1QualityData.h5ad")
adata.var_names_make_unique()

geneMarkers = pd.read_csv('./dataSaveOriginal/cellTypeMarkers.csv')
geneNames = geneMarkers["gene"].tolist()

####################################################################################################
# Clustering 
####################################################################################################
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.9)


####################################################################################################
# Clustering Plots
####################################################################################################

print("Saving Clustering Plots")
with PdfPages(outputDirectory+'Clustering Plots '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.pca(adata, show=showPlots)
    sc.pl.pca_loadings(adata, show=showPlots)
    sc.pl.pca_variance_ratio(adata, show=showPlots)
    sc.pl.umap(adata, color='leiden', show=showPlots)
    sc.pl.umap(adata, color='sample', show=showPlots)
    sc.pl.umap(adata, color='group', groups="DS", show=showPlots)
    sc.pl.umap(adata, color='group', show=showPlots)
    sc.pl.umap(adata, color='n_genes_by_counts', show=showPlots)
    sc.pl.umap(adata, color='pct_counts_mt', show=showPlots)
    sc.pl.pca(adata, color='sample', components='1,2', show=showPlots)
    sc.pl.pca(adata, color='group', components='1,2', show=showPlots)
    sc.pl.pca(adata, color='sample', components='3,4', show=showPlots)
    sc.pl.pca(adata, color='group', components='3,4', show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_DS2U", show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_H9", show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_IMR", show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_ihtc", show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_2DS3", show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_DS1", show=showPlots)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_DSP", show=showPlots)
    sc.pl.umap(adata, color=geneNames, show=showPlots)
    sc.pl.dotplot(adata, geneNames, groupby='leiden', show=showPlots)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP
