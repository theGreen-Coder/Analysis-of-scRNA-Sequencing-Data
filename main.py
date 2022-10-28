from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import loadScanpy as loadScanpy
import classifyClusters as classify
import os
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


####################################################################################################
# Global Settings 
####################################################################################################
np.set_printoptions(threshold=sys.maxsize)
pd.options.display.max_columns = sys.maxsize
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

results_file = './dataSave/savedData.h5ad'
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
    if word != "main.py" and word != "-showPlots":
        arg = arg + word + " "
        # if "-" not in word:
        #     OUTPUT_DIR = outputDirectory+word+"/"
        #     CHECK_FOLDER = os.path.isdir(OUTPUT_DIR)
        #     if not CHECK_FOLDER:
        #         os.makedirs(OUTPUT_DIR)
        #     outputDirectory = OUTPUT_DIR
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

####################################################################################################
# Import Data 
####################################################################################################
print(bcolors.FAIL + "Importing Data..." + bcolors.ENDC)
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
print(bcolors.FAIL + "Performing Quality Control..." + bcolors.ENDC)
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
print(bcolors.FAIL + "Applying Quality Control Filters..." + bcolors.ENDC)
if "-filter" in sysArgs:
    print("Filtering Damaged Cells")
    adata = adata[adata.obs.n_genes_by_counts > 1000, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]
else:
    print("Not Filtering Damaged Cells")

####################################################################################################
# QC Plots 
####################################################################################################
print(bcolors.FAIL + "Saving Quality Control Plots..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Quality Control Plots '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    violinPlots = sc.pl.violin(adata,['nCount_RNA', 'nFeature_RNA', 'pct_counts_mt', "percentageTop50"], groupby="sample", jitter=0.4, multi_panel=True, show=showPlots)
    scatterPlots = sc.pl.scatter(adata, x='pct_counts_mt', y='percentageTop50', show=showPlots)
    scatterPlots = sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=showPlots)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

####################################################################################################
# Clustering 
####################################################################################################
print(bcolors.FAIL + "Performing Clustering..." + bcolors.ENDC)
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
sc.tl.leiden(adata, resolution=0.4)

####################################################################################################
# Clustering Plots
####################################################################################################
print(bcolors.FAIL + "Saving Clustering Plots..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Clustering Plots '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.pca(adata, show=showPlots)
    sc.pl.pca_loadings(adata, show=showPlots)
    sc.pl.pca_variance_ratio(adata, show=showPlots)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
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
    sc.pl.umap(adata, color=["DLX2", "SOX2"], show=showPlots)
    sc.pl.umap(adata, color=["MKI67", "EOMES", "DLX2", "GLI3", "NEUROD6", "AQP4", "MEF2C"], title=["Marker Gene Expression"], show=showPlots) #Marker Gene Expression
    sc.pl.umap(adata, color=["PGK1", "ALDOA", "ARCN1", "GORASP2"], title=["Gene Stress Markers"], show=showPlots) #Marker Stress Genes
    sc.pl.umap(adata, color=["APP", "DYRK1A"], groups="CON", show=showPlots)
    sc.pl.umap(adata, groups="DS", color=["APP", "DYRK1A"],  show=showPlots)
    sc.pl.umap(adata, groups="CON", color=["APP", "DYRK1A"], show=showPlots)
    sc.pl.dotplot(adata, geneNames, groupby='leiden', show=showPlots)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP
# DS CON

####################################################################################################
# Finding Marker Genes
####################################################################################################
print(bcolors.FAIL + "Finding Marker Genes..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Marker Genes '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False, show=showPlots)

    # sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=showPlots)
    # sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
    # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=showPlots)

    getColumns = classify.findCellType(["VIM"])
    dfObj = pd.DataFrame(columns = getColumns["Cell Type"].to_list())

    topGenesPerCluster = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
    newClusterNames = []
    rowDataset = 0
    for column in topGenesPerCluster:
        genesList = topGenesPerCluster[column].to_list()
        genesResult = classify.findCellType(genesList)
        
        dfObj.loc[rowDataset] = genesResult["Total Found Genes"].to_list()
        rowDataset +=1

        genesResult = genesResult.sort_values(['Total Found Genes', 'Total Average Difference Genes'], ascending = [False, False])
        isInList = True
        i = 0
        while isInList:
            if genesResult.iloc[i,0] in newClusterNames:
                i += 1
            else:
                newClusterNames.append(genesResult.iloc[i,0])
                isInList = False
                print("Cluster : "+genesResult.iloc[i,0])
        
        # isInList = True
        # i = 0
        # while isInList:
        #     if genesResult.iloc[0,0] in newClusterNames:
        #         i += 1
        #         genesResult.iloc[0,0] += str(i) 
        #     else:
        #         newClusterNames.append(genesResult.iloc[0,0])
        #         isInList = False
        #         print("Cluster : "+genesResult.iloc[0,0])
        
    adata.rename_categories('leiden', newClusterNames)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')

    

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')


N=200
THE_FIGURE = plt.figure(figsize=(8.27, 6), dpi=300)
ax = plt.subplot(1, 1, 1)
sns.heatmap(dfObj, annot=True)
THE_FIGURE.savefig('image.pdf', bbox_inches='tight', pad_inches=0.1)

