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

inputData = "./output/savedData.h5ad" # ./dataSaveOriginal/rawDataset5000.h5ad
results_file = './output/savedData.h5ad'
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
    adata = sc.read("./dataSaveOriginal/rawDataset.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

print(adata)

geneMarkers = pd.read_csv('./dataSaveOriginal/cellTypeMarkers.csv')
geneNames = geneMarkers["gene"].tolist()

adataCON_DS2U =  adata[adata.obs["sample"] == "CON_DS2U"]
adataCON_H9 =  adata[adata.obs["sample"] == "CON_H9"]
adataCON_IMR =  adata[adata.obs["sample"] == "CON_IMR"]
adataCON_ihtc =  adata[adata.obs["sample"] == "CON_ihtc"]

adataDS_2DS3 =  adata[adata.obs["sample"] == "DS_2DS3"]
adataDS_DS1 =  adata[adata.obs["sample"] == "DS_DS1"]
adataDS_DSP =  adata[adata.obs["sample"] == "DS_DSP"]

adata = adataCON_DS2U.concatenate(adataCON_H9)
adata = adata.concatenate(adataCON_IMR)
# adata = adata.concatenate(adataCON_ihtc)
adata = adata.concatenate(adataDS_2DS3)
adata = adata.concatenate(adataDS_DSP)
# adata = adata.concatenate(adataDS_DS1)

# print(adata)
# adata.write(results_file)

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP

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
    plt.suptitle("Violin Plots QC Metrics", y=1, fontsize=25)
    scatterPlots = sc.pl.scatter(adata, x='pct_counts_mt', y='percentageTop50', show=showPlots)
    plt.suptitle("% MT vs % Top 50 Genes", y=1, fontsize=25)
    scatterPlots = sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=showPlots)
    plt.suptitle("Total Counts vs Nº of Genes by Count", y=1, fontsize=25)

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
sc.tl.leiden(adata, resolution=2)
# adata.write(results_file)

####################################################################################################
# Integration
####################################################################################################
print(bcolors.FAIL + "Integrating Datasets..." + bcolors.ENDC)
if "-noIntegration" not in sysArgs and "-ingestIntegration" not in sysArgs:
    print("Performing BBKNN Integration...")
    with PdfPages(outputDirectory+'Integration Plots '+arg+'.pdf') as pdf:
        dummyPlot = plt.plot([1, 2])
        figureNum = plt.gcf().number

        sc.pl.umap(adata, color=['sample', "group"], show=showPlots)
        plt.suptitle("Non-Integrated Dataset UMAP", y=1, fontsize=25)
        sc.external.pp.bbknn(adata, batch_key='sample')
        sc.tl.umap(adata)
        plt.suptitle("Integrated Dataset UMAP", y=1, fontsize=25)
        sc.pl.umap(adata, color='sample', show=showPlots)
        plt.suptitle("Integrated Dataset by Sample", y=1, fontsize=25)
        sc.pl.umap(adata, color='group', show=showPlots)
        plt.suptitle("Integrated Dataset by Group", y=1, fontsize=25)

        for fig in range(figureNum+1,  plt.gcf().number+1):
            pdf.savefig(figure=fig, bbox_inches='tight')
elif "-ingestIntegration" in sysArgs:
    print("Performing Ingest Integration...")
    with PdfPages(outputDirectory+'Integration Plots '+arg+'.pdf') as pdf:
        dummyPlot = plt.plot([1, 2])
        figureNum = plt.gcf().number

        adataCON_DS2U =  adata[adata.obs["sample"] == "CON_DS2U"]
        print(adataCON_DS2U)
        adataCON_H9 =  adata[adata.obs["sample"] == "CON_H9"]
        print(adataCON_H9)
        adataCON_IMR =  adata[adata.obs["sample"] == "CON_IMR"]
        print(adataCON_IMR)
        adataDS_2DS3 =  adata[adata.obs["sample"] == "DS_2DS3"]
        print(adataDS_2DS3)
        adataDS_DSP =  adata[adata.obs["sample"] == "DS_DSP"]
        print(adataDS_DSP)

        sc.pl.umap(adata, color=['sample', "group"], show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="CON_DS2U", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="CON_H9", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="CON_IMR", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="DS_2DS3", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="DS_DSP", show=showPlots)

        sc.tl.pca(adataCON_IMR, svd_solver='arpack')
        sc.pp.neighbors(adataCON_IMR, n_neighbors=10, n_pcs=40)

        sc.tl.ingest(adataCON_DS2U, adataCON_IMR, obs="leiden")
        sc.tl.ingest(adataCON_H9, adataCON_IMR, obs="leiden")
        sc.tl.ingest(adataDS_2DS3, adataCON_IMR, obs="leiden")
        sc.tl.ingest(adataDS_DSP, adataCON_IMR, obs="leiden")

        adata = adataCON_IMR.concatenate(adataCON_H9)
        adata = adata.concatenate(adataCON_DS2U)
        adata = adata.concatenate(adataDS_2DS3)
        adata = adata.concatenate(adataDS_DSP)

        sc.pl.umap(adata, color=['sample', "group"], show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="CON_DS2U", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="CON_H9", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="CON_IMR", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="DS_2DS3", show=showPlots)
        sc.pl.umap(adata, color='sample', na_color="white",groups="DS_DSP", show=showPlots)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

        for fig in range(figureNum+1,  plt.gcf().number+1):
            pdf.savefig(figure=fig, bbox_inches='tight')
####################################################################################################
# Clustering Plots
####################################################################################################
print(bcolors.FAIL + "Saving Clustering Plots..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Clustering Plots '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.pca(adata, show=showPlots)
    plt.suptitle("PC1 & PC2", y=1, fontsize=25)
    sc.pl.pca_loadings(adata, show=showPlots)
    plt.suptitle("PCA Loadings", y=1, fontsize=25)
    sc.pl.pca_variance_ratio(adata, show=showPlots)
    plt.suptitle("PCA Variance Ratio", y=1, fontsize=25)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("Leiden Clusters", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', show=showPlots)
    plt.suptitle("UMAP by Sample", y=1, fontsize=25)
    sc.pl.umap(adata, color='group', groups="DS", show=showPlots)
    plt.suptitle("UMAP of DS", y=1, fontsize=25)
    sc.pl.umap(adata, color='group', show=showPlots)
    plt.suptitle("UMAP by Group", y=1, fontsize=25)
    sc.pl.umap(adata, color='n_genes_by_counts', cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=showPlots)
    plt.suptitle("Nº of Genes by Count", y=1, fontsize=25)
    sc.pl.umap(adata, color='pct_counts_mt', cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=showPlots)
    plt.suptitle("% MT Count", y=1, fontsize=25)
    sc.pl.pca(adata, color='sample', components='1,2', show=showPlots)
    plt.suptitle("PCA Components by Sample", y=1, fontsize=25)
    sc.pl.pca(adata, color='group', components='1,2', show=showPlots)
    plt.suptitle("PCA Components by Group", y=1, fontsize=25)
    sc.pl.pca(adata, color='sample', components='3,4', show=showPlots)
    plt.suptitle("PCA Components by Sample", y=1, fontsize=25)
    sc.pl.pca(adata, color='group', components='3,4', show=showPlots)
    plt.suptitle("PCA Components by Group", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_DS2U", show=showPlots)
    plt.suptitle("UMAP for CON_DS2U", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_H9", show=showPlots)
    plt.suptitle("UMAP for CON_H9", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_IMR", show=showPlots)
    plt.suptitle("UMAP for CON_IMR", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_2DS3", show=showPlots)
    plt.suptitle("UMAP for DS_2DS3", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_DSP", show=showPlots)
    plt.suptitle("UMAP for DS_DSP", y=1, fontsize=25)
    sc.pl.umap(adata, color=geneNames, cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=showPlots)
    plt.suptitle("UMAP for Marker Genes", y=1, fontsize=25)
    sc.pl.umap(adata, color=["DLX2", "SOX2"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), palette="tab20", show=showPlots)
    plt.suptitle("UMAP for Marker Genes", y=1, fontsize=25)
    sc.pl.umap(adata, color=["MKI67", "EOMES", "DLX2", "GLI3", "NEUROD6", "AQP4", "MEF2C"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=showPlots) #Marker Gene Expression
    plt.suptitle("Marker Gene Expression", y=1, fontsize=25)
    sc.pl.umap(adata, color=["PGK1", "ALDOA", "ARCN1", "GORASP2"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=showPlots) #Marker Stress Genes
    plt.suptitle("Gene Stress Markers", y=1, fontsize=25)
    sc.pl.umap(adata, color=["APP", "DYRK1A"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), groups="CON", show=showPlots)
    plt.suptitle("Overexpressed Chr21 Genes", y=1, fontsize=25)
    sc.pl.dotplot(adata, geneNames, groupby='leiden', show=showPlots)
    plt.suptitle("Dot Plot for Marker Genes", y=1, fontsize=25)

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

    getColumns = classify.findCellTypeIndividualCellTypes(["VIM"], [0])
    dfObj = pd.DataFrame(columns = getColumns["Cell Type"].to_list())

    topGenesPerCluster = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    topScores = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'scores']}).head(50)

    print("Printing topGenesPerCluster")
    print(topGenesPerCluster)
    newClusterNames = []
    newGroupClusters1 = []
    newGroupClusters2 = []
    newClusterNamesScore = []
    newGroupClusters1Score = []
    newGroupClusters2Score = []
    for column in topGenesPerCluster:
        genesList = topScores[str(column)+"_n"].to_list()
        scoreList = topScores[str(column)+"_s"].to_list()

        sc.pl.umap(adata, color=genesList[0:8], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=showPlots)
        plt.suptitle("Highly Expressed For Cluster "+str(column), y=1, fontsize=25)

        genesResult = classify.findCellTypeIndividualCellTypes(genesList, scoreList)
        groupResult = classify.findCellTypesGroup(genesResult, '% Genes')
        
        genesResult = genesResult.sort_values(['% Genes'], ascending = [False])
        newClusterNames.append(genesResult['Cell Type'].iloc[0])
        print(genesResult)

        groupResult = groupResult.sort_values(['Score'], ascending = [False])
        newGroupClusters1.append(groupResult['Cell Type'].iloc[0])
        print(groupResult)
        
        groupResult = groupResult.sort_values(['Sum'], ascending = [False])
        newGroupClusters2.append(groupResult['Cell Type'].iloc[0])
        print(groupResult)

        genesResult = classify.findCellTypeIndividualCellTypes(genesList, scoreList)
        groupResult = classify.findCellTypesGroup(genesResult, 'Marker Score')
        
        genesResult = genesResult.sort_values(['Marker Score'], ascending = [False])
        newClusterNamesScore.append(genesResult['Cell Type'].iloc[0])
        print(genesResult)

        groupResult = groupResult.sort_values(['Score'], ascending = [False])
        newGroupClusters1Score.append(groupResult['Cell Type'].iloc[0])
        print(groupResult)
        
        groupResult = groupResult.sort_values(['Sum'], ascending = [False])
        newGroupClusters2Score.append(groupResult['Cell Type'].iloc[0])
        print(groupResult)
        
    
    sc.pl.dotplot(adata, geneNames, groupby='leiden', show=showPlots)
    plt.suptitle("Dot Plot for Marker Genes"+str(column), y=1, fontsize=25)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP for Leiden"+str(column), y=1, fontsize=25)

    newClusterNames = classify.renameList(newClusterNames)
    newGroupClusters1 = classify.renameList(newGroupClusters1)
    newGroupClusters2 = classify.renameList(newGroupClusters2)
    newClusterNamesScore = classify.renameList(newClusterNamesScore)
    newGroupClusters1Score = classify.renameList(newGroupClusters1Score)
    newGroupClusters2Score = classify.renameList(newGroupClusters2Score)
    
    adata.rename_categories('leiden', newClusterNames)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP by Total Genes Found", y=1, fontsize=25)

    adata.rename_categories('leiden', newGroupClusters1)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP by Group Score ", y=1, fontsize=25)

    adata.rename_categories('leiden', newGroupClusters2)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP by Group Sum", y=1, fontsize=25)

    adata.rename_categories('leiden', newClusterNamesScore)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP by Total Genes Found Marker Score", y=1, fontsize=20)

    adata.rename_categories('leiden', newGroupClusters1Score)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP by Group Score Marker Score TRUST", y=1, fontsize=20)

    adata.rename_categories('leiden', newGroupClusters2Score)
    sc.pl.umap(adata, color='leiden', show=showPlots, legend_loc='on data')
    plt.suptitle("UMAP by Group Sum Marker Score", y=1, fontsize=20)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')


THE_FIGURE = plt.figure(figsize=(10, 10), dpi=500)
ax = plt.subplot(1, 1, 1)
sns.heatmap(dfObj, annot=False)
THE_FIGURE.savefig('image.pdf', bbox_inches='tight', pad_inches=0.1)

adata.write(results_file)