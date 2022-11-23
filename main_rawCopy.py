from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import tests.loadScanpy as loadScanpy
import modules.classifyClusters.classifyClusters as classify
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
adata = sc.read("./dataSaveOriginal/rawDataset.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')
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
adata = adata.concatenate(adataDS_2DS3)
adata = adata.concatenate(adataDS_DSP)

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP

####################################################################################################
# QC Calculations 
####################################################################################################
print(bcolors.FAIL + "Performing Quality Control..." + bcolors.ENDC)
# Mark Mitochondrial Genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Mark Ribosomal Genes
riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)

# Calculate QC Metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
adata = adata[adata.obs.n_genes_by_counts < upper_lim]
adata = adata[adata.obs.pct_counts_mt < 10, :]

####################################################################################################
# Clustering 
####################################################################################################
print(bcolors.FAIL + "Performing Clustering..." + bcolors.ENDC)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.25)
sc.external.pp.bbknn(adata, batch_key='sample')
sc.tl.umap(adata)

####################################################################################################
# Clustering Plots
####################################################################################################
print(bcolors.FAIL + "Saving Clustering Plots..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Clustering Plots Main RawCopy.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.pca(adata, show=False)
    plt.suptitle("PC1 & PC2", y=1, fontsize=25)
    sc.pl.pca_loadings(adata, show=False)
    plt.suptitle("PCA Loadings", y=1, fontsize=25)
    sc.pl.pca_variance_ratio(adata, show=False)
    plt.suptitle("PCA Variance Ratio", y=1, fontsize=25)
    sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
    plt.suptitle("Leiden Clusters", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', show=False)
    plt.suptitle("UMAP by Sample", y=1, fontsize=25)
    sc.pl.umap(adata, color='group', groups="DS", show=False)
    plt.suptitle("UMAP of DS", y=1, fontsize=25)
    sc.pl.umap(adata, color='group', show=False)
    plt.suptitle("UMAP by Group", y=1, fontsize=25)
    sc.pl.umap(adata, color='n_genes_by_counts', cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
    plt.suptitle("Nº of Genes by Count", y=1, fontsize=25)
    sc.pl.umap(adata, color='pct_counts_mt', cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
    plt.suptitle("% MT Count", y=1, fontsize=25)
    sc.pl.pca(adata, color='sample', components='1,2', show=False)
    plt.suptitle("PCA Components by Sample", y=1, fontsize=25)
    sc.pl.pca(adata, color='group', components='1,2', show=False)
    plt.suptitle("PCA Components by Group", y=1, fontsize=25)
    sc.pl.pca(adata, color='sample', components='3,4', show=False)
    plt.suptitle("PCA Components by Sample", y=1, fontsize=25)
    sc.pl.pca(adata, color='group', components='3,4', show=False)
    plt.suptitle("PCA Components by Group", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_DS2U", show=False)
    plt.suptitle("UMAP for CON_DS2U", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_H9", show=False)
    plt.suptitle("UMAP for CON_H9", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="CON_IMR", show=False)
    plt.suptitle("UMAP for CON_IMR", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_2DS3", show=False)
    plt.suptitle("UMAP for DS_2DS3", y=1, fontsize=25)
    sc.pl.umap(adata, color='sample', na_color="white",groups="DS_DSP", show=False)
    plt.suptitle("UMAP for DS_DSP", y=1, fontsize=25)
    sc.pl.umap(adata, color=geneNames, cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
    plt.suptitle("UMAP for Marker Genes", y=1, fontsize=25)
    sc.pl.umap(adata, color=["DLX2", "SOX2"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), palette="tab20", show=False)
    plt.suptitle("UMAP for Marker Genes", y=1, fontsize=25)
    sc.pl.umap(adata, color=["MKI67", "EOMES", "DLX2", "GLI3", "NEUROD6", "AQP4", "MEF2C"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False) #Marker Gene Expression
    plt.suptitle("Marker Gene Expression", y=1, fontsize=25)
    sc.pl.umap(adata, color=["PGK1", "ALDOA", "ARCN1", "GORASP2"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False) #Marker Stress Genes
    plt.suptitle("Gene Stress Markers", y=1, fontsize=25)
    sc.pl.umap(adata, color=["APP", "DYRK1A"], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), groups="CON", show=False)
    plt.suptitle("Overexpressed Chr21 Genes", y=1, fontsize=25)
    sc.pl.dotplot(adata, geneNames, groupby='leiden', show=False)
    plt.suptitle("Dot Plot for Marker Genes", y=1, fontsize=25)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP
# DS CON

####################################################################################################
# Finding Marker Genes
####################################################################################################
print(bcolors.FAIL + "Finding Marker Genes..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Marker Genes Raw.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=50, sharey=False, show=False)

    # sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
    # sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
    # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)

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
    counter = 0
    for column in topGenesPerCluster:
        genesList = topScores[str(column)+"_n"].to_list()
        scoreList = topScores[str(column)+"_s"].to_list()

        sc.pl.umap(adata, color=genesList[0:8], cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
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

        groupResult = groupResult.sort_values(['Cell Type'], ascending = [False])
        if(counter == 0):
            dfObj = pd.DataFrame(columns = groupResult["Cell Type"].to_list())
            dfObj.loc[counter] = groupResult["Score"].to_list()
        else:
            dfObj.loc[counter] = groupResult["Score"].to_list()

        counter += 1
        
        genesResult = genesResult.sort_values(['Marker Score'], ascending = [False])
        newClusterNamesScore.append(genesResult['Cell Type'].iloc[0])
        print(genesResult)

        groupResult = groupResult.sort_values(['Score'], ascending = [False])
        newGroupClusters1Score.append(groupResult['Cell Type'].iloc[0])
        print(groupResult)
        
        groupResult = groupResult.sort_values(['Sum'], ascending = [False])
        newGroupClusters2Score.append(groupResult['Cell Type'].iloc[0])
        print(groupResult)
    
    sc.pl.dotplot(adata, geneNames, groupby='leiden', show=False)
    plt.suptitle("Dot Plot for Marker Genes"+str(column), y=1, fontsize=25)
    sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
    plt.suptitle("UMAP for Leiden"+str(column), y=1, fontsize=25)

    newClusterNames = classify.renameList(newClusterNames)
    newGroupClusters1 = classify.renameList(newGroupClusters1)
    newGroupClusters2 = classify.renameList(newGroupClusters2)
    newClusterNamesScore = classify.renameList(newClusterNamesScore)
    newGroupClusters1Score = classify.renameList(newGroupClusters1Score)
    newGroupClusters2Score = classify.renameList(newGroupClusters2Score)
    if "-minimal" not in sysArgs:
        adata.rename_categories('leiden', newClusterNames)
        sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
        plt.suptitle("UMAP by Total Genes Found", y=1, fontsize=25)

        adata.rename_categories('leiden', newGroupClusters1)
        sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
        plt.suptitle("UMAP by Group Score ", y=1, fontsize=25)

        adata.rename_categories('leiden', newGroupClusters2)
        sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
        plt.suptitle("UMAP by Group Sum", y=1, fontsize=25)

        adata.rename_categories('leiden', newClusterNamesScore)
        sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
        plt.suptitle("UMAP by Total Genes Found Marker Score", y=1, fontsize=20)

        adata.rename_categories('leiden', newGroupClusters1Score)
        sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
        plt.suptitle("UMAP by Group Score Marker Score TRUST", y=1, fontsize=20)

        adata.rename_categories('leiden', newGroupClusters2Score)
        sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
        plt.suptitle("UMAP by Group Sum Marker Score", y=1, fontsize=20)

    adata.rename_categories('leiden', newGroupClusters1Score)
    plt.figure(figsize = (10,4))
    ax = plt.subplot(1, 1, 1)
    print(dfObj)
    sns.heatmap(dfObj, annot=True)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

print("Saved File to Results_File!")
adata.write(results_file)

####################################################################################################
# Differential Expression Using Scanpy
####################################################################################################
print(bcolors.FAIL + "Performing Differential Expression..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Differential Expression Raw.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups_heatmap(adata, groupby="leiden", n_genes=10, show=False)
    sc.pl.rank_genes_groups_heatmap(adata, groupby="group", n_genes=10, show=False)
    sc.tl.rank_genes_groups(adata, 'group', method='t-test')
    sc.pl.rank_genes_groups_heatmap(adata, groupby="group", n_genes=10, show=False)
    plt.suptitle("UMAP by Group Score Marker Score TRUST", y=1, fontsize=20)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

####################################################################################################
# Frequency Analysis 
####################################################################################################
print(bcolors.FAIL + "Frequency Analysis..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Frequency Analysis Raw.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    num_tot_cells = adata.obs.groupby(['group']).count()
    num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.nFeature_RNA))

    cell_type_counts = adata.obs.groupby(['group', 'sample','leiden']).count()
    cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
    cell_type_counts = cell_type_counts[cell_type_counts.columns[0:6]]

    cell_type_counts['total_cells'] = cell_type_counts.group.map(num_tot_cells).astype(int)

    cell_type_counts['frequency'] = cell_type_counts.nFeature_RNA / cell_type_counts.total_cells

    cell_type_counts = cell_type_counts.sort_values(by=['leiden'])
    
    plt.figure(figsize = (10,4))

    ax = sns.boxplot(data = cell_type_counts, x = 'leiden', y = 'nCount_RNA', hue = 'group')

    plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')


    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')