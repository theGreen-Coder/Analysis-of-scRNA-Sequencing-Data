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

# Arguments Settings
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
print(adata)

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

print("Hello There!")
print(adata)

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

print(adata.obs)
print(adata.var)

####################################################################################################
# QC Filtering 
####################################################################################################
print(bcolors.FAIL + "Applying Quality Control Filters..." + bcolors.ENDC)
print("Filtering Damaged Cells")
# sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
adata = adata[adata.obs.n_genes_by_counts < upper_lim]
adata = adata[adata.obs.pct_counts_mt < 10, :]
# adata = adata[adata.obs.n_genes_by_counts > 1000, :]
####################################################################################################

####################################################################################################
# Clustering 
####################################################################################################
print(bcolors.FAIL + "Performing Clustering..." + bcolors.ENDC)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
sc.tl.umap(adata)
sc.external.pp.bbknn(adata, batch_key='sample')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.7)
sc.pl.umap(adata, color="leiden", legend_loc='on data', show=True)

####################################################################################################
# Finding Marker Genes
####################################################################################################
print(bcolors.FAIL + "Finding Marker Genes..." + bcolors.ENDC)

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
topGenesPerCluster = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
topScores = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'scores']}).head(50)
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


newClusterNames = classify.renameList(newClusterNames)
newGroupClusters1 = classify.renameList(newGroupClusters1)
newGroupClusters2 = classify.renameList(newGroupClusters2)
newClusterNamesScore = classify.renameList(newClusterNamesScore)
newGroupClusters1Score = classify.renameList(newGroupClusters1Score)
newGroupClusters2Score = classify.renameList(newGroupClusters2Score)
adata.rename_categories('leiden', newGroupClusters1Score)
sc.pl.umap(adata, color="leiden", legend_loc='on data', show=True)

print(bcolors.FAIL + "Final Restuls..." + bcolors.ENDC)
outputAdata = sc.read("./output/rawDatasetSamples.h5ad")
print("After Reading!")
print(outputAdata)
outputAdata.obs["leiden"] = adata.obs["leiden"]
print("Adata")
print(adata)
print("OutputData")
print(outputAdata)
outputAdata.write("./output/savedDataClustersFinal.h5ad")
print("Saved File to Results_File!")