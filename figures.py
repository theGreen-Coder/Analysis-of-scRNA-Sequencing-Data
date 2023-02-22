from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import modules.classifyClusters.classifyClusters as classify
import os
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


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
adata = sc.read("./output/savedDataClustersFinal.h5ad")
adata.var_names_make_unique()
showPlots = False

####################################################################################################
# QC Calculations 
####################################################################################################
print(bcolors.FAIL + "Performing Quality Control..." + bcolors.ENDC)
adata.var['mt'] = adata.var_names.str.startswith('MT-')

riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 10, :]
adata = adata[adata.obs.n_genes_by_counts > 1000, :]

####################################################################################################
# Clustering 
####################################################################################################
print(bcolors.FAIL + "Performing Clustering..." + bcolors.ENDC)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.var.at["OLIG2", "highly_variable"] = True

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
sc.tl.umap(adata)

adataDS =  adata[adata.obs["group"] == "DS"]
adataCON =  adata[adata.obs["group"] == "CON"]

sc.external.pp.bbknn(adata, batch_key='sample')
sc.tl.umap(adata)

####################################################################################################
# Plot Umap
####################################################################################################
print(bcolors.FAIL + "Plotting UMAP..." + bcolors.ENDC)

with PdfPages('./Figures/Figure 1.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.umap(adata, color="leiden", show=showPlots)
    sc.pl.umap(adata, color='group', groups="CON", show=showPlots)
    plt.suptitle("UMAP of DS", y=1, fontsize=25)
    sc.pl.umap(adata, color='group', groups="DS", show=showPlots)
    plt.suptitle("UMAP of DS", y=1, fontsize=25)
    sc.pl.umap(adata, color='group', show=showPlots)
    plt.suptitle("UMAP by Group", y=1, fontsize=25)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

####################################################################################################
# Frequency Analysis 
####################################################################################################
print(bcolors.FAIL + "Frequency Analysis..." + bcolors.ENDC)
with PdfPages('./Figures/Frequency Analysis.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    ####################################################################################################
    # Grouped Frequency Analysis 
    ####################################################################################################
    leidenNames = list(adata.obs["leiden"])
    leidenNewNames = []

    for name in leidenNames:
        if "_" in name:
            newName = name.split("_")[0]
            leidenNewNames.append(newName)
        else:
            leidenNewNames.append(name)
    adata.obs["leiden"] = leidenNewNames

    ####################################################################################################
    # Grouped Frequency Analysis 
    ####################################################################################################
    num_tot_cells = adata.obs.groupby(['group']).count()
    num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.nFeature_RNA))

    cell_type_counts = adata.obs.groupby(['group', 'sample','leiden']).count()
    cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
    cell_type_counts = cell_type_counts[cell_type_counts.columns[0:6]]

    cell_type_counts['total_cells'] = cell_type_counts.group.map(num_tot_cells).astype(int)

    cell_type_counts['frequency'] = cell_type_counts.nFeature_RNA / cell_type_counts.total_cells

    cell_type_counts = cell_type_counts.sort_values(by=['leiden'])
    print(cell_type_counts)
    
    plt.figure(figsize = (10,4))

    print(cell_type_counts)
    ax = sns.boxplot(data = cell_type_counts, x = 'leiden', y = 'frequency', hue = 'group')

    plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')


    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')