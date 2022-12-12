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

####################################################################################################
# Global Settings 
####################################################################################################
np.set_printoptions(threshold=sys.maxsize)
pd.options.display.max_columns = sys.maxsize
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

inputData = "./dataSaveOriginal/rawDataset.h5ad" # ./dataSaveOriginal/rawDataset5000.h5ad

adata = sc.read("./output/savedDataClustersFinal.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

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

print(adata)

geneMarkers = pd.read_csv('./dataSaveOriginal/cellTypeMarkers.csv')
geneNames = geneMarkers["gene"].tolist()

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP

####################################################################################################
# QC Calculations 
####################################################################################################
# Mark Mitochondrial Genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Mark Ribosomal Genes
riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)

# Calculate QC Metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)

with PdfPages('./outputPDFs/Doublets Plots savedData.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    print(adata)
    sc.external.pp.scrublet(adata)
    sc.external.pl.scrublet_score_distribution(adata, show=False)

    print(adata.obs)
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
    sc.external.pp.bbknn(adata, batch_key='sample')
    sc.tl.umap(adata)

    # Pre-Filtering
    sc.pl.umap(adata, color='doublet_score', show=False)
    adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(str).astype('category')
    sc.pl.umap(adata, color='predicted_doublet', show=False)

    # Filtering
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    adata = adata[adata.obs.n_genes_by_counts > 1000, :]

    # Post-filtering
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='doublet_score', show=False)
    sc.pl.umap(adata, color='predicted_doublet', show=False)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

print(adata)