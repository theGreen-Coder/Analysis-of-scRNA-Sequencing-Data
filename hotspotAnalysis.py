import hotspot
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

inputDataAdata = "./output/savedDataClustersFinal.h5ad"
inputDataAdata = "./dataSaveOriginal/rawDataset1000.h5ad"

adata = sc.read(inputDataAdata) 
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Mark Ribosomal Genes
riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)

# Calculate QC Metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_genes(adata, min_counts=3)

sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.3)

adataHotspot = sc.read(inputDataAdata)
adataHotspot.var_names_make_unique()
adataHotspot.X = adataHotspot.X.astype('float64')
adataHotspot.obsm["X_pca"] = adata.obsm["X_pca"]
print(adataHotspot)
sc.pp.filter_genes(adataHotspot, min_counts=50)
print(adataHotspot)

hs = hotspot.Hotspot(
    adataHotspot,
    model='danb',
    latent_obsm_key="X_pca",
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

hs_results = hs.compute_autocorrelations()

print(hs.results)

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index # Select genes
hs.compute_local_correlations(hs_genes) # jobs

print("Outputing HS_GENE!")
print(hs_genes)

modules = hs.create_modules(
    min_gene_threshold=30, core_only=True, fdr_threshold=0.05
)

hs.plot_local_correlations()

print(hs.plot_local_correlations())