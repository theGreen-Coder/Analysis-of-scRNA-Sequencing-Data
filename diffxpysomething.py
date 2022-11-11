import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import pandas as pd
import scipy.stats
import scanpy as sc
from scipy.sparse import csr_matrix
import diffexpr as de


adata = sc.read("./dataSaveOriginal/rawDataset.h5ad") # ./dataSaveOriginal/rawDataset5000.h5ad
realAdata = sc.read("./dataSaveOriginal/rawDataset5000.h5ad") # ./dataSaveOriginal/rawDataset5000.h5ad
adata.var_names_make_unique()

adata.X = csr_matrix(adata.X)

adata.X = adata.X.todense()

adata.obs
adata.obs = adata.obs.rename(columns = {'group':'condition'})
adata.obs

realAdata.X = csr_matrix(realAdata.X)

realAdata.X = realAdata.X.todense()

realAdata.obs
realAdata.obs = realAdata.obs.rename(columns = {'group':'condition'})
realAdata.obs

data = anndata.AnnData(
  X=adata.X,
  var=adata.var,
  obs=adata.obs
)

realData = anndata.AnnData(
  X=realAdata.X,
  var=realAdata.var,
  obs=realAdata.obs
)

sc.pp.highly_variable_genes(data, n_top_genes=255, subset = True,
                           flavor = "seurat_v3", batch_key="sample")

test = de.test.wald(
    data=data,
    formula_loc="~ 1 + condition",
    factor_loc_totest="condition",
)

test.qval[:10]

sumar = test.summary().iloc[:10,:]

sumar.sort_values(['pval'], ascending = [True])
print(sumar)