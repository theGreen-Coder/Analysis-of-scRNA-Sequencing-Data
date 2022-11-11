import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import pandas as pd
import scipy.stats
import scanpy as sc
from scipy.sparse import csr_matrix
import diffxpy.api as de

def main():
    print("Reading!")
    adata = sc.read("./dataSaveOriginal/rawDataset5000.h5ad") # ./dataSaveOriginal/rawDataset5000.h5ad
    adata.var_names_make_unique()

    adata.X = csr_matrix(adata.X)
    print("To Dense!")
    adata.X = adata.X.todense()

    print(adata.obs)
    adata.obs = adata.obs.rename(columns = {'group':'condition'})
    print(adata.obs)

    print("Data!")
    data = anndata.AnnData(
    X=adata.X,
    var=adata.var,
    obs=adata.obs
    )

    print("Highly Variable Genes!")
    sc.pp.highly_variable_genes(data, n_top_genes=1000, subset = True, flavor = "seurat_v3", batch_key="sample")


    print("Wald Test!")
    test = de.test.wald(
        data=data,
        formula_loc="~ 1 + condition",
        factor_loc_totest="condition",
    )
    sumar = test.summary()
    print(test.qval[:10])
    print(sumar)
    sumar.to_csv("./diffxpyLogs.csv", sep='\t')


if __name__ == '__main__':
    main()