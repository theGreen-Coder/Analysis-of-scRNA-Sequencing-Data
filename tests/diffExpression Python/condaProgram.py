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
from matplotlib.backends.backend_pdf import PdfPages

def main():
    print("Reading!")
    adata = sc.read("./dataSaveOriginal/rawDataset5000.h5ad") # ./dataSaveOriginal/rawDataset5000.h5ad
    adata.var_names_make_unique()

    genesHSA21 = pd.read_csv("./dataSaveOriginal/HSA21genes.csv")
    genesHSA21List = []

    for gene in genesHSA21["hgnc_symbol"].tolist():
        if type(gene) == str:
            genesHSA21List.append(gene)

    adata.X = adata.raw.X

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
    sc.pp.highly_variable_genes(data, n_top_genes=10000, subset = True, flavor = "seurat_v3", batch_key="sample")


    print("Wald Test!")
    test = de.test.wald(
        data=data,
        formula_loc="~ 1 + condition",
        factor_loc_totest="condition",
    )
    sumar = test.summary().sort_values('log2fc', ascending = False).reset_index(drop = True)
    sumar['log2fc'] = sumar['log2fc']*-1
    sumar = sumar.sort_values('log2fc', ascending = False).reset_index(drop = True)
    sumar = sumar[(sumar.qval < 0.05) & (abs(sumar.log2fc) > .5)]
    sumar = sumar[sumar['mean'] > 0.15]

    sumar.to_csv("./diffxpyLogs.csv")

    genes_to_show = sumar[-25:].gene.tolist() + sumar[:25].gene.tolist() #top 25 and bottom 25 from sorted df

    print(len(genes_to_show))
    

    with PdfPages('./outputPDFs/Frequency Analysis diffxpy.pdf') as pdf:
        dummyPlot = plt.plot([1, 2])
        figureNum = plt.gcf().number

        sc.pl.heatmap(data, genes_to_show, groupby='condition', vmax=5, swap_axes=True, show=False)
        sc.pl.heatmap(data, genes_to_show, groupby='sample', vmax=5, swap_axes=True, show=False)
        sc.pl.heatmap(adata, genesHSA21List, groupby="condition", vmax=5, show=False)
        sc.pl.heatmap(adata, genesHSA21List, groupby="sample", vmax=5, show=False)

        for fig in range(figureNum+1,  plt.gcf().number+1):
            pdf.savefig(figure=fig, bbox_inches='tight')

    


if __name__ == '__main__':
    main()