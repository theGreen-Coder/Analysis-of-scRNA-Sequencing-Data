import scanpy as sc
import numpy as np
adata = sc.read("./rawDataset.h5ad")
adata.var_names_make_unique()

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

# compute the number of genes per cell (computes ‘n_genes' column) 
sc.pp.filter_cells(adata, min_genes=0)
# mito and genes/counts cuts
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes 
adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
# add the total counts per cell as observations-annotation to adata 
adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))

sc.pp.filter_cells(adata, min_genes=200) 
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['n_genes'] < 4000, :] 
adata = adata[adata.obs['percent_mito'] < 0.15, :]

print(adata)

import pandas as pd
significantGenes = list(pd.read_csv("./significantGenes.csv")["genes"])
significantGenes
varGenes = list(adata.var.index.values)
significantGeneList = []
for gene in varGenes:
    if gene in significantGenes:
        significantGeneList.append(True)
    else:
        significantGeneList.append(False)

adata.var["significantGenes"] = significantGeneList
adata = adata[:, adata.var.significantGenes]

print(adata)

import loompy as lp
row_attrs = {
"Gene": np.array(adata.var_names),
}
col_attrs = {
"CellID": np.array(adata.obs_names),
"nGene": np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(), "nUMI": np.array(np.sum(adata.X.transpose(),axis=0)).flatten(),
}

lp.create("rawDataset_filtered_significant.loom",adata.X.transpose(),row_attrs, col_attrs)