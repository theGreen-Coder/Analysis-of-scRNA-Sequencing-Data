import scanpy as sc
import numpy as np
adata = sc.read("./dataSaveOriginal/rawDataset.h5ad")
adata.var_names_make_unique()
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

import loompy as lp
row_attrs = {
"Gene": np.array(adata.var_names),
}
col_attrs = {
"CellID": np.array(adata.obs_names),
"nGene": np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(), "nUMI": np.array(np.sum(adata.X.transpose(),axis=0)).flatten(),
}

lp.create("rawDataset_filtered.loom",adata.X.transpose(),row_attrs, col_attrs)
