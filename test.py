import pandas as pd
import scanpy as sc

print("Used Adata for DE:")
adataDE = sc.read_h5ad('./output/savedDataClustersFinal.h5ad')
sc.pp.normalize_total(adataDE, target_sum=10000)
sc.pp.log1p(adataDE)
sc.pp.highly_variable_genes(adataDE, min_mean=0.0125, max_mean=3, min_disp=0.5)

print(adataDE.var)

adataDE.var.at["OR4F5", "highly_variable"] = True

print(adataDE.var)

