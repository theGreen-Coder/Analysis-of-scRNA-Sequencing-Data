import pandas as pd
import scanpy as sc

adata = sc.read("../output/savedDataClustersFinal.h5ad")
print(adata.var)

geneNames = adata.var["features"].to_list()

output = ""
output = '\n'.join(geneNames)

with open('allGenesNames.txt', 'w') as f:
    f.write(output)