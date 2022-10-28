import pandas as pd
import scanpy as sc
import numpy as np

####################################################################################################
# SetUp AnnData from Data 
####################################################################################################

def readData():
    adata = sc.read_10x_mtx(
        'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        cache=True)                              # write a cache file for faster subsequent reading
    
    return adata