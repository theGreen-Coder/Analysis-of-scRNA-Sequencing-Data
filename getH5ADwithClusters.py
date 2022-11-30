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
results_file = './output/savedData.h5ad'
outputDirectory = "./outputPDFs/"

# Arguments Settings

sysArgs = sys.argv

class bcolors:
    FAIL = '\033[91m'
    ENDC = '\033[0m'  

####################################################################################################
# Import Data 
####################################################################################################
print(bcolors.FAIL + "Importing Data..." + bcolors.ENDC)
adata = sc.read("./dataSaveOriginal/rawDataset.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

print(adata)

# Remove Outlier Samples from the Data
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

####################################################################################################
# QC Calculations 
####################################################################################################
print(bcolors.FAIL + "Performing Quality Control..." + bcolors.ENDC)
# Mark Mitochondrial Genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Mark Ribosomal Genes
riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)

# Calculate QC Metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)

####################################################################################################
# QC Filtering 
####################################################################################################
print(bcolors.FAIL + "Applying Quality Control Filters..." + bcolors.ENDC)
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.pct_counts_mt < 10, :]
adata = adata[adata.obs.n_genes_by_counts > 1000, :]

####################################################################################################
# Output Raw File with Cluster Column "Leiden" 
####################################################################################################
print(bcolors.FAIL + "Output Raw File with Cluster Column Leiden..." + bcolors.ENDC)
savedData = sc.read("./output/savedData.h5ad")
print("SAVED DATA CELLS:")
print(savedData)
print("OUTPUT DATA CELLS:")
print(adata)

adata.obs["leiden"] = savedData.obs["leiden"]
print("OUTPUT DATA CELLS after ADDING COLUMN:")
print(adata)
adata.write("./output/savedDataClustersFinal.h5ad")
print("Saved File to Results_File!")

####################################################################################################
# Output Raw File with Cluster Column "Leiden" Grouped
####################################################################################################
print(bcolors.FAIL + "Output Raw File with Cluster Column Leiden..." + bcolors.ENDC)
savedDataGrouped = sc.read("./output/savedDataGrouped.h5ad")
print("SAVED DATA CELLS:")
print(savedDataGrouped)
print("OUTPUT DATA CELLS:")
print(adata)

adata.obs["leiden"] = savedDataGrouped.obs["leiden"]
print("OUTPUT DATA CELLS after ADDING COLUMN:")
print(adata)
adata.write("./output/savedDataGroupedClustersFinal.h5ad")
print("Saved File to Results_File!")