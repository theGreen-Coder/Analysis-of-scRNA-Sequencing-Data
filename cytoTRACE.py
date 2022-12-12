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
results_file_group = './output/savedDataGrouped.h5ad'
outputDirectory = "./outputPDFs/"

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'    

####################################################################################################
# Import Data 
####################################################################################################
print(bcolors.FAIL + "Importing Data..." + bcolors.ENDC)
adata = sc.read("./output/savedDataClustersFinal.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

print(adata)

# Generate Expression Matrix
expressionMatrix = pd.DataFrame.sparse.from_spmatrix(adata.X)

cellsIDs = list(adata.obs.index)

geneIDs = list(adata.var.index)

expressionMatrix = expressionMatrix.T

expressionMatrix.columns = cellsIDs
expressionMatrix.index = geneIDs

expressionMatrix.to_csv("./output/cytoTRACE/expressionMatrix.csv")

# Phenotype and Batch
phenotypeDataframe = pd.DataFrame(adata.obs["leiden"])
sampleDataframe = pd.DataFrame(adata.obs["sample"])
phenotypeDataframe.to_csv("./output/cytoTRACE/phenotypeDataframe.csv")
sampleDataframe.to_csv("./output/cytoTRACE/sampleDataframe.csv")