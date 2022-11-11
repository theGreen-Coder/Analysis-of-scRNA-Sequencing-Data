from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import tests.loadScanpy as loadScanpy
import classifyClusters as classify
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

inputData = "./output/savedData.h5ad" # ./dataSaveOriginal/rawDataset5000.h5ad
results_file = './output/savedData.h5ad'
outputDirectory = "./outputPDFs/"

# Arguments Settings

sysArgs = sys.argv

showPlots = False
if "-showPlots" in sysArgs:
    showPlots = True
else:
    showPlots = False

arg = ""
for word in sysArgs:
    if word != "main.py" and word != "-showPlots":
        arg = arg + word + " "
        # if "-" not in word:
        #     OUTPUT_DIR = outputDirectory+word+"/"
        #     CHECK_FOLDER = os.path.isdir(OUTPUT_DIR)
        #     if not CHECK_FOLDER:
        #         os.makedirs(OUTPUT_DIR)
        #     outputDirectory = OUTPUT_DIR
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
if "-loadGenomics" in sysArgs:
    adata = loadScanpy.readData()
    print(adata)
else:
    adata = sc.read(inputData)
adata.var_names_make_unique()

geneMarkers = pd.read_csv('./dataSaveOriginal/cellTypeMarkers.csv')
geneNames = geneMarkers["gene"].tolist()

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP

####################################################################################################
# Trajectory Analysis 
####################################################################################################
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden', legend_loc='on data')

sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, threshold=0.03, show=True)

sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color="leiden", legend_loc='on data')



