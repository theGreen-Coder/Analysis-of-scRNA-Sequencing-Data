from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import tests.loadScanpy as loadScanpy
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

inputData = "./output/savedData.h5ad" # ./dataSaveOriginal/rawDataset5000.h5ad
results_file = './output/savedData.h5ad'
outputDirectory = "./outputPDFs/"

adata = sc.read("./dataSaveOriginal/rawDatasetClusters.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')
print("Printing Raw Data")
print(adata.X)

buckets = [0] * 36601

def getColSample(sample, leiden):
    if leiden == "":
        subset =  adata[adata.obs["sample"] == sample]
        returnArray = [sum(x) for x in zip(*subset.X)]
        return returnArray[0].toarray()[0] #returnArray[0].toarray()[0]
    else:
        subset =  adata[adata.obs["sample"] == sample]
        print(subset)
        subset = subset[subset.obs["leiden"] == leiden]
        print(subset)
        returnArray = [sum(x) for x in zip(*subset.X)]
        print(returnArray)
        if len(returnArray) == 0:
            print(buckets)
            return buckets
        print(len(returnArray[0].toarray()[0]))
        return returnArray[0].toarray()[0] #

rows = adata.var["features"].values

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP
listLeiden = adata.obs["leiden"].unique().tolist()

df = pd.DataFrame(rows, columns=['id'])

for item in listLeiden:
    print("1")
    df["CON_DS2U"+"."+item] = getColSample("CON_DS2U", item)
    print("2")
    df["CON_H9"+"."+item] = getColSample("CON_H9", item)
    print("3")
    df["CON_IMR"+"."+item] = getColSample("CON_IMR", item)
    print("4")
    df["DS_2DS3"+"."+item] = getColSample("DS_2DS3", item)
    print("5")
    df["DS_DSP"+"."+item] = getColSample("DS_DSP", item)


df = df.set_index('id')
print(df)

df.to_csv('./output/rawCountsDESeqPerCluster.csv')