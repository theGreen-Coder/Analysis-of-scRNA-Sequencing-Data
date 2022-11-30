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

inputData = "./output/savedData.h5ad" # ./dataSaveOriginal/rawDataset5000.h5ad
results_file = './output/savedData.h5ad'
outputDirectory = "./outputPDFs/"

def getCountsSample():
    adata = sc.read("./output/savedDataClustersFinal.h5ad")
    adata.var_names_make_unique()
    adata.X = adata.X.astype('float64')
    print("Printing Raw Data")
    print(adata.X)

    buckets = [0] * 36601

    rows = adata.var["features"].values

    # CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP
    listLeiden = adata.obs["leiden"].unique().tolist()

    df = pd.DataFrame(rows, columns=['id'])

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

    for item in listLeiden:
        print("1")
        df[item+".CON_DS2U"] = getColSample("CON_DS2U", item)
        print("2")
        df[item+".CON_H9"] = getColSample("CON_H9", item)
        print("3")
        df[item+".CON_IMR"] = getColSample("CON_IMR", item)
        print("4")
        df[item+".DS_2DS3"] = getColSample("DS_2DS3", item)
        print("5")
        df[item+".DS_DSP"] = getColSample("DS_DSP", item)


    df = df.set_index('id')
    print(df)

    df.to_csv('./output/rawCountsDESeqPerClusterNew.csv')

def getCountsCONvsDS():
    filePath = './outputPDFs/DiffExpOut/normTransformed.csv'
    if os.path.exists(filePath):
        data = pd.read_csv(filePath)
        data = data.set_index(["Unnamed: 0"])

        print(data)

        controlData = pd.DataFrame(index=data.index)
        downData = pd.DataFrame(index=data.index)

        for column in data:
            if("CON" in column):
                controlData[column] = data[column]
            elif("DS" in column):
                downData[column] = data[column]

        print(controlData)
        print(downData)

        controlAvgData = pd.DataFrame(index=data.index)
        downAvgData = pd.DataFrame(index=data.index)

        count = 0
        for column in controlData:
            count += 1
            if count == 1:
                sampleDF = pd.DataFrame([], index=data.index)
                sampleDF[column] = controlData[column]
            elif count == 2:
                sampleDF[column] = controlData[column]
            elif count == 3:
                sampleDF[column] = controlData[column]
                sampleDF['avg'] = sampleDF.mean(axis=1)
                colName = column.split(".")[0]
                controlAvgData[colName] = sampleDF['avg']

                sampleDF = sampleDF.iloc[0:0]
                count=0
            
        count = 0
        sampleDF = sampleDF.iloc[0:0]
        for column in downData:
            count += 1
            if count == 1:
                sampleDF = pd.DataFrame([], index=data.index)
                sampleDF[column] = downData[column]
            elif count == 2:
                sampleDF[column] = downData[column]
                sampleDF['avg'] = sampleDF.mean(axis=1)
                colName = column.split(".")[0]
                downAvgData[colName] = sampleDF['avg']

                sampleDF = sampleDF.iloc[0:0]
                count=0

        print(controlAvgData) 
        print(downAvgData)

        controlAvgData.to_csv("./output/NormAvgCountsCONTROL.csv")
        downAvgData.to_csv("./output/NormAvgCountsDOWN.csv")

        diffAvgData = controlAvgData

        for rowIndex, row in controlAvgData.iterrows(): #iterate over rows
            for columnIndex, value in row.items():
                valueCON = controlAvgData.loc[rowIndex, columnIndex]
                valueDS = downAvgData.loc[rowIndex, columnIndex]
                # print(valueCON)
                # print(valueDS)
                if valueCON > valueDS:
                    diffAvgData.loc[rowIndex, columnIndex] = valueCON - valueDS
                elif valueCON < valueDS:
                    diffAvgData.loc[rowIndex, columnIndex] = valueCON - valueDS
                else:
                    print("WTF just happened?!")

        print(diffAvgData)
        diffAvgData.to_csv("./output/NormDiffCountsCONvsDS.csv")

getCountsSample()
getCountsCONvsDS()