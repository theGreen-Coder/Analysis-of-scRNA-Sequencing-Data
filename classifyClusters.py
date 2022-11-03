import pandas as pd
import numpy as np
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
pd.options.display.max_columns = sys.maxsize

# Genes For Cluster 9
genes = ["CENPF", "NUSAP1", "HMGB2", "PTTG1", "TOP2A", "TPX2", "VIM", "CCNB1", "DLGAP5", "SMC4", "H2AFX", "HMGN2", "CCNA2", "CKS2", "MKI67", "BIRC5", "ASPM", "CENPE", "PBK", "NUF2", "GTSE1", "PRC1", "SGO2", "CKAP2", "CCNB2"]

def renameList(listEx):
    renamedList = []
    for i in listEx:
        if i in renamedList:
            isInList = True
            count = 1
            newName = i
            while isInList:
                if newName in renamedList:
                    newName = i+"_"+str(count)
                    print(newName)
                    count +=1
                else:
                    renamedList.append(newName)
                    isInList = False
        else:
            renamedList.append(i)
    return renamedList

def findCellTypeIndividualCellTypes(listGenes):
    results = pd.DataFrame(data={'Cell Type': [], 'Scoring': [], 'Total Average Difference Genes': [], 'Row Number': [], 'Average Row Number': [], 'Score Number': [], 'Genes Found': []})
    row = 0
    for filename in os.listdir("./clusterGeneNames/files/"):
        df = pd.read_csv('./clusterGeneNames/files/'+filename, sep='\t', header=0)
        df = df.head(500)
        lenDF = len(df.index)
        i = 0
        j = 1
        avgDiff = 0
        totalRow = 0
        score = 0
        for gene in listGenes:
            if gene in df["id"].to_list():
                rowDataframe = df.loc[df['id'] == gene]
                avgDiff = float(rowDataframe["avg_diff|float"])
                totalRow += j
                i+=1
            j+=1
        if(i!=0.0): 
            averageRow = float(totalRow/i)
            score = float((i)/averageRow)
        else: 
            averageRow = 0
            score = 0
        results.loc[row] = [filename.replace('.tsv', ''), float((i/lenDF)*100), avgDiff, totalRow, averageRow, score, i]
        row +=1
    return results

def findCellTypesGroup(pdData):
    results = pd.DataFrame(data={'Cell Type': [], 'Sum': [], 'Num': [], 'Score': [],})
    groups = pd.read_csv('./clusterGeneNames/cellTypeToGroup.csv', header=0)
    for index, row in pdData.iterrows():
        groupID = groups.loc[groups["cellType"] == row["Cell Type"]]["group"].to_list()[0]
        if(groupID in results.values):
            indexValue = results.loc[results['Cell Type'] == groupID].index.values[0]

            sumValue = results.loc[results['Cell Type'] == groupID]["Sum"].values
            totalValue = sumValue + row["Scoring"]

            numValue = results.loc[results['Cell Type'] == groupID]["Num"].values
            results.at[indexValue, "Sum"] = totalValue
            results.at[indexValue, "Num"] = numValue+1
            results.at[indexValue, "Score"] = totalValue/float(numValue+1)
        else:
            results.loc[len(results.index)] = [groupID, row["Scoring"], 1, row["Scoring"]]
    results = results.sort_values(['Sum'], ascending = [False])
    return results
        

# genesResult = findCellTypeIndividualCellTypes(genes)
# genesResult = genesResult.sort_values(['Scoring'], ascending = [False])
# print(genesResult)
# print(genesResult['Scoring'].iloc[0])
# groupResult = findCellTypesGroup(genesResult)
# print(groupResult)
# print(groupResult['Sum'].iloc[0])

# sampleList = ["EN", "EN", "EN", "EN", "IN", "AS", "IN", "EN", "EN", "EN", "EN", "IN", "AS", "IN"]
# # print(sampleList)
# # print(renameList(sampleList))