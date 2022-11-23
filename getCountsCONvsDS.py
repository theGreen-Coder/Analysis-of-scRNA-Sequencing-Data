import pandas as pd
import scanpy as sc

data = pd.read_csv("./outputPDFs/DiffExpOut/matScaledFullDataset.csv")
data = data.set_index(["Unnamed: 0"])

controlData = pd.DataFrame(index=data.index)
downData = pd.DataFrame(index=data.index)

for column in data:
    if(column.startswith("CON")):
        controlData[column] = data[column]
    elif(column.startswith("DS")):
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
        colName = column.split(".")[1]
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
        colName = column.split(".")[1]
        downAvgData[colName] = sampleDF['avg']

        sampleDF = sampleDF.iloc[0:0]
        count=0

# print(controlAvgData) 
# print(downAvgData)

controlAvgData.to_csv("./output/AvgCountsCONTROL.csv")
downAvgData.to_csv("./output/AvgCountsDOWN.csv")

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
diffAvgData.to_csv("./output/DiffCountsCONvsDS.csv")