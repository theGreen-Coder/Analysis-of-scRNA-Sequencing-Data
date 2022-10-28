import pandas as pd
import numpy as np
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
pd.options.display.max_columns = sys.maxsize

# Genes For Cluster 9
# genes = ["CENPF", "NUSAP1", "HMGB2", "PTTG1", "TOP2A", "TPX2", "VIM", "CCNB1", "DLGAP5", "SMC4", "H2AFX", "HMGN2", "CCNA2", "CKS2", "MKI67", "BIRC5", "ASPM", "CENPE", "PBK", "NUF2", "GTSE1", "PRC1", "SGO2", "CKAP2", "CCNB2"]

def findCellType(listGenes):
    results = pd.DataFrame(data={'Cell Type': [], 'Total Found Genes': [], 'Total Average Difference Genes': []})
    row = 0
    for filename in os.listdir("./clusterGeneNames/files/"):
        df = pd.read_csv('./clusterGeneNames/files/'+filename, sep='\t', header=0)
        i = 0
        avgDiff = 0
        for gene in listGenes:
            if gene in df["id"].to_list():
                rowDataframe = df.loc[df['id'] == gene]
                avgDiff = float(rowDataframe["avg_diff|float"])
                i+=1
        results.loc[row] = [filename.replace('.tsv', ''), i , avgDiff]
        row +=1

    results = results.sort_values(['Total Found Genes', 'Total Average Difference Genes'], ascending = [False, False])
    return results
