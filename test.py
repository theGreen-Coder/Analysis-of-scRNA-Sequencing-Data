from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
import tests.loadScanpy as loadScanpy
import modules.classifyClusters.classifyClusters as classify
import os
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import diffxpy.api as de
from modules.diffexpr.diffexpr.py_deseq import py_DESeq2
import pandas as pd 
import numpy as np


####################################################################################################
# Global Settings 
####################################################################################################
np.set_printoptions(threshold=sys.maxsize)
pd.options.display.max_columns = sys.maxsize
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

inputData = "./dataSaveOriginal/rawDataset5000.h5ad" # ./dataSaveOriginal/rawDataset5000.h5ad
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

genesHSA21 = pd.read_csv("./dataSaveOriginal/HSA21genes.csv")
genesHSA21List = []

for gene in genesHSA21["hgnc_symbol"].tolist():
    if type(gene) == str:
        genesHSA21List.append(gene)

# ####################################################################################################
# # Pre-Processing
# ####################################################################################################
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
adata.raw = adata
# adata.var['mt'] = adata.var_names.str.startswith('MT-')
# riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
# riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
# adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)
# sc.pp.normalize_total(adata, target_sum=10000)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# sc.pp.scale(adata, max_value=10)
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
# sc.tl.umap(adata)
# sc.tl.leiden(adata, resolution=0.5)
# sc.external.pp.bbknn(adata, batch_key='sample')

# ####################################################################################################
# # Differential Expression Using DESeq2
# ####################################################################################################
print(bcolors.FAIL + "Performing Differential Expression using DESeq2..." + bcolors.ENDC)
def getColSample(string):
    subset =  adata[adata.obs["sample"] == string]
    returnArray = [sum(x) for x in zip(*subset.X)]
    print(string)
    print(returnArray[0])
    return returnArray[0].toarray()[0]

rows = adata.var["features"].values

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP

df = pd.DataFrame(rows, columns=['id'])
df["N_1"] = getColSample("CON_DS2U")
df["N_2"] = getColSample("CON_H9")
df["N_3"] = getColSample("CON_IMR")
df["S_1"] = getColSample("DS_2DS3")
df["S_2"] = getColSample("DS_DSP")

print(df.head(10))

df.N_1 = df.N_1.astype(int)
df.N_2 = df.N_2.astype(int)
df.N_3 = df.N_3.astype(int)
df.S_1 = df.S_1.astype(int)
df.S_2 = df.S_2.astype(int)
print(df.head(50))

print("Setting Sample Dataframe...")
sample_df = pd.DataFrame({'samplename': df.columns}) \
        .query('samplename != "id"')\
        .assign(sample = lambda d: d.samplename.str.extract('([NS])_', expand=False)) \
        .assign(replicate = lambda d: d.samplename.str.extract('_([123])', expand=False)) 
sample_df.index = sample_df.samplename

print("Running DESeq2...")
dds = py_DESeq2(count_matrix = df,
               design_matrix = sample_df,
               design_formula = '~ replicate + sample',
               gene_column = 'id') # <- telling DESeq2 this should be the gene ID column
    
dds.run_deseq() 
dds.get_deseq_result(contrast = ['sample','S','N'])
res = dds.deseq_result 

res = res.sort_values(['log2FoldChange'], ascending = [True])
res = res[res["log2FoldChange"] < 0]
print(res.head(50))

geneList = res["id"].tolist()[0:50]
print(geneList)

# geneList = ["LOXL3", "M1AP"]

plotDF = pd.DataFrame(columns=["CON", "DS"]) # index = geneList
i = 0
for resultGene in geneList:
    rowResult = df.loc[df['id'] == resultGene]
    averageCON = float((rowResult["N_1"]+rowResult["N_2"]+rowResult["N_3"])/3)
    averageDS = float((rowResult["S_1"]+rowResult["S_2"])/2)
    sumTotal = averageCON + averageDS
    normalCON = float(averageCON/sumTotal)
    normalDS = float(averageDS/sumTotal)
    plotDF.loc[resultGene] = [averageCON, averageDS]
    i +=1

genesHSA21data = pd.DataFrame()
s = 0
for geneHSA in genesHSA21List:
    rowResultHSA = res.loc[res['id'] == geneHSA]
    genesHSA21data.loc[geneHSA] = rowResultHSA
    s+=1

print("genesHSA21data")

print("asdfasdf")
print(plotDF)
print(df.loc[df['id'] == "NLGN4Y"]) 

presMarkers = ["MITF", "EMX2", "ZNF558", "KLF3", "ZNF589", "ZNF337", "ZNF16", "HMGA1", "ZIC2", "MECOM", "ESR2", "ZNF536", "RFXANK", "DLX1", "DLX2", "RORB", "OLIG2", "ZBTB20", "CREB3L2", "NPAS2", "ZNF491", "PBX3", "AFF2", "ZNF425", "ZHX3", "PKNOX1", "ARNT", "FOXN3", "GATAD2B", "ASH1L"]

with PdfPages(outputDirectory+'Differential Expression DESeq2 Test.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    # sns.heatmap(plotDF, annot=True, cmap=sns.blend_palette(["black", "green"], as_cmap=True))
    # plt.show()

    genesHSA21List.remove('H2BS1')
    genesHSA21List.remove('CFAP298-TCP10L')
    genesHSA21List.remove('GET1-SH3BGR')
    genesHSA21List.remove('SLX9')


    sc.pl.heatmap(adata, geneList, groupby="group", vmax=5, show=False)
    sc.pl.heatmap(adata, geneList, groupby="sample", vmax=5, show=False)
    sc.pl.heatmap(adata, genesHSA21List, groupby="group", vmax=5, show=False)
    sc.pl.heatmap(adata, genesHSA21List, groupby="sample", vmax=5, show=False)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')



# # from the last cell, we see the arrangement of coefficients, 
# # so that we can now use "coef" for lfcShrink
# # the comparison we want to focus on is 'sample_B_vs_A', so coef = 4 will be used
# lfc_res = dds.lfcShrink(coef=4, method='apeglm')
# lfc_res.head()













# ####################################################################################################
# # Differential Expression Using DIFFXpy
# ####################################################################################################
# print(bcolors.FAIL + "Performing Differential Expression using DIFFXpy..." + bcolors.ENDC)
# # # print(adata.X)
# # # print("Hello there!")
# # # print(adata.X.toarray())
# # # print("Hello there 2!")

# adata.X = csr_matrix(adata.X)

# adata.X = adata.X.toarray()

# # adata = adata.raw.to_adata()
# print(adata.obs)

# adata.obs = adata.obs.rename(columns = {'group':'condition'})
# print(adata.obs)

# print("Starting Wald Test...")
# test = de.test.wald(
#     data=adata.copy,
#     formula_loc="~ 1 + condition",
#     factor_loc_totest="condition"
# )


# dedf = res.summary().sort_values('log2fc', ascending = False).reset_index(drop = True)
# print(dedf)
# adata.obs.cell_type.unique()
# print(adata.obs.cell_type.unique())

# dedf['log2fc'] = dedf['log2fc']*-1
# dedf = dedf.sort_values('log2fc', ascending = False).reset_index(drop = True)
# print(dedf)

# dedf = dedf[(dedf.qval < 0.05) & (abs(dedf.log2fc) > .5)]
# print(dedf)

# dedf = dedf[dedf['mean'] > 0.15]
# print(dedf)

# genes_to_show = dedf[-25:].gene.tolist() + dedf[:25].gene.tolist() #top 25 and bottom 25 from sorted df
# sc.pl.heatmap(adata, genes_to_show, groupby='cell_type', swap_axes=True)