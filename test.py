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

# Arguments Settings

sysArgs = ["-filter"]

showPlots = True

arg = ""
for word in sysArgs:
    if word != "main.py" and word != "-showPlots":
        arg = arg + word + " "

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
adata = sc.read("./dataSaveOriginal/rawDataset.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

print(adata)

geneMarkers = pd.read_csv('./dataSaveOriginal/cellTypeMarkers.csv')
geneNames = geneMarkers["gene"].tolist()

# Remove Outlier Samples from the Data
adataCON_DS2U =  adata[adata.obs["sample"] == "CON_DS2U"]
adataCON_H9 =  adata[adata.obs["sample"] == "CON_H9"]
adataCON_IMR =  adata[adata.obs["sample"] == "CON_IMR"]
adataCON_ihtc =  adata[adata.obs["sample"] == "CON_ihtc"]

adataDS_2DS3 =  adata[adata.obs["sample"] == "DS_2DS3"]
adataDS_DS1 =  adata[adata.obs["sample"] == "DS_DS1"]
adataDS_DSP =  adata[adata.obs["sample"] == "DS_DSP"]

if "-fullDataset" not in sysArgs:
    adata = adataCON_DS2U.concatenate(adataCON_H9)
    adata = adata.concatenate(adataCON_IMR)
    adata = adata.concatenate(adataDS_2DS3)
    adata = adata.concatenate(adataDS_DSP)
else:
    adata = adataCON_DS2U.concatenate(adataCON_H9)
    adata = adata.concatenate(adataCON_IMR)
    adata = adata.concatenate(adataCON_ihtc)
    adata = adata.concatenate(adataDS_2DS3)
    adata = adata.concatenate(adataDS_DSP)
    adata = adata.concatenate(adataDS_DS1)

# CON_DS2U CON_H9 CON_IMR CON_ihtc DS_2DS3 DS_DS1 DS_DSP

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

print(adata.obs)
print(adata.var)

####################################################################################################
# QC Filtering 
####################################################################################################
print(bcolors.FAIL + "Applying Quality Control Filters..." + bcolors.ENDC)
if "-filter" in sysArgs:
    print("Filtering Damaged Cells")
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    adata = adata[adata.obs.n_genes_by_counts > 1000, :]
else:
    print("Not Filtering Damaged Cells")

totalGenes = adata.var.index.to_list()

####################################################################################################
# Clustering 
####################################################################################################
print(bcolors.FAIL + "Performing Clustering..." + bcolors.ENDC)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

neuroDevelopmentalGenes = [
    "ID1",
    "ID2",
    "ID4",
    "ID3",
    "HES1",
    "HES2",
    "HES6",
    "HES3",
    "HES5",
    "HEY1",
    "HEY2",
    "TCF3",
    "TCF12",
    "TCF4",
    "ATOH7",
    "ATOH1",
    "NEUROD1",
    "NEUROD2",
    "NEUROD6",
    "NEUROD4",
    "NEUROG1",
    "NEUROG3",
    "NEUROG2",
    "OLIG2",
    "OLIG3",
    "BHLHE22",
    "BHLHE23",
    "OLIG1",
    "FER1",
    "PTF1A",
    "FER2",
    "FER3",
    "FERD3L",
    "TWIST1",
    "TWIST2",
    "ASCL3",
    "ASCL4",
    "ASCL1",
    "ASCL2",
    "NHLH1",
    "NHLH2",
    "DLX1",
    "DLX2",
    "DLX3",
    "DLX4",
    "DLX5",
    "DLX6",
    "DLX7",
    "DLX8",
    "DLX9",
    "SOX1",
    "SOX2",
    "SOX3",
    "SOX14",
    "SOX21",
    "SOX4",
    "SOX11",
    "SOX12",
    "SOX5",
    "SOX6",
    "SOX13",
    "SOX8",
    "SOX9",
    "SOX10",
    "SOX7",
    "SOX17",
    "SOX18",
    "SOX15",
    "SOX30",
    "NFIA",
    "NFIB",
    "NFIC",
    "NFIX",
    "PAX1",
    "PAX2",
    "PAX3",
    "PAX4",
    "PAX5",
    "PAX6",
    "PAX7",
    "PAX8",
    "PAX9",
    "LHX1",
    "LHX2",
    "LHX3",
    "LHX6",
]
plotDevGenes = []
adata.var
for devGene in neuroDevelopmentalGenes:
    if devGene in totalGenes:
        try:
            adata.var.at[devGene, "highly_variable"] = True
            plotDevGenes.append(devGene)
            print("Added Gene")
        except:
            print("Skipped Gene!")

# adata.var.at["OLIG2", "highly_variable"] = True

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.3)
# adata.write(results_file)

adataDS =  adata[adata.obs["group"] == "DS"]
adataCON =  adata[adata.obs["group"] == "CON"]

####################################################################################################
# Integration
####################################################################################################
print(bcolors.FAIL + "Integrating Datasets..." + bcolors.ENDC)
if "-noIntegration" not in sysArgs and "-ingestIntegration" not in sysArgs:
    print("Performing BBKNN Integration...")
    with PdfPages(outputDirectory+'Integration Plots '+arg+'.pdf') as pdf:
        dummyPlot = plt.plot([1, 2])
        figureNum = plt.gcf().number

        # Plotting UMAP before integration
        sc.pl.umap(adata, color=['sample', "group"], show=showPlots)
        plt.suptitle("Non-Integrated Dataset UMAP", y=1, fontsize=25)
        # BBKNN Integration
        sc.external.pp.bbknn(adata, batch_key='sample')
        # Plotting UMAP after integration
        sc.tl.umap(adata)
        plt.suptitle("Integrated Dataset UMAP", y=1, fontsize=25)
        sc.pl.umap(adata, color='sample', show=showPlots)
        plt.suptitle("Integrated Dataset by Sample", y=1, fontsize=25)
        sc.pl.umap(adata, color='group', show=showPlots)
        plt.suptitle("Integrated Dataset by Group", y=1, fontsize=25)

        for fig in range(figureNum+1,  plt.gcf().number+1):
            pdf.savefig(figure=fig, bbox_inches='tight')

adataDS =  adata[adata.obs["group"] == "DS"]
adataCON =  adata[adata.obs["group"] == "CON"]

####################################################################################################
# Clustering Plots
####################################################################################################
print(bcolors.FAIL + "Saving Clustering Plots ADATA..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Clustering Plots Dev Genes '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.umap(adata, color=plotDevGenes, cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
    plt.suptitle("Developmental Genes", y=1, fontsize=25)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

print(bcolors.FAIL + "Saving Clustering Plots CON..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Clustering Plots Dev Genes CON '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.umap(adataCON, color=plotDevGenes, cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
    plt.suptitle("Developmental Genes CON", y=1, fontsize=25)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

print(bcolors.FAIL + "Saving Clustering Plots DS..." + bcolors.ENDC)
with PdfPages(outputDirectory+'Clustering Plots Dev Genes DS '+arg+'.pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.umap(adataDS, color=plotDevGenes, cmap=sns.blend_palette(["lightgray", "green"], as_cmap=True), show=False)
    plt.suptitle("Developmental Genes DS", y=1, fontsize=25)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')
