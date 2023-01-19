
import hotspot
from itertools import groupby
from json import load
import pandas as pd
import scanpy as sc
import numpy as np
import sys
sys.path.insert(1, '/Users/greencode/Documents/Coding/Analysis-of-scRNA-Sequencing-Data')
import modules.GO.geneOntology as GOAnalysis
import os
import math
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pickle
import random
import gseapy as gp
import json
from matplotlib.backends.backend_pdf import PdfPages
pd.set_option('display.max_columns', None)

minGenes = 10
outputDirectory = "./results/"+str(minGenes)+"genes/"

MYDIR = ("./results/"+str(minGenes)+"genes/")
CHECK_FOLDER = os.path.isdir(MYDIR)

# If folder doesn't exist, then create it.
if not CHECK_FOLDER:
    os.makedirs(MYDIR)
    print("created folder : ", MYDIR)

else:
    print(MYDIR, "folder already exists.")

significantGenes = list(pd.read_csv("./significantGenes.csv")["genes"])
HSA21genesDataframe = pd.read_csv("./HSA21_genes_biomaRt_conversion.csv")
HSA21genes = [x for x in HSA21genesDataframe["hgnc_symbol"] if str(x) != 'nan']

genesInSignificant = 0
significantHSA21genes = []

for gene in significantGenes:
    if gene in HSA21genes:
        genesInSignificant +=1
        significantHSA21genes.append(gene)

print(genesInSignificant)

file = open(outputDirectory+'significantHSA21genes.txt','w')
for item in significantHSA21genes:
	file.write(item+"\n")
file.close()


with open('./hotspotObjectSignificantGenes.pkl', 'rb') as inp:
    hs = pickle.load(inp)

hs.compute_local_correlations(significantGenes)

# Creating Modules
modules = hs.create_modules(
    min_gene_threshold=minGenes, core_only=True, fdr_threshold=0.05
)
# Plot Local Correlations
plotLocalCorrelations = hs.plot_local_correlations()
plt.savefig(outputDirectory+'LocalCorrelationsPlot.png')

modules = modules[~modules.index.duplicated(keep="first")]
modulesGeneDict = modules.to_dict()

modulesSignificantDictionary = {}
for i in range(1, modules.max()+1):
    modulesSignificantDictionary[i] = []
modulesSignificantDictionary

for gene in significantGenes:
    moduleGene = modulesGeneDict[gene]
    if moduleGene != -1:
        modulesSignificantDictionary[moduleGene].append(gene)
modulesSignificantDictionary

with open(outputDirectory+"modulesSignificantGenesDictionary.json", "w") as outfile:
    json.dump(modulesSignificantDictionary, outfile, indent=4)

modulesHSA21Dictionary = {}
for i in range(1, modules.max()+1):
    modulesHSA21Dictionary[i] = []
modulesHSA21Dictionary

for gene in significantHSA21genes:
    moduleGene = modulesGeneDict[gene]
    if moduleGene != -1:
        modulesHSA21Dictionary[moduleGene].append(gene)
modulesHSA21Dictionary

with open(outputDirectory+"modulesHSA21Dictionary.json", "w") as outfile:
    json.dump(modulesHSA21Dictionary, outfile, indent=4)

hs.results = hs.results.loc[list(modules.index)]

# Convert Modules into Pandas Dataframe
modulesDataframe = pd.DataFrame(modules)
modulesDataframe["Z"] = hs.results["Z"].tolist()
modulesDataframe = modulesDataframe[modulesDataframe["Module"] != -1]
maximumModule = max(list(modulesDataframe["Module"]))

# Generate Gene Set Enrichment Analysis
for i in range(1, maximumModule+1):
    moduleRank = modulesDataframe[modulesDataframe["Module"] == i]
    moduleRank = moduleRank.rename_axis("Gene").reset_index()
    ranking = moduleRank[['Gene', 'Z']]
    rankingGenes = list(moduleRank["Gene"])

    results = GOAnalysis.geneOntologyAnalysis(rankingGenes)
    results.to_csv(outputDirectory+"GSEA_Module_"+str(i)+".csv")
    print(results)

    # pre_res = gp.prerank(rnk = ranking, gene_sets = 'GO_Biological_Process_2021', min_size=10, seed = 6, permutation_num = 100)
    # out = []
    # for term in list(pre_res.results):
    #     out.append([term,
    #             pre_res.results[term]['fdr'],
    #             pre_res.results[term]['es'],
    #             pre_res.results[term]['nes']])

    # out_df = pd.DataFrame(out, columns = ['Term','fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop = True)
    # out_df = out_df[out_df["fdr"] < 0.05]
    # out_df.to_csv("./output/GSEA_Module_"+str(i)+".csv")

# Set List of Prefered GO Terms
selectedTerms = [0,0,0,0,0,0,0]


# Calculating Module Scores
module_scores = hs.calculate_module_scores()

def pp(adataFile):
    inputDataAdata = adataFile
    adata = sc.read(inputDataAdata)
    adata.var_names_make_unique()
    adata.X = adata.X.astype('float64')
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    riboURL = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    riboGenes = pd.read_table(riboURL, skiprows=2, header = None)
    adata.var['ribo'] = adata.var_names.isin(riboGenes[0].values)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40) # n_neighbors=15 is default
    sc.tl.umap(adata)
    sc.external.pp.bbknn(adata, batch_key='sample')
    # Plotting UMAP after integration
    sc.tl.umap(adata)
    return adata
adata = pp("../output/savedDataClustersFinal.h5ad")

modulesDictionary = {}
for i in range(1, maximumModule+1):
    modulesDictionary[i] = 'module'+str(i)
modulesDictionary


if(list(adata.obs.index) == list(module_scores.index)):
    pdConcat = pd.concat([adata.obs,module_scores], axis=1, ignore_index=False)

pdConcat = pdConcat.rename(columns=modulesDictionary)

adata.obs = pdConcat

with PdfPages(outputDirectory+'GO Analysis of Modules.pdf') as pdf:
    # I have too include this dummy plot everytime I want to export something into PDF. 
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    sc.pl.umap(adata, color=list(modulesDictionary.values()), show=False)

    adataCON = adata[adata.obs["group"] == "CON"]
    sc.pl.umap(adataCON, color=list(modulesDictionary.values()), show=False)


    adataDS = adata[adata.obs["group"] == "DS"]
    sc.pl.umap(adataDS, color=list(modulesDictionary.values()), show=False)


    goRenameDictionary = {}
    for i in range(1, maximumModule+1):
        fileName = pd.read_csv(outputDirectory+"GSEA_Module_"+str(i)+".csv")
        if(len(fileName.index)!=0):
            goTerm = fileName["term"][selectedTerms[i-1]]
        else:
            goTerm = "NO SIGNIFICANT GO TERM"+str(i)
        goRenameDictionary['module'+str(i)] = goTerm

    renamedAdata = adata
    renamedAdata.obs = renamedAdata.obs.rename(columns=goRenameDictionary)

    list(goRenameDictionary.values())

    sc.pl.umap(renamedAdata, color=list(goRenameDictionary.values()), show=False)

    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')




