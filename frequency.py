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

inputData = "./output/savedDataFilterBadScore.h5ad" # "./output/savedDataFilterBadScore.h5ad"
results_file = './output/savedData.h5ad'
outputDirectory = "./outputPDFs/"

adata = sc.read(inputData)
adata.var_names_make_unique()

####################################################################################################
# Frequency Analysis 
####################################################################################################
with PdfPages(outputDirectory+'Frequency Analysis .pdf') as pdf:
    dummyPlot = plt.plot([1, 2])
    figureNum = plt.gcf().number

    num_tot_cells = adata.obs.groupby(['group']).count()
    num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.nFeature_RNA))

    cell_type_counts = adata.obs.groupby(['group', 'sample','leiden']).count()
    print(cell_type_counts)
    cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
    print(cell_type_counts)
    cell_type_counts = cell_type_counts[cell_type_counts.columns[0:6]]
    print(cell_type_counts)

    cell_type_counts['total_cells'] = cell_type_counts.group.map(num_tot_cells).astype(int)
    print(cell_type_counts)

    cell_type_counts['frequency'] = cell_type_counts.nFeature_RNA / cell_type_counts.total_cells
    print(cell_type_counts)

    cell_type_counts = cell_type_counts.sort_values(by=['leiden'])
    
    plt.figure(figsize = (10,4))

    ax = sns.boxplot(data = cell_type_counts, x = 'leiden', y = 'nCount_RNA', hue = 'group')

    plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')


    for fig in range(figureNum+1,  plt.gcf().number+1):
        pdf.savefig(figure=fig, bbox_inches='tight')

