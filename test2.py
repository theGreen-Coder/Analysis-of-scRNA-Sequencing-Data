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

adata = sc.read("./output/savedDataClustersFinal.h5ad")
print(adata)