import pandas as pd
import json
from tqdm.auto import tqdm
import labelCytoscape

networkDS = pd.read_csv('./data/adjDatasetFullGenes.tsv', sep='\t')

genesToLabelled = networkDS["TF"].to_list() + networkDS["target"].to_list()
genesToLabelled = list(dict.fromkeys(genesToLabelled))

labelCytoscape.createLabellingTable(genesToLabelled, "./adjDatasetFullGenes2000Labels.csv")