import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib.pyplot import rc_context
from goatools.cli.ncbi_gene_results_to_python import ncbi_tsv_to_py

from genes_ncbi_9606_proteincoding import GENEID2NT as GeneID2nt_mus
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap

ncbi_tsv = 'gene_result.txt'
output_py = 'genes_ncbi_9606_proteincoding.py'

def geneOntologyAnalysis(listOfGenes):
    # run one time to initialize
    # obo_fname = download_go_basic_obo()
    fin_gene2go = download_ncbi_associations()
    obodag = GODag("go-basic.obo")

    #run one time to initialize
    mapper = {}

    for key in GeneID2nt_mus:
        mapper[GeneID2nt_mus[key].Symbol] = GeneID2nt_mus[key].GeneID
        
    inv_map = {v: k for k, v in mapper.items()}

    #run one time to initialize

    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, taxids=[9606])
    # Get namespace2association where:
    #    namespace is:
    #        BP: biological_process               
    #        MF: molecular_function
    #        CC: cellular_component
    #    assocation is a dict:
    #        key: NCBI GeneID
    #        value: A set of GO IDs associated with that gene
    ns2assoc = objanno.get_ns2assc()

    #run one time to initialize
    goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_mus.keys(), # List of mouse protein-coding genes
            ns2assoc, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method

    #run one time to initialize
    GO_items = []

    temp = goeaobj.ns2objgoea['BP'].assoc
    for item in temp:
        GO_items += temp[item]
        

    temp = goeaobj.ns2objgoea['CC'].assoc
    for item in temp:
        GO_items += temp[item]
        

    temp = goeaobj.ns2objgoea['MF'].assoc
    for item in temp:
        GO_items += temp[item]

    #pass list of gene symbols
    def go_it(test_genes):
        print(f'input genes: {len(test_genes)}')
        
        mapped_genes = []
        for gene in test_genes:
            try:
                mapped_genes.append(mapper[gene])
            except:
                pass
        print(f'mapped genes: {len(mapped_genes)}')
        print(mapped_genes)
        
        goea_results_all = goeaobj.run_study(mapped_genes)
        goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
        GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh,\
                    x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO), list(map(lambda y: inv_map[y], x.study_items)),\
                    ], goea_results_sig)), columns = ['GO', 'term', 'class', 'p', 'p_corr', 'n_genes',\
                                                        'n_study', 'n_go', 'study_genes'])

        GO = GO[GO.n_genes > 1]
        return GO

    # genesLists = ["UTY", "ZNF736", "AC109439.2", "LINC02506", "AC097654.1", "AL162493.1", "AC005863.1", "AL009179.2", "MMP20", "LINCO2511", "LINCO1060", "NPY", "HCG24", "AC090138.1", "AC012645.1", "IL33", "GATD3A", "BCAN", "AL691447.2", "LINC01717", "C80rf34", "DEC1", "LINC02712", "AC016590.", "TNC", "AL355338.1", "NHLH2", "EBF2", "AL049637.2", "WNT7B", "NPIPA8", "ACO25508.1", "LHX1-DT", "RSP01", "VSX1", "CDC20B", "RSP03", "AC116609.3", "b39 10", "LHX5", "C6orf141", "LHX5-AS1", "ILHX1", "FAM135B", "WNT8B", "GMNC", "AC090348.1", "AL354863.1", "COL22A1", "CCNO", "CHCHD2"]
    genesLists = ["CHCHD2", "ZNF736", "AC109439.2", "LINC02506", "AC097654.1", "AL162493.1", "AC005863.1",]

    df = go_it(listOfGenes) #listOfGenes
    print(df)

    return df

    # df['per'] = df.n_genes/df.n_go

    # df = df[0:10]
    # print(df)

    # fig, ax = plt.subplots(figsize = (0.5, 2.75))

    # cmap = mpl.cm.bwr_r
    # norm = mpl.colors.Normalize(vmin = df.p_corr.min(), vmax = df.p_corr.max())

    # mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)

    # cbl = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
