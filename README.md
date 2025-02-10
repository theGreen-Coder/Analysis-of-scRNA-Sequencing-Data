# Analysis of scRNA Sequencing Data

This repository contains the code used for analyzing single-cell RNA sequencing (scRNA-seq) data to investigate neurodevelopmental defects in Down syndrome (DS). The analysis is based on data from Tang et al. and focuses on identifying key altered neurodevelopmental genes and regulatory interactions.

## Overview

The analysis follows a structured pipeline:

1. **Data Preprocessing**: 
   - Raw scRNA-seq data is processed using the **10X Genomics Cell Ranger pipeline**.
   - Quality control is performed using **Scanpy** to remove low-quality cells.
   - Data normalization and log transformation are applied.
   - Highly variable genes are selected for further analysis.

2. **Cell Clustering and Labeling**:
   - **UMAP** is used for dimensionality reduction.
   - **Batch-balanced KNN (BBKNN)** is used to reduce batch effects.
   - Cells are clustered using the **Leiden algorithm**.
   - Clusters are labeled using a **custom algorithm** based on gene expression comparison.

3. **Differential Expression Analysis**:
   - **DESeq2** is used to identify differentially expressed genes (DEGs) between DS and control organoids.
   - Significant genes are mapped to their corresponding chromosomes.
   - Gene ontology (GO) analysis is conducted to identify enriched biological processes.

4. **Regulatory Network Inference**:
   - **pySCENIC** is used to infer gene regulatory networks.
   - Motif enrichment analysis is performed to identify key transcription factor (TF) targets.
   - **Cytoscape** is used for visualization of regulatory interactions.

5. **Post-Translational Interaction Analysis**:
   - **OmniPath** is used to explore indirect interactions between HSA21 genes and neurodevelopmental regulators.
   - NetworkX is used to find the shortest interaction paths.

## Requirements

The analysis requires the following dependencies:

- Python (>=3.8)
- R (>=4.0)
- Scanpy (`pip install scanpy`)
- DESeq2 (R package)
- pySCENIC (`pip install pyscenic`)
- Cytoscape (for visualization)
- OmniPath (`pip install omnipath`)
- Other standard scientific Python libraries (NumPy, Pandas, Matplotlib, Seaborn)

## Usage

This repository was primarily created for my own analysis rather than for public use, so the code might be somewhat convoluted and not fully optimized for ease of execution. Some modifications may be required to adapt it to different datasets or environments.

## Citation
If you use this code, please cite the original study:

- Tang et al. (2021) - [[Link to dataset](https://www.jci.org/articles/view/135763)]
- This repository - [[GitHub link](https://github.com/theGreen-Coder/Analysis-of-scRNA-Sequencing-Data)]
