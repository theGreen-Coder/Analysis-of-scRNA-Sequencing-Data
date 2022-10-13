library(scater)
library(Seurat)
library(cowplot)
library(reticulate)

ad <- import("anndata", convert = FALSE)
pbmc_ad <- ad$read_h5ad("./R_in/07_seurat_1000_cells.rda")
pbmc3k <- Convert(pbmc_ad, to = "seurat")
p1 <- TSNEPlot(pbmc3k, group.by = "louvain", do.return = TRUE)
p2 <- VlnPlot(pbmc3k, c("CST3", "NKG7", "PPBP"), group.by = "louvain", do.return = TRUE)
plot_grid(p1, p2)

pbmc_ad <- Convert(from = pbmc, to = "anndata", filename = "./pbmc3k.h5ad")
pbmc_ad