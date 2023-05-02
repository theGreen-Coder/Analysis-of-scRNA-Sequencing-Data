### analyse differential abundance and gene expression by cell class (DESeq2 pseudobulk analysis)
#       also plot UMAP etc with cluster and cell class labels

message("####### Start R05 (by cluster): ", Sys.time())

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(DESeq2)
library(colorRamps)
library(viridis)
library(lmerTest)
library(pheatmap)
library(clusterProfiler)
# library(DOSE)
library(org.Hs.eg.db)

# define output folder
out_dir = "./output/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

out_dir2 = "./output/"
if (!dir.exists(out_dir2)) {dir.create(out_dir2)}

#load group and file info
gr_tab = read_csv("input/group_tab_rcs.csv")

#get marker gene panels
marker_tab = read_csv("input/cell type markers 22-02-23.csv")
GOI = list()
t1 = read_csv("input/Transcription Factors hg19 - Fantom5_21-12-21.csv")
GOI$TF = t1$Symbol 
t1 = read_csv("input/HSA21_genes_biomaRt_conversion.csv")
GOI$HSA21 = na.omit(t1$hgnc_symbol)

#load processed/integrated dataset
# load(file = "./output/R04_seur_clean_clustered.rda") 
# mSeur = seur 
load(file = "./rawDataset.rda") 
seur$leiden = as.factor(seur$leiden)
seur$seurat_clusters = seur$leiden
seur$cell_class = seur$leiden

# load cluster assignment (clusters not labelled here)
clust_assign_tab = read_csv("./input/R04_cluster_assignment_cleaned.csv")



###########################################################
# plot with updated cluster_name and cell_class (from clust_assign_tab)
###########################################################


#add cluster_name, cell_class
t1 = clust_assign_tab
seur$cluster_name = seur$leiden
seur$cell_class = seur$leiden

# define color scales for groups and samples (avoid white color usage for matlab.like(3/5))
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
gr_colors = matlab.like(6)
if (length(gr)>5){gr_colors = matlab.like(length(gr))} 
sample_colors = matlab.like(6)
if (length(samples)>5){sample_colors = matlab.like(length(samples))}

clusters = unique(clust_assign_tab$cell_class)
cluster_names = unique(clust_assign_tab$cell_class)
clust_colors = matlab.like(length(clusters))

cell_classes = unique(clust_assign_tab$cell_class)
class_colors = matlab.like(6)
if (length(cell_classes)>5){class_colors = matlab.like(length(cell_classes))}

# clusters <- factor(clusters, levels = c(clust_assign_tab$cell_class))
# cluster_names <- factor(cluster_names, levels = c(clust_assign_tab$cell_class))

#############################################################################
# cluster quantification and differential abundance analysis
#############################################################################

meta = seur@meta.data
meta  = meta[order(match(meta$sample, gr_tab$sample)),] #order to match sample order in gr_tab


# create table with one row for each cluster for each sample

cl = unique(clust_assign_tab$cell_class)
samples = unique(meta$sample)

stat_tab = tibble(cluster = unlist(lapply(cl, rep, length.out = length(samples))), 
                  sample = rep(samples, length(cl)))
stat_tab$group = gr_tab$group[match(stat_tab$sample, gr_tab$sample)] #add group assignment

#count cells per cluster per sample, add to stat_tab
t1 = meta %>% group_by(cluster_name, sample) %>% summarize(N_cells = n())
stat_tab$N_cells = t1$N_cells[match( paste0(stat_tab$cluster,stat_tab$sample), paste0(t1$cluster_name,t1$sample) )]
stat_tab$N_cells[is.na(stat_tab$N_cells)] = 0

#add total cells per sample and per cluster, fraction of cluster, fraction of sample
N_sample = stat_tab %>% group_by(sample) %>% summarize(N = sum(N_cells))
stat_tab$N_sample = N_sample$N[match(stat_tab$sample, N_sample$sample)]
stat_tab$fract_sample = stat_tab$N_cells / stat_tab$N_sample
N_cluster = stat_tab %>% group_by(cluster) %>% summarize(N = sum(N_cells))
stat_tab$N_cluster = N_cluster$N[match(stat_tab$cluster, N_cluster$cluster)]
stat_tab$fract_cluster = stat_tab$N_cells / stat_tab$N_cluster

write_csv(stat_tab, file = paste0(out_dir2,"05_cell_abundance_by_sample_cluster.csv"))



#crossbar-dotplot quantification of cluster contribution fraction of sample

t1 = stat_tab
t2 = t1 %>% group_by(cluster, group) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster, y = mean_fract, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(data = t2, aes(x = cluster,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               color = group),
                position = position_dodge(), width = 0.3, lwd = 0.4)+
  geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group, shape = sample), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  geom_hline(yintercept = 0)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_fill_manual(limits = gr, values = gr_colors)+
  scale_shape_manual(limits = unique(stat_tab$sample), 
                     values = rep_len(c(1:4), length.out = length(unique(stat_tab$sample))))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf(file = paste0(out_dir2,"R05_fract_sample_in_cluster_barplot.pdf"), width = 9, height = 4)
plot(p1)
dev.off()



## statistical analysis with GLMM (lmerTest algorithm) approach from by Yu et al., 2022(Neuron) 

t1 = stat_tab
cl = unique(t1$cluster)
samples = unique(t1$sample)
gr = unique(t1$group)


# generate logistic regression models for all clusters in vs not in cluster

lmerTest_list = lapply(cl, function(cl) {
  t2 = meta
  t2$in_clustX = t2$cluster_name == cl
  mod1 = lmerTest::lmer(in_clustX~group+(1|sample), data=t2)
  return(mod1)
})

names(lmerTest_list) = cl


# extract test statistics from models, perform multiple testing correction (BH)

t3 = NULL

for (i in cl){
  t4 = as.data.frame(summary(lmerTest_list[[i]])$coefficients)
  t5 = as_tibble(cbind(comp = rownames(t4), t4))
  t3 = as_tibble(rbind(t3, cbind(cluster = i, t5)))
}

t4 = t3[t3$comp != "(Intercept)",]
t4$padj = p.adjust(t4$`Pr(>|t|)`, method = "BH")
lmerTest_out = t4

write_csv(t4, paste0(out_dir2,"/R05_fract_sample_in_cluster_lmerTest.csv"))



## standard t-tests

t1 = stat_tab

cl = unique(t1$cluster)

t2 = tibble(cluster = cl, p = NA, padj = NA)

for (i in cl){
  t4 = t1[t1$cluster == i,]
  t5 = pairwise.t.test(t4$fract_sample, t4$group, p.adjust.method = "none")
  t2$p[t2$cluster == i] = t5$p.value[1,1]
}

t2$padj = p.adjust(t2$p, method = "BH")

write_csv(t2, paste0(out_dir2,"/R05_fract_sample_in_cluster_t-Test.csv"))






################################################################
# extract pseudobulk count data for DEG analysis
################################################################


###create pseudobulk dataset

create_pseudobulk = function(seur, sample, clust){
  
  #extract cells of sample and closest to stage_center
  
  meta = seur@meta.data
  t1 = meta[meta$sample == sample & meta$cluster_name == clust,]
  bulk_cells = rownames(t1)
  
  ### extract count matrix and SCT normalised matrix for stage cells and create pseudobulk
  
  m1 = seur@assays$RNA@counts
  m2 = m1[,bulk_cells]
  pseudobulk_count = apply(as.matrix(m2), 1, sum)
  
  
  ### create list object with count vector, cell_ids, sample, pseudotime stage 
  
  pseudobulk = list(sample = sample, 
                    clust = clust, 
                    count = pseudobulk_count, 
                    cell_ids = bulk_cells)
  return(pseudobulk)
  
}




### extract pseudobulk data for expression for all samples and clusters

l1 = lapply(cluster_names, function(x){
  
  l2 = lapply(samples, function(y){
    
    l3 = create_pseudobulk(seur = seur, sample = y, 
                           clust = x)
    return(l3)
    
  })
  names(l2) = samples
  message("Processed cluster ", x)
  return(l2)
})
names(l1) = cluster_names

#flatten list for easier data extraction
pseudobulk_list = unlist(l1, recursive = F)

save(pseudobulk_list, file = paste0(out_dir2, "R05_pseudobulk_list.rda"))


### create pseudobulk metadata table

meta_b = tibble(pseudobulk = names(pseudobulk_list), sample = NA, group = NA, clust = NA, N_cells = NA)

for (i in names(pseudobulk_list)){
  l1 = pseudobulk_list[[i]]
  meta_b$sample[meta_b$pseudobulk == i] = l1$sample
  meta_b$clust[meta_b$pseudobulk == i] = l1$clust
  meta_b$N_cells[meta_b$pseudobulk == i] = length(l1$cell_ids)
}

meta_b$group = gr_tab$group[match( meta_b$sample, gr_tab$sample)]
# meta_b = meta_b[order(meta_b$clust, meta_b$group, meta_b$sample),]


### create count matrix from pseudobulk list

#remove pseudobulks with <20 cells and clusters containing less than 2 samples per group
t1 = meta_b[meta_b$N_cells>=1,] 
t2 = t1 %>% group_by(clust, group) %>% summarise(N_samples = n())
t3 = t2 %>% group_by(clust) %>% summarise(N_groups = n(), min_samples = min(N_samples))
t4 = t3[t3$N_groups>=2 & t3$min_samples>=2,]
meta_b_cleaned = t1[t1$clust %in% t4$clust,]
pseudobulk_list_cleaned = pseudobulk_list[names(pseudobulk_list)%in% meta_b_cleaned$pseudobulk] 

# create combined pseudobulk count table and table with average counts/cell in each pseudobulk

m1 = matrix(nrow = length(pseudobulk_list_cleaned[[1]]$count), ncol = length(pseudobulk_list_cleaned),
            dimnames = list(gene = names(pseudobulk_list_cleaned[[1]]$count), meta_b_cleaned$pseudobulk) )
m_fract_exp = m1

for (i in names(pseudobulk_list_cleaned)){
  v1 = pseudobulk_list_cleaned[[i]]$count
  m1[,i] = v1
  m_fract_exp[,i] = v1/meta_b_cleaned$N_cells[meta_b_cleaned$pseudobulk == i]
}

# remove genes with low expression (expressed in <20% of cells in all clusters/samples )
v2 = apply(m_fract_exp, 1, max)
m2 = m1[v2 >= 0.2,]

count_mat_b = m2


### merge pseudobulk expression data in single dataset

bulk_data = list(meta = meta_b_cleaned, count_mat = count_mat_b)

save(bulk_data, file = paste0(out_dir2,"R05_pseudobulk_dataset_cleaned.rda"))





#################################################
# DEG analysis for count data with DESeq2
#################################################

t1 = bulk_data$meta
m1 = bulk_data$count_mat

rownames(t1) = t1$pseudobulk
m2 = m1[,t1$pseudobulk]

# merge cluster+group for DESeq model factor (comp)
t1$comp = paste0(t1$clust,"_", t1$group)

dds = DESeqDataSetFromMatrix(m2, colData = t1, design = ~ comp)
dds = DESeq(dds)

bulk_data$dds = dds


# extract DEseq results stats for DS vs CON for each cluster

clusts = unique(t1$clust)

res_list = lapply(clusts, function(x){
  v1 = unique(t1$comp[t1$clust == x])
  dds_res1 = as.data.frame(results(dds, contrast=c("comp", v1)))
  return(dds_res1)
})
names(res_list) = clusts

bulk_data$deseq_results = res_list

# extract up/down genes for DS vs CON for each cluster

genes_up = lapply(res_list, function(res){
  t1 = res[!is.na(res$padj),]
  t1 = t1[t1$log2FoldChange<=-0.3 & t1$padj <=0.1, ]
  return(rownames(t1))
})
names(genes_up) = paste0(names(genes_up) , "_DS_up")

genes_down = lapply(res_list, function(res){
  t1 = res[!is.na(res$padj),]
  t1 = t1[t1$log2FoldChange>=0.3 & t1$padj <=0.1, ]
  return(rownames(t1))
})
names(genes_down) = paste0(names(genes_down) , "_DS_down")

l1 = c(genes_up, genes_down)
l1 = l1[order(names(l1))]

bulk_data$DEG = l1


#extract mean-centered log2 transformed expr data

ntd = normTransform(dds)
m1 = assay(ntd)
m2 = m1 - apply(m1, 1, mean)
bulk_data$deseq_log_matrix = m2

save(bulk_data, file = paste0(out_dir2,"R05_pseudobulk_dataset_cleaned.rda"))




### create DEG stats and save individual genes (includeing separately TFs and HSA21 genes)

l1 = bulk_data$DEG

l2 = lapply(l1, intersect, GOI$TF)
names(l2) = paste0(names(l1), "_TFs")
l3 = lapply(l1, intersect, GOI$HSA21)
names(l3) = paste0(names(l1), "_HSA21")

l1 = c(l1, l2, l3)

t1 = tibble(cluster = names(l1), N_genes = lengths(l1))

write_csv(t1, file = paste0(out_dir2, "R05_diff_genes_DS_vs_CON_N_genes.csv"))

m1 = matrix(nrow = max(lengths(l1)), ncol = length(l1))
colnames(m1) = names(l1)

for (i in names(l1)){
  v1 = l1[[i]]
  if (length(v1) > 0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""
t1 = as_tibble(m1)


write_csv(t1, file = paste0(out_dir2, "R05_diff_genes_DS_vs_CON_genes.csv"))



#############################################################
### plot pseudobulk expression heatmaps of selected genes
#############################################################
#     mean centered expression by sample and cluster
#         group mean for each cluster
#         relative expression DS vs CON
#         p value, adjusted p value with 0.05 cutoff
#############################################################


### calculate expression values for DESeq expression

t1 = bulk_data$meta

clusts = unique(t1$clust)
gr = unique(t1$group)

l1 = lapply(gr, function(gr){
  m1 =  bulk_data$deseq_log_matrix[,t1$pseudobulk[t1$group == gr]]
  l2 = lapply(clusts, function(x){
    m2 = m1[,intersect(colnames(m1), t1$pseudobulk[t1$clust == x])]
    m2 = apply(m2, 1, mean)
    return(m2)
  })
  names(l2) = clusts
  t3 = as.data.frame(l2)
  return(t3)
})
names(l1) = gr

l1$CON = as.matrix(l1$CON)
l1$DS = as.matrix(l1$DS)
l1$DELTA = l1$DS-l1$CON

bulk_data$deseq_by_group_clust = l1

### CUSTOM!!!!!!!!!!!!

t1 = bulk_data$meta

clusts = unique(t1$clust)
gr = unique(t1$group)
samples = unique(t1$sample)

test1 = lapply(samples, function(samples){
  m1 =  bulk_data$deseq_log_matrix[,t1$pseudobulk[t1$sample == samples]]
  l2 = lapply(samples, function(x){
    m2 = m1[,intersect(colnames(m1), t1$pseudobulk[t1$sample == x])]
    m2 = apply(m2, 1, mean)
    return(m2)
  })
  names(l2) = samples
  t3 = as.data.frame(l2)
  return(t3)
})
names(test1) = samples

test1$CON_DS2U = as.matrix(test1$CON_DS2U)
test1$CON_H9 = as.matrix(test1$CON_H9)
test1$CON_IMR = as.matrix(test1$CON_IMR)
test1$DS_2DS3 = as.matrix(test1$DS_2DS3)
test1$DS_DSP = as.matrix(test1$DS_DSP)

dfSamples <- data.frame(test1$CON_DS2U, test1$CON_H9, test1$CON_IMR, test1$DS_2DS3, test1$DS_DSP)

bulk_data$deseq_by_samples = test1

### create stat matrices for DESeq analysis

# -log10(p) matrix

l1 = lapply(clusts, function(x){
  x = as.character(x)
  v1 = bulk_data$deseq_results[[x]]$pvalue
  names(v1) = rownames(bulk_data$deseq_results[[x]])
  return(v1)
})
names(l1) = clusts

m1 = as.matrix(as.data.frame(l1))

bulk_data$deseq_by_group_clust$log10p = -log10(m1)


# -log10(padj) matrix (cutoff < -log10(0.05) and NAs as 0)

l1 = lapply(clusts, function(x){
  x = as.character(x)
  v1 = bulk_data$deseq_results[[x]]$padj
  names(v1) = rownames(bulk_data$deseq_results[[x]])
  return(v1)
})
names(l1) = clusts

m1 = as.matrix(as.data.frame(l1))
m2 = -log10(m1)
m2[m2 <= -log10(0.05)] = 0
m2[is.na(m2)] = 0

bulk_data$deseq_by_group_clust$log10padj = m2



save(bulk_data, file = paste0(out_dir2,"R05_pseudobulk_dataset_cleaned.rda"))



#################################################
# plots for gene expression by cluster
#################################################


### function to plot nice heatmaps
#     return pheatmap object to allow to reconstruct gene order from clustering

plot_heatmap = function(mat, genes = rownames(mat), cluster_rows = FALSE,
                        color = colorRampPalette(c("blue", "white", "red"))(250),
                        fontsize = 10, cellheight = 10, cellwidth = 10, lim = NULL,
                        annotation_col = NA, main = ""){
  
  m2 = mat[match(genes, rownames(mat), nomatch = 0),]
  
  #define heatmap layout parameters
  
  if (is.null(lim)){
    lim = range(na.omit(m2))
    if (lim[1] == lim[2]){lim = c(0, 1)} #else can't plot matrices with all values equal
  }
  
  #plot 
  
  pl1 = pheatmap(m2, show_rownames=TRUE, cluster_rows = cluster_rows,
                 cluster_cols = FALSE, show_colnames = TRUE, 
                 clustering_distance_rows = "euclidean",
                 clustering_method = "ward.D2",
                 treeheight_row = 50,
                 color = color,
                 breaks = seq(lim[1], lim[2], length.out = length(color)+1),
                 border_color = NA, fontsize = fontsize,
                 cellwidth = cellwidth, cellheight = cellheight,
                 annotation_col = annotation_col,
                 main = main
  )
  
  return(pl1)
  
}


### function to plot gene expression from pseudobulk: 
#       individual samples, control group mean, deviation DS vs CON, -log10(p), -log10(padj)

plot_cluster_diff_expr = function(bulk_data, pl_genes, 
                                  file = "clust_diff_expr.pdf"){
  
  # extract metadata and sct data
  
  meta = bulk_data$meta
  
  m_sample = bulk_data$deseq_log_matrix[intersect(rownames(bulk_data$deseq_log_matrix), pl_genes),]
  l1 = bulk_data$deseq_by_group_clust
  
  # plot settings (plot and cell dimensions, common range for expression values, annotations)
  
  ng = length(pl_genes)
  
  if (ng<30){cellheight = 10}else if (ng<60){cellheight = 5} else {cellheight = 300/ng}
  if (ng<30){fontsize = 10}else if (ng<60){fontsize = 5} else {fontsize = 10}
  
  v1 = max(abs(m_sample))
  lim_expr = c(-v1, v1)*0.8
  
  annotation_col = as.data.frame(meta[,c("group", "clust")])
  rownames(annotation_col) = meta$pseudobulk
  
  pdf(file = file, 
      width = ncol(m_sample)/3+2, height = nrow(m_sample)*cellheight/40+2)
  {
    p1 = plot_heatmap(m_sample, genes = pl_genes, cluster_rows = TRUE,
                      color = viridis_pal(option = "viridis")(250),
                      fontsize = fontsize, cellheight = cellheight, cellwidth = 10, lim = lim_expr,
                      annotation_col = annotation_col,
                      main = "Norm Expr by Sample")
    
    pl_genes_clust =p1$tree_row$labels[p1$tree_row$order]
    write.csv(pl_genes_clust, "./significantGenes-padj<0.1.csv", row.names=FALSE)
    
    plot_heatmap(l1$CON, genes = pl_genes_clust, cluster_rows = FALSE,
                 color = viridis_pal(option = "viridis")(250),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, lim = lim_expr,
                 main = "CON Mean Expr by cluster")
    
    plot_heatmap(l1$DS, genes = pl_genes_clust, cluster_rows = FALSE,
                 color = viridis_pal(option = "viridis")(250),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, lim = lim_expr,
                 main = "DS Mean Expr by cluster")
    
    plotClusterDifference = plot_heatmap(l1$DELTA, genes = pl_genes_clust, cluster_rows = TRUE,
                                         color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                                         fontsize = fontsize, cellheight = cellheight, cellwidth = 10, lim = c(-2,2),
                                         main = "DS vs CON Expr by cluster")
    
    dataframePlot = l1$DELTA
    write.csv(dataframePlot, "./geneOrderDataframe-padj<0.1.csv", row.names=TRUE)
    plotDifferenceOrder = plotClusterDifference$tree_row$labels[plotClusterDifference$tree_row$order]
    # write.csv(plotDifferenceOrder, "./geneOrderDifference.csv", row.names=FALSE)
    
    plot_heatmap(l1$log10p, genes = pl_genes_clust, cluster_rows = FALSE,
                 color = viridis_pal(option = "magma")(250),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, 
                 main = "-log10(p) Expr by cluster")
    
    plot_heatmap(l1$log10padj, genes = pl_genes_clust, cluster_rows = FALSE,
                 color = viridis_pal(option = "magma")(250),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, 
                 main = "-log10(padj) Expr by cluster")
    genesToPlot = c("DYRK1A", "S100B", "OLIG2", "C21orf5", "DOPEY2", "DSCAM", "SYNJ1", "ITSN1", "SOD1", "ERG", "APP", "DSCR1", "RCAN1", "DSCR5", "RIPPLY3", "DSCR6", "SIM2", "DNMT3L", "PKNOX1", "DSCR4", "KNCJ6", "GIKR2", "RUNX1")
    significantHSA21genes = c("DYRK1A", "APP", "RCAN1", "HMGN1", "SYNJ1", "TIAM1", "MCM3AP", "DOP1B", "TTC3", "PRMT2", "GET1", "HSPA13", "RRP1B", "ATP5PF", "CCT8", "ITSN1", "CRYZL1", "PDE9A", "PAXBP1", "MORC3", "HLCS", "LSS", "PKNOX1", "BRWD1", "IFNAR2", "GART", "DONSON", "GATD3A")
    keyDevelopmentalGenes = c("EOMES", "SOX5", "PROX1", "LHX5", "LHX5-AS1", "DLX6", "TNC", "VIM", "SOX4", "TNIK", "LRRC4C", "UNC5C", "SHTN1", "EGFR")
    nonSignificantGenes = c("DLX1", "DLX2", "DLX5", "PAX6", "NEUROD1", "NEUROD2", "NEUROD6", "SOX2", "RBFOX2", "GLI3", "CELF4", "AUTS2", "ZBTB20", "ZNF558")
    top25Up = c('COL24A1','ADAMTSL1','EOMES','KHDRBS3','DLK1','DNAJC15','AL049637.2','CCND1','UNC5C','LEF1','PCDHB5','FAM89A','DMRT3','TMEM132D','DMRTA2','MYO1E','ARHGAP25','LHX5','WNT7B','TTR','LHX5-AS1','LHX1','AP000936.1','ZNF558','CHCHD2')
    top25Down = c('UTY','LINC02506','AL162493.1','GATD3A','AC012645.1','EGFR','LHFPL3','SCRG1','SST','TNC','DLX6','SERTAD4','SH3RF3','GAD1','LINC02223','LRRC4C','AP001347.1','DNM3','SHTN1','PHACTR3','SP9','LSS','MCM3AP-AS1','CD24','CACNG2')
    top50 = c(top25Down, top25Up)
    allGenesToPlot = c(keyDevelopmentalGenes, nonSignificantGenes)
    plot_heatmap(dfSamples, genes = significantHSA21genes, cluster_rows = FALSE,
                 color = viridis_pal(option = "viridis")(250),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, lim = c(-0.6,0.6),
                 main = "K")
    plot_heatmap(dfSamples, genes = significantHSA21genes, cluster_rows = FALSE,
                 color = colorRampPalette(c("magenta", "black", "yellow"))(250), lim = c(-0.6,0.6),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, 
                 main = "K")
    plot_heatmap(dfSamples, genes = keyDevelopmentalGenes, cluster_rows = TRUE,
                 color = colorRampPalette(c("magenta", "black", "yellow"))(250), lim = c(-0.6,0.6),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, 
                 main = "K")
    plot_heatmap(dfSamples, genes = allGenesToPlot, cluster_rows = TRUE,
                 color = colorRampPalette(c("magenta", "black", "yellow"))(250), lim = c(-0.6,0.6),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 10, 
                 main = "K")
    plot_heatmap(dfSamples, genes = top50, cluster_rows = FALSE,
                 color = colorRampPalette(c("magenta", "black", "yellow"))(250), lim = c(-1.5, 1.5),
                 fontsize = fontsize, cellheight = cellheight, cellwidth = 20, 
                 main = "K")
    
  }
  
  dev.off()
  
}



### load cleaned pseudobulk data with stats

load(file = paste0(out_dir2,"R05_pseudobulk_dataset_cleaned.rda"))

### plot selected gene sets
# all differentially expressed genes (and top50 by p value)
# all and all diff HSA21 genes
# all diff TFs (and top30 by p value)

t1 = bulk_data$DEG

deg = unique(unlist(t1))
plot_cluster_diff_expr(bulk_data, deg,
                       file = paste0(out_dir2,"05_heatmaps_all_diff_genes_16Apr.pdf"))


v1 = intersect(deg, GOI$HSA21)
plot_cluster_diff_expr(bulk_data, v1,
file = paste0(out_dir2,"05_heatmaps_HSA21_diff_genesHeatmapFigures_16Apr.pdf"))


# v1 = GOI$HSA21
# plot_cluster_diff_expr(bulk_data, v1,
#                      file = paste0(out_dir2,"05_heatmaps_HSA21_genesTest.pdf"))
# 
# v1 = intersect(deg, GOI$TF)
# plot_cluster_diff_expr(bulk_data, v1,
#                        file = paste0(out_dir2,"05_heatmaps_diff_TFs.pdf"))
# 
# # extract genes with top p values
# t1 = bulk_data$deseq_by_group_clust$log10p
# t2 = apply(t1,1,max)
# max_log10p = t2[order(-t2)]
# 
# v1 = names(max_log10p)[1:50]
# plot_cluster_diff_expr(bulk_data, v1,
#                        file = paste0(out_dir2,"05_heatmaps_diff_genes_top50.pdf"))
# 
# v1 = intersect(names(max_log10p), GOI$TF)[1:30]
# plot_cluster_diff_expr(bulk_data, v1,
#                        file = paste0(out_dir2,"05_heatmaps_diff_TFs_top30.pdf"))


#################################################
### GO analysis of DEG
#################################################

GO_dir = paste0(out_dir2,"/05_GO_DEG/")
if (!dir.exists(GO_dir)) {dir.create(GO_dir)}

#get DEG
GO_genes_list = bulk_data$DEG

comp = names(GO_genes_list)


#perform GO analysis

ego_list = lapply(comp, function(comp){
  ego = enrichGO(gene         = GO_genes_list[[comp]],
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
  
  message(paste0("completed GO analysis ", comp))
  
  return(ego)
  
})

names(ego_list) = comp

save(ego_list, file = paste0(out_dir2,"05_DEG_GO_analysis.rda")) 


#filter and save GO results as tables (only tables with sign GO terms)

t1 = lapply(comp, function(comp){
  ego = ego_list[[comp]]
  
  if (!is.null(ego)){
    # ego = dropGO(ego, level = c(1,2))
    t1 = ego@result
    t1 = t1[t1$Count > 2 & t1$p.adjust <= 0.05,]
    
    if (nrow(t1)>0){
      write_csv(t1, file = paste0(GO_dir,"/GO_DEG_",comp,".csv"))
    }
    
  }
  
  
})



message("####### End R05 (by cluster): ", Sys.time())