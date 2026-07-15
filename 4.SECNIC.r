################################################# Step1 ####################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)

data <- readRDS(file = "Primary_clonocluster.rds")
write.csv(t(as.matrix(data@assays$RNA$counts)),file = "/Volumes/herlynm/linux/ychen/Haiyin/metastasis/SECNIC/Tumor/sce_exp_tumor.csv")

data <- readRDS(file = "Lung_clonocluster.rds")
write.csv(t(as.matrix(data@assays$RNA$counts)),file = "/Volumes/herlynm/linux/ychen/Haiyin/metastasis/SECNIC/Lung/sce_exp_lung.csv")

data <- readRDS(file = "Liver_clonocluster.rds")
write.csv(t(as.matrix(data@assays$RNA$counts)),file = "/Volumes/herlynm/linux/ychen/Haiyin/metastasis/SECNIC/Liver/sce_exp_liver.csv")

################################################# Step2 ####################################################

import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp
import numpy as np
import scanpy as sc
x = sc.read_csv("sce_exp.csv"); 
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sce.loom", x.X.transpose(), row_attrs, col_attrs)

################################################# Step3 ####################################################

pyscenic grn --num_workers 16 --sparse --output sce.adj.tsv --method grnboost2 sce.loom hs_hgnc_tfs.txt
pyscenic ctx --num_workers 16 sce.adj.tsv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname sce.loom --output sce.regulon.csv --all_modules --mask_dropouts --mode "dask_multiprocessing" --min_genes 10
pyscenic aucell --num_workers 16 sce.loom sce.regulon.csv --output sce_SCENIC.loom

################################################# Step4 ####################################################

setwd("~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/SECNIC/")
library(Seurat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)

data <- readRDS("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/Primary_clonocluster.rds")

sce_SCENIC <- open_loom("/Volumes/herlynm/linux/ychen/Haiyin/metastasis/SECNIC/Tumor/sce_tumor_SCENIC.loom")
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = data@meta.data$barcode_group)
rss <- na.omit(rss)
rssPlot <- plotRSS(rss, zThreshold = 1.1, cluster_columns = FALSE, order_rows = TRUE, thr = 0.01,
                   col.low = '#330066', col.mid = '#66CC66', col.high = '#FFCC33')
rssPlot

regulon_AUC <- regulonAUC@NAMES
data@meta.data = cbind(data@meta.data, t(SummarizedExperiment::assay(regulonAUC[regulon_AUC,])))

for(i in colnames(data@meta.data)[10:275]){
  tiff(filename = paste0("UMAP_",i,".tif"), width = 572, height = 383, res = 100, compression = "lzw")
  p <- FeaturePlot(data, features = i, cols = c("lightgrey", "#FFD700", "#FF6B6B"))
  print(p)
  dev.off()
}

FeaturePlot(data, features = "RB1(+)", cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, features = "LEF1(+)", cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, features = "JUND(+)", cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, features = "MITF(+)", cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, features = "TFAP2A(+)", cols = c("lightgrey", "#E99F9F", "#FF6B6B"))



