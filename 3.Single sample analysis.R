setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)
library(dplyr)
library(ggplot2)
`%ni%` <- Negate(`%in%`)

data <- readRDS("Liver_clonocluster.rds")
data <- subset(data, barcode_group %in% c("barcode group 1", "barcode group 2"))

res1 <- FindAllMarkers(data, logfc.threshold = 0.5, min.pct = 0.3, only.pos = T, group.by = "barcode_group")
WriteXLS::WriteXLS(x = res1, ExcelFileName = "~/Documents/Wistar/Haiyin/Metastasis/results/Lung_barcode_group_DEGs.xlsx")

res2 <- FindAllMarkers(data, logfc.threshold = 0.5, min.pct = 0.3, only.pos = T, group.by = "barcode_group")
WriteXLS::WriteXLS(x = res2, ExcelFileName = "~/Documents/Wistar/Haiyin/Metastasis/results/Liver_barcode_group_DEGs.xlsx")

res3 <- FindAllMarkers(data, logfc.threshold = 0.5, min.pct = 0.3, only.pos = T, group.by = "barcode_group")
WriteXLS::WriteXLS(x = res3, ExcelFileName = "~/Documents/Wistar/Haiyin/Metastasis/results/Primary_barcode_group_DEGs_no_g3.xlsx")


res1 <- readxl::read_xlsx("~/Documents/Wistar/Haiyin/Metastasis/results/Lung_barcode_group_DEGs.xlsx")
res2 <- readxl::read_xlsx("~/Documents/Wistar/Haiyin/Metastasis/results/Liver_barcode_group_DEGs.xlsx")
res3 <- readxl::read_xlsx("~/Documents/Wistar/Haiyin/Metastasis/results/Primary_barcode_group_DEGs.xlsx")

res1 <- subset(res1, avg_log2FC >= 1.5 & pct.1 >= 0.4 & p_val_adj < 0.05)
res2 <- subset(res2, avg_log2FC >= 1.5 & pct.1 >= 0.4 & p_val_adj < 0.05)
res3 <- subset(res3, avg_log2FC >= 1.5 & pct.1 >= 0.4 & p_val_adj < 0.05)

library(UpSetR)

upset_data <- fromList(list(
  Primary_group1 = res3$gene[which(res3$cluster == "barcode group 1")], 
  Primary_group2 = res3$gene[which(res3$cluster == "barcode group 2")], 
  #Primary_group3 = res3$gene[which(res3$cluster == "barcode group 3")], 
  Liver_group1 = res2$gene[which(res2$cluster == "barcode group 1")], 
  Liver_group2 = res2$gene[which(res2$cluster == "barcode group 2")], 
  Lung_group1 = res1$gene[which(res1$cluster == "barcode group 1")], 
  Lung_group2 = res1$gene[which(res1$cluster == "barcode group 2")]
))

upset(
  data = upset_data,
  nsets = 6,
  nintersects = NA,
  sets.bar.color = "grey40",
  matrix.color = "black",
  main.bar.color = "grey20",
  point.size = 4,
  line.size = 1.5,
  text.scale = c(2, 1.6, 2, 1.6, 2, 2.4),
  sets.x.label = "Barcode count",
  order.by = "freq",
  decreasing = c(TRUE, TRUE),
  keep.order = TRUE
)

genes <- intersect(res2$gene[which(res2$cluster == "barcode group 1")], res1$gene[which(res1$cluster == "barcode group 1")])
write.table(x = genes, file = "../../../results/group1_genes.txt", quote = F, row.names = F, col.names = F)

genes <- intersect(res2$gene[which(res2$cluster == "barcode group 2")], res1$gene[which(res1$cluster == "barcode group 2")])
write.table(x = genes, file = "../../../results/group2_genes.txt", quote = F, row.names = F, col.names = F)

library(msigdbr)
library(AUCell)
counts <- GetAssayData(object = data, assay = "RNA", layer = "data")
genesets <- msigdbr(species = "Homo sapiens", collection = "H")
gene_sets <- split(x = genesets$gene_symbol, f = genesets$gs_name)
pathway <- gene_sets$HALLMARK_TNFA_SIGNALING_VIA_NFKB
data$pathway <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = pathway)))
#data$pathway <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = intersect(pathway, VariableFeatures(data)))))
FeaturePlot(data, reduction = "umap_wf6", features = "pathway", cols = c("lightgrey", "#EF3B2C"))

Signature <- read_excel(path = "../../../../MeRLin/Data/Source data/MeRLin_signatures.xlsx", sheet = "Identified signatures")
Stress <- as.character(na.omit(Signature$Stress_like_signature))
NC <- as.character(na.omit(Signature$NC_like_signature))
Lipid <- as.character(na.omit(Signature$Lipid_metabolism))
PI3K <- as.character(na.omit(Signature$PI3K_signaling))
ECM <- as.character(na.omit(Signature$ECM_remodeling))
ME <- as.character(na.omit(Signature$Melanocytic_markers))

data$`Stress-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Stress)))
data$`Neural crest-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NC)))
data$`Lipid metabolism` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Lipid)))
data$`PI3K signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = PI3K)))
data$`ECM remodeling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ECM)))
data$`Melanocytic markers` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ME)))

FeaturePlot(data, reduction = "umap_wf6", features = "Stress-like signature", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, reduction = "umap_wf6", features = "Neural crest-like signature", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, reduction = "umap_wf6", features = "Lipid metabolism", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, reduction = "umap_wf6", features = "PI3K signaling", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, reduction = "umap_wf6", features = "ECM remodeling", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, reduction = "umap_wf6", features = "Melanocytic markers", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))

#######################################################################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)
library(msigdbr)
library(AUCell)
library(ggplot2)
`%ni%` <- Negate(`%in%`)

data <- readRDS("Liver_clonocluster.rds")

genesets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
gene_sets <- split(x = genesets$gene_symbol, f = genesets$gs_name)

expr_matrix <- GetAssayData(data, assay = "RNA", layer = "data")
gsva_scores <- getAUC(AUCell_run(exprMat = expr_matrix, geneSets = gene_sets))
data[['GOMF']] <- CreateAssayObject(data = gsva_scores)

cluster_pathways <- FindAllMarkers(object = data, 
                                   assay = "GOMF", 
                                   logfc.threshold = 1, 
                                   min.pct = 0.5,
                                   return.thresh = 0.01,
                                   group.by = "barcode_group",
                                   only.pos = T,
                                   test.use = "wilcox")

WriteXLS::WriteXLS(x = as.data.frame(cluster_pathways), ExcelFileName = "../../../results/GSEA_GOMF.xlsx", row.names = T)

VlnPlot(object = data, features = "GOBP-POSITIVE-REGULATION-OF-SECONDARY-METABOLITE-BIOSYNTHETIC-PROCESS", group.by = "barcode_group", slot = "data")
FeaturePlot(object = data, features = "GOBP-POSITIVE-REGULATION-OF-SECONDARY-METABOLITE-BIOSYNTHETIC-PROCESS", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))

VlnPlot(object = data, features = "GOBP-DITERPENOID-BIOSYNTHETIC-PROCESS", group.by = "barcode_group", slot = "data")
FeaturePlot(object = data, features = "GOBP-DITERPENOID-BIOSYNTHETIC-PROCESS", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))


#######################################################################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)
library(ggplot2)
library(readxl)
library(Cairo)

data <- readRDS("Primary_clonocluster.rds")
list <- read_excel(path = "~/Desktop/temp.xlsx")

gene_list <- c("DCT", "MITF", "PMEL", "POMC", "EDNRB", "JUN", "TYRP1")
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/")

gene_list <- na.omit(list$Group1_DEGs)
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/Group1_DEGs/")

gene_list <- na.omit(list$Group1_DEGs)
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/Group1_DEGs/")

gene_list <- na.omit(list$`Old Lipid metabolism`)
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/Old Lipid Signature/")

gene_list <- na.omit(list$`New Lipid signature`)
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/New Lipid Signature/")

gene_list <- na.omit(list$`Old NCSC signature`)
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/Old NCSC Signature/")

gene_list <- na.omit(list$`New NCSC signature`)
plot_every_gene(object = data, gene_list = gene_list, 
                dir = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/genes/Primary/New NCSC Signature/")


plot_every_gene <- function(object, gene_list, dir){
  for(i in seq_along(gene_list)){
    gene <- gene_list[i]
    p <- FeaturePlot(object, reduction = "umap_wf6", features = gene, pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
    
    figure_name <- paste0(dir, gene, ".pdf")
    CairoPDF(file = figure_name, width = 30, height = 21, dpi = 300)
    print(p)
    dev.off()
  }
}


#######################################################################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(readxl)
library(Seurat)
library(AUCell)

data <- readRDS("Liver_clonocluster.rds")
expr_matrix <- GetAssayData(data, assay = "RNA", layer = "data")

Signature <- read_excel(path = "../../../2025-09-10 MET_DEG/2025-09-04 Group Signature.xlsx", sheet = "Group 2_Lipid")
Lipid <- as.character(na.omit(Signature$Gene))
Signature <- read_excel(path = "../../../2025-09-10 MET_DEG/2025-09-04 Group Signature.xlsx", sheet = "Group 1_NC")
NC <- as.character(na.omit(Signature$Gene))

NC <- c("ALDH1A3","AXL","EBF1","EFNA5","FGF13","NCALD","NEDD9","NOX4","OLFML3","OLIG1",
        "PXN","SEMA3A","TSPAN7","ZMYND8","ZNF521","DOCK4","EMILIN1","SFRP1","ZEB1",
        "MEF2C","SPARC")

Lipid <- c("HACD1","ST8SIA1","ABCA5","PLPP1","ASAH1","LPCAT2","PLA1A","PLAAT3","PTGR1",
           "TECR","UGCG","ACAA1","PAFAH1B3","PEBP1","PLTP","SCARB1","SMS","AGPAT2",
           "CEBPD","PLD1","ACSL1","HACL1")

data$`Lipid signature` <- as.numeric(getAUC(AUCell_run(exprMat = expr_matrix, geneSets = Lipid)))
data$`Neural crest-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = expr_matrix, geneSets = NC)))

FeaturePlot(data, reduction = "umap_wf6", features = "Lipid signature", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))
FeaturePlot(data, reduction = "umap_wf6", features = "Neural crest-like signature", pt.size = 0.4, cols = c("lightgrey", "#E99F9F", "#FF6B6B"))

#######################################################################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(readxl)
library(Seurat)

data <- readRDS("Liver_clonocluster.rds")
common_bc <- read.table("../../../results/common_bc.txt")

data$Common_barcodes <- "Other"
data$Common_barcodes[which(data$Barcode == common_bc$V1[1])] <- common_bc$V1[1]
data$Common_barcodes[which(data$Barcode == common_bc$V1[2])] <- common_bc$V1[2]
data$Common_barcodes[which(data$Barcode == common_bc$V1[3])] <- common_bc$V1[3]
data$Common_barcodes[which(data$Barcode == common_bc$V1[4])] <- common_bc$V1[4]
data$Common_barcodes[which(data$Barcode == common_bc$V1[5])] <- common_bc$V1[5]
data$Common_barcodes[which(data$Barcode == common_bc$V1[6])] <- common_bc$V1[6]
data$Common_barcodes[which(data$Barcode == common_bc$V1[7])] <- common_bc$V1[7]
data$Common_barcodes[which(data$Barcode == common_bc$V1[8])] <- common_bc$V1[8]
data$Common_barcodes[which(data$Barcode == common_bc$V1[9])] <- common_bc$V1[9]
data$Common_barcodes[which(data$Barcode == common_bc$V1[10])] <- common_bc$V1[10]
data$Common_barcodes[which(data$Barcode == common_bc$V1[11])] <- common_bc$V1[11]
data$Common_barcodes[which(data$Barcode == common_bc$V1[12])] <- common_bc$V1[12]
data$Common_barcodes[which(data$Barcode == common_bc$V1[13])] <- common_bc$V1[13]
data$Common_barcodes[which(data$Barcode == common_bc$V1[14])] <- common_bc$V1[14]
data$Common_barcodes <- factor(x = data$Common_barcodes, 
                               levels = c("TCCTGCAGTA", "ACATGGTCAA", "TCTTGAAGAA", "AGATCTTGTA", "AGGTGCTGGA",
                                          "TCATGTACGT", "TGAACGTGCA", "TGGACAAGCT", "TCCAGATCAT", "AGATCAAGGA", 
                                          "TGGTCTTCAA", "AGTAGGTGCA", "AGCTGTACGT", "TGGTGAAGCA", "Other"))

CairoPDF(file = "~/Documents/Wistar/Haiyin/Metastasis/figure/scRNAseq/Liver_common_bc.pdf", width = 32, height = 21, dpi = 300)
DimPlot(data, reduction = "umap_wf6", group.by = "Common_barcodes") +
  scale_color_manual(values = c(
    "TCCTGCAGTA" = "#F8766D",
    "ACATGGTCAA" = "#E38900",
    "TCTTGAAGAA" = "#C49A00",
    "AGATCTTGTA" = "#99A800",
    "AGGTGCTGGA" = "#53B400",
    "TCATGTACGT" = "#00BC56",
    "TGAACGTGCA" = "#00C094",
    "TGGACAAGCT" = "#00BFC4",
    "TCCAGATCAT" = "#00B6EB",
    "AGATCAAGGA" = "#06A4FF",
    "TGGTCTTCAA" = "#A58AFF",
    "AGTAGGTGCA" = "#DF70F8",
    "AGCTGTACGT" = "#FB61D7",
    "TGGTGAAGCA" = "#FF66A8",
    "Other" = "lightgrey"
  ))
dev.off()

