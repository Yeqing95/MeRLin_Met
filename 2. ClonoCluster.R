setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)
library(dplyr)
`%ni%` <- Negate(`%in%`)

data <- readRDS("Liver_filtered.rds")
data <- subset(data, cells = names(which(!is.na(data$Barcode))))

data <- NormalizeData(data)
data <- ScaleData(data)

cm <- t(as.matrix(data@assays$RNA$scale.data)) %>% data.table::as.data.table(keep.rownames = TRUE)
data.table::fwrite(cm, "primary_cm.tsv", sep = "\t")

bt <- data.frame(rn = colnames(data), Barcode = data$Barcode)
data.table::fwrite(bt, "primary_bt.tsv", sep = "\t")

library(magrittr)
library(ClonoCluster)
library(ggplot2)
library(Cairo)

cm <- data.table::fread("liver_cm.tsv")
cm %<>% dt2m
pca <- irlba_wrap(cm, npc = 10)

bt <- data.table::fread("liver_bt.tsv")

clust <- clonocluster(pca, bt, alpha = seq(0, 1, by = 0.1), beta = 0.1, res = 1)

wfs <- seq(0, 10, by = 1)
umaps <- lapply(wfs, function(s){
  uws <- engage_warp(pca, bt, s)
  return(uws)
}) %>% data.table::rbindlist()

umaps <- merge(umaps, bt, by = "rn")
umaps[, Barcode :=
        ifelse(rn %>% unique %>% length > 1, Barcode, "Singlet"),
      by = "Barcode"]

ggplot(umaps, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.3, alpha = 0.5) +
  facet_wrap(~warp) +
  theme_void()

umap_wf6 <- engage_warp(pca, bt, s = 6) %>% select(c(rn, UMAP_1, UMAP_2)) %>% tibble::column_to_rownames(var = "rn")
data@reductions$umap_wf6 <- CreateDimReducObject(embeddings = as.matrix(umap_wf6), key = "UMAP_", assay = DefaultAssay(data))

bt1 <- read.table("met_to_both_bc.txt")$V1
bt2 <- read.table("met_to_liver_bc.txt")$V1
bt3 <- read.table("met_to_lung_bc.txt")$V1
bt4 <- read.table("primary_only_bc.txt")$V1

data$Barcode_type <- NA
data$Barcode_type[which(data$Barcode %in% bt1)] <- "Met to both"
data$Barcode_type[which(data$Barcode %in% bt2)] <- "Met to liver"
data$Barcode_type[which(data$Barcode %in% bt3)] <- "Met to lung"
data$Barcode_type[which(data$Barcode %in% bt4)] <- "Primary only"

DimPlot(data, reduction = "umap_wf6", group.by = "Barcode_type") + 
  scale_color_manual(values = c("Primary only" = "#3C5488FF",
                                "Met to both" = "#E64B35FF", 
                                "Met to lung" = "#00A087FF", 
                                "Met to liver" = "#7E6148FF"))

data <- FindNeighbors(object = data, reduction = "umap_wf6", dims = 1:2)
data <- FindClusters(object = data, resolution = 0.001)
DimPlot(object = data, reduction = "umap_wf6", group.by = "seurat_clusters")

df <- table(data$Barcode, data$seurat_clusters)
max_clusters <- apply(df, 1, function(x) {
  max_indices <- which(x == max(x))
  if (length(max_indices) > 1) {
    return(NA)
  } else {
    return(names(x)[max_indices])
  }
})

highlight_cells <- WhichCells(data, expression = Barcode %in% names(which(max_clusters == "0")))
DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()

group1 <- names(which(max_clusters == "0"))
group2 <- names(which(max_clusters == "1"))

group11 <- read.table("group1_bc.txt")$V1
group22 <- read.table("group2_bc.txt")$V1

write.table(x = union(group1, group11), file = "group1_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = union(group2, group22), file = "group2_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)

data$barcode_group <- NA
data$barcode_group[which(data$Barcode %in% group2)] <- "barcode group 1"
data$barcode_group[which(data$Barcode %in% group1)] <- "barcode group 2"

DimPlot(data, reduction = "umap_wf6", group.by = "barcode_group") +
  scale_fill_manual(values = c("barcode group 1" = "#E64B35FF",
                               "barcode group 2" = "#4DBBD5FF",
                               "barcode group 3" = "#00A087FF"))

View(as.data.frame(table(data$Barcode)))

highlight_cells <- WhichCells(data, expression = Barcode == "TCCTGCAGTA")
DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#D33343")) + NoLegend() + ggtitle("Top 1 BC: TCCTGCAGTA")

highlight_cells <- WhichCells(data, expression = Barcode == "TGGTGTTGAT")
DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#F19939")) + NoLegend() + ggtitle("BC: TGGTGTTGAT")

saveRDS(object = data, file = "Liver_clonocluster.rds")

FeaturePlot(data, feature = "AXL", cols = c("lightgrey", "#FF6B6B"))
FeaturePlot(data, feature = "DCT", cols = c("lightgrey", "#FF6B6B"))
FeaturePlot(data, feature = "PMEL", cols = c("lightgrey", "#FF6B6B"))
FeaturePlot(data, feature = "TFAP2A", cols = c("lightgrey", "#FF6B6B"))









