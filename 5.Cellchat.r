setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)
library(CellChat)

primary <- readRDS("Primary_clonocluster.rds")
cellchat_primary <- createCellChat(object = primary, meta = primary@meta.data, group.by = "barcode_group", datatype = "RNA")
cellchat_primary@DB <- CellChatDB.human
cellchat_primary <- subsetData(cellchat_primary)
cellchat_primary <- identifyOverExpressedGenes(cellchat_primary)
cellchat_primary <- identifyOverExpressedInteractions(cellchat_primary)
cellchat_primary <- computeCommunProb(cellchat_primary, type = "triMean")
cellchat_primary <- filterCommunication(cellchat_primary, min.cells = 10)
cellchat_primary <- computeCommunProbPathway(cellchat_primary)
cellchat_primary <- aggregateNet(cellchat_primary)
saveRDS(object = cellchat_primary, file = "cellchat/Primary_cellchat.rds")

lung <- readRDS("Lung_clonocluster.rds")
cellchat_lung <- createCellChat(object = lung, meta = lung@meta.data, group.by = "barcode_group", datatype = "RNA")
cellchat_lung@DB <- CellChatDB.human
cellchat_lung <- subsetData(cellchat_lung)
cellchat_lung <- identifyOverExpressedGenes(cellchat_lung)
cellchat_lung <- identifyOverExpressedInteractions(cellchat_lung)
cellchat_lung <- computeCommunProb(cellchat_lung, type = "triMean")
cellchat_lung <- filterCommunication(cellchat_lung, min.cells = 10)
cellchat_lung <- computeCommunProbPathway(cellchat_lung)
cellchat_lung <- aggregateNet(cellchat_lung)
saveRDS(object = cellchat_lung, file = "cellchat/Lung_cellchat.rds")

liver <- readRDS("Liver_clonocluster.rds")
cellchat_liver <- createCellChat(object = liver, meta = liver@meta.data, group.by = "barcode_group", datatype = "RNA")
cellchat_liver@DB <- CellChatDB.human
cellchat_liver <- subsetData(cellchat_liver)
cellchat_liver <- identifyOverExpressedGenes(cellchat_liver)
cellchat_liver <- identifyOverExpressedInteractions(cellchat_liver)
cellchat_liver <- computeCommunProb(cellchat_liver, type = "triMean")
cellchat_liver <- filterCommunication(cellchat_liver, min.cells = 10)
cellchat_liver <- computeCommunProbPathway(cellchat_liver)
cellchat_liver <- aggregateNet(cellchat_liver)
saveRDS(object = cellchat_liver, file = "cellchat/Liver_cellchat.rds")

setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/cellchat/")
library(CellChat)

cellchat_primary <- readRDS("Primary_cellchat.rds")
cellchat_lung <- readRDS("Lung_cellchat.rds")
cellchat_liver <- readRDS("Liver_cellchat.rds")

net_primary <- as.data.frame(as.table(cellchat_primary@netP$prob))
net_lung <- as.data.frame(as.table(cellchat_lung@netP$prob))
net_liver <- as.data.frame(as.table(cellchat_liver@netP$prob))

colnames(net_primary) <- c("source", "target", "pathway", "prob")
colnames(net_lung) <- c("source", "target", "pathway", "prob")
colnames(net_liver) <- c("source", "target", "pathway", "prob")

net_primary$sample <- "Primary"
net_lung$sample <- "Lung"
net_liver$sample <- "Liver"

net_df <- rbind(net_primary, net_lung, net_liver)

library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

net_df2 <- net_df %>%
  mutate(prob = as.numeric(prob),
         prob = ifelse(is.na(prob), 0, prob),
         row_id = paste0(sample, " | ", source, " \u2192 ", target))

wide_df <- net_df2 %>%
  select(row_id, sample, pathway, prob) %>%
  pivot_wider(names_from = pathway, values_from = prob, values_fill = 0) %>%
  as.data.frame(check.names = FALSE)          # ← 转成 data.frame

rownames(wide_df) <- wide_df$row_id           # ← 现在安全
ann_row <- data.frame(sample = wide_df$sample, row.names = rownames(wide_df))

mat <- as.matrix(wide_df[, setdiff(colnames(wide_df), c("row_id", "sample"))])
mat[which(mat < 0.2)] <- 0
keep_rows <- rowSums(mat) > 0; keep_cols <- colSums(mat) > 0
mat_plot <- mat[keep_rows, keep_cols, drop = FALSE]
ann_plot <- ann_row[keep_rows, , drop = FALSE]

sample_levels <- sort(unique(ann_plot$sample))
sample_colors <- setNames(
  RColorBrewer::brewer.pal(max(3, length(sample_levels)), "Set2")[seq_along(sample_levels)],
  sample_levels
)

pheatmap(mat_plot,
         annotation_row = ann_plot,
         annotation_colors = list(sample = sample_colors),
         cluster_rows = FALSE, cluster_cols = TRUE,
         scale = "none", fontsize_row = 6, fontsize_col = 7, border_color = NA)
