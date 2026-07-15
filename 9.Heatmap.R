setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch2/")
source("~/Documents/Wistar/Haiyin/barcode extractor/utils.R")

# ---- Step 0: 读取 barcode 向量列表 ----
m287_Tumor <- extract_barcode(sam = "sam/ND238_287_Tumor_Unmapped.sam", len = 10)
m287_Lung <- extract_barcode(sam = "sam/ND238_287_Lung_Unmapped.sam", len = 10)
m292_Tumor <- extract_barcode(sam = "sam/ND238_292_Tumor_Unmapped.sam", len = 10)
m292_Lung <- extract_barcode(sam = "sam/ND238_292_Lung_Unmapped.sam", len = 10)
m292_Liver <- extract_barcode(sam = "sam/ND238_292_Liver_Unmapped.sam", len = 10)
m292_Spleen <- extract_barcode(sam = "sam/ND238_292_Spleen_Unmapped.sam", len = 10)
m294_Tumor <- extract_barcode(sam = "sam/ND238_294_Tumor_Unmapped.sam", len = 10)

library(dplyr)
library(tidyr)
library(openxlsx)
library(pheatmap)

# ---- Step 1: barcode 向量列表 ----
barcode_list <- list(
  m287_Tumor = m287_Tumor,
  m287_Lung = m287_Lung,
  m292_Tumor = m292_Tumor,
  m292_Lung = m292_Lung,
  m292_Liver = m292_Liver,
  m292_Spleen = m292_Spleen,
  m294_Tumor = m294_Tumor
)

# ---- Step 2: 转换为比例矩阵 ----
barcode_counts <- lapply(barcode_list, function(vec) {
  tab <- table(vec)
  prop <- tab / sum(tab)   # 计算比例
  return(as.data.frame(prop))
})

# 合并所有 sample
merged_counts <- Reduce(function(x, y) full_join(x, y, by = c("vec" = "vec")),
                        barcode_counts)
colnames(merged_counts) <- c("Barcode", names(barcode_list))
merged_counts[is.na(merged_counts)] <- 0

# 设置 rownames
rownames(merged_counts) <- merged_counts$Barcode
merged_counts <- merged_counts[,-1]

# ---- Step 3: annotation_col (鼠 & 组织) ----
annotation_col <- data.frame(
  Mouse = gsub("(_.*)", "", colnames(merged_counts)),
  Tissue = gsub(".*_", "", colnames(merged_counts))
)
rownames(annotation_col) <- colnames(merged_counts)

# ---- Step 4: annotation_row (barcode 分类) ----
met_to_both <- read.table("../batch1/processed/met_to_both_bc.txt", stringsAsFactors = FALSE)$V1
met_to_liver <- read.table("../batch1/processed/met_to_liver_bc.txt", stringsAsFactors = FALSE)$V1
met_to_lung <- read.table("../batch1/processed/met_to_lung_bc.txt", stringsAsFactors = FALSE)$V1
primary_only <- read.table("../batch1/processed/primary_only_bc.txt", stringsAsFactors = FALSE)$V1

barcode_type <- data.frame(Barcode = rownames(merged_counts))
barcode_type$Type <- "Other"
barcode_type$Type[barcode_type$Barcode %in% met_to_both] <- "Met_to_both"
barcode_type$Type[barcode_type$Barcode %in% met_to_liver] <- "Met_to_liver"
barcode_type$Type[barcode_type$Barcode %in% met_to_lung] <- "Met_to_lung"
barcode_type$Type[barcode_type$Barcode %in% primary_only] <- "Primary_only"
rownames(barcode_type) <- barcode_type$Barcode
barcode_type <- barcode_type[,-1, drop=FALSE]

# ---- Step 5: pheatmap ----
pheatmap(merged_counts,
         annotation_col = annotation_col,
         annotation_row = barcode_type,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         scale = "none")

# 设置比例阈值
threshold <- 0.025

# 保留全部的类型
always_keep_types <- c("Met_to_both", "Met_to_liver", "Met_to_lung")
filter_types <- c("Primary_only", "Other")

# 找到满足条件的 barcode
keep_barcodes <- rownames(barcode_type)[barcode_type$Type %in% always_keep_types]
filter_barcodes <- rownames(barcode_type)[barcode_type$Type %in% filter_types &
                                            apply(merged_counts, 1, function(x) any(x > threshold))]

# 合并需要保留的 barcode
final_barcodes <- c(keep_barcodes, filter_barcodes)

# ❌ 去掉指定强势 barcode
final_barcodes <- setdiff(final_barcodes, "TCCTGCAGTA")

# 筛选 merged_counts 和 barcode_type
merged_counts_filtered <- merged_counts[final_barcodes, ]
barcode_type_filtered <- barcode_type[final_barcodes, , drop = FALSE]

# 画热图
pheatmap(merged_counts_filtered[,c(1,3,7,2,4,5,6)],
         annotation_col = annotation_col,
         annotation_row = barcode_type_filtered,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "column")









library(dplyr)
library(ggplot2)
library(Cairo)
library(ggsci)

plot_top10_barcode_pie <- function(barcode_vector, sample_name){
  # 统计 barcode 次数并排序
  bc_table <- as.data.frame(table(barcode_vector)) %>%
    arrange(desc(Freq))
  
  # 取前10个，其余合并成 Other
  top10 <- bc_table[1:10, ]
  other_sum <- sum(bc_table$Freq[11:nrow(bc_table)])
  
  plot_data <- rbind(
    top10,
    data.frame(barcode_vector = "Other", Freq = other_sum)
  )
  
  plot_data$barcode_vector <- factor(
    plot_data$barcode_vector,
    levels = unique(plot_data$barcode_vector)
  )
  
  # 计算比例
  plot_data <- plot_data %>%
    slice(n():1) %>%
    mutate(prop = Freq / sum(Freq) * 100,
           ypos = cumsum(prop) - 0.5 * prop)
  
  # 设置颜色（前10用 palette，Other 用灰色）
  my_colors <- c(pal_frontiers("default", alpha = 0.8)(10), "#D9D9D9CC")

  # 输出 PDF
  ggplot(plot_data, aes(x = "", y = prop, fill = barcode_vector)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    scale_fill_manual(values = my_colors) +
    coord_polar("y", start = 0) +
    theme_void() + 
    geom_text(aes(y = ypos, label = paste0(round(prop,1), "%")),
              color = "white", size = 5,
              position = position_nudge(x = 0.2)) +
    theme(text = element_text(size = 14),
          legend.box.margin = margin(-5,-5,-5,-5)) +
    labs(title = paste("Top 10 barcodes -", sample_name), fill = "Barcode")
}

# 使用示例
plot_top10_barcode_pie(
  barcode_vector = m287_Tumor,
  sample_name = "m287_Tumor")

plot_top10_barcode_pie(
  barcode_vector = m287_Lung,
  sample_name = "m287_Lung")

plot_top10_barcode_pie(
  barcode_vector = m292_Tumor,
  sample_name = "m292_Tumor")

plot_top10_barcode_pie(
  barcode_vector = m292_Liver,
  sample_name = "m292_Liver")

plot_top10_barcode_pie(
  barcode_vector = m292_Lung,
  sample_name = "m292_Lung")

plot_top10_barcode_pie(
  barcode_vector = m292_Spleen,
  sample_name = "m292_Spleen")

plot_top10_barcode_pie(
  barcode_vector = m294_Tumor,
  sample_name = "m294_Tumor")

plot_top10_barcode_pie(
  barcode_vector = CC,
  sample_name = "Carbon Copy")



cc.tumor <- read.table("../batch1/original/Primary/Primary_bc.txt", stringsAsFactors = FALSE)$V2
cc.lung <- read.table("../batch1/original/Lung/Lung_bc.txt", stringsAsFactors = FALSE)$V2
cc.liver <- read.table("../batch1/original/Liver/Liver_bc.txt", stringsAsFactors = FALSE)$V2

plot_top10_barcode_pie(
  barcode_vector = cc.tumor,
  sample_name = "scRNAseq Tumor")

plot_top10_barcode_pie(
  barcode_vector = cc.lung,
  sample_name = "scRNAseq Lung")

plot_top10_barcode_pie(
  barcode_vector = cc.liver,
  sample_name = "scRNAseq Liver")



cc.tumor <- read.table("../../../MeRLin/Data/Original data/In-vitro/", stringsAsFactors = FALSE)$V2
cc.lung <- read.table("../batch1/original/Lung/Lung_bc.txt", stringsAsFactors = FALSE)$V2
cc.liver <- read.table("../batch1/original/Liver/Liver_bc.txt", stringsAsFactors = FALSE)$V2

plot_top10_barcode_pie(
  barcode_vector = cc.tumor,
  sample_name = "scRNAseq Tumor")

plot_top10_barcode_pie(
  barcode_vector = cc.lung,
  sample_name = "scRNAseq Lung")






plot_top10_barcode_pie(
  barcode_vector = CC,
  sample_name = "Res CC")

plot_top10_barcode_pie(
  barcode_vector = D0_126,
  sample_name = "Res D0_126")

plot_top10_barcode_pie(
  barcode_vector = D0_130,
  sample_name = "Res D0_130")

plot_top10_barcode_pie(
  barcode_vector = D0_134,
  sample_name = "Res D0_134")

plot_top10_barcode_pie(
  barcode_vector = CTRL,
  sample_name = "Res CTRL")


