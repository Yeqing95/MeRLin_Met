################################################# Step1 ####################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)

data <- readRDS(file ="Primary_clonocluster.rds")

M <- as.matrix(GetAssayData(object = data, assay = "RNA", layer = "count"))
saveRDS(object = M, file = "inferCNV/Lung_Matrix.rds")

G <- data.frame(Group = data$barcode_group, row.names = colnames(data))
saveRDS(object = G, file = "inferCNV/Lung_Group.rds")

################################################# Step2 ####################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/inferCNV/")
library(infercnv)

M <- readRDS("Lung_Matrix.rds")
G <- readRDS("Lung_Group.rds")
P <- readRDS("GRCh38_gene_pos_gene_id.rds")

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = M, 
                                    annotations_file = G, 
                                    gene_order_file = P, 
                                    ref_group_names = NULL)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = ".",
                             cluster_by_groups = FALSE,
                             denoise = TRUE,
                             HMM = TRUE,
                             analysis_mode = "subclusters", 
                             resume_mode = TRUE,
                             num_threads = 4)

################################################# Step3 ####################################################
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/inferCNV/")
library(infercnv)
library(Seurat)

list_to_vec <- function(idx_list){
  group_vec <- rep(NA_character_, max(unlist(idx_list)))
  for (grp in names(idx_list)) {
    group_vec[idx_list[[grp]]] <- grp
  }
  return(group_vec)
}

P <- readRDS("../Primary_clonocluster.rds")
Primary <- readRDS("Primary_run.final.infercnv_obj")
Primary_cnv <- round(Primary@expr.data, 6) - 0.002755
Primary_meta <- data.frame(cell = colnames(Primary_cnv),
                           barcode = P$Barcode,
                           group = list_to_vec(Primary@observation_grouped_cell_indices))
rm(P)

Lu <- readRDS("../Lung_clonocluster.rds")
Lung <- readRDS("Lung_run.final.infercnv_obj")
Lung_cnv <- round(Lung@expr.data, 6) - 0.002743
Lung_meta <- data.frame(cell = colnames(Lung_cnv),
                           barcode = Lu$Barcode,
                           group = list_to_vec(Lung@observation_grouped_cell_indices))
rm(Lu)

Li <- readRDS("../Liver_clonocluster.rds")
Liver <- readRDS("Liver_run.final.infercnv_obj")
Liver_cnv <- round(Liver@expr.data, 6) - 0.002454
Liver_meta <- data.frame(cell = colnames(Liver_cnv),
                        barcode = Li$Barcode,
                        group = list_to_vec(Liver@observation_grouped_cell_indices))
rm(Li)

stopifnot(all(colnames(Primary_cnv) == Primary_meta$cell))
stopifnot(all(colnames(Lung_cnv) == Lung_meta$cell))
stopifnot(all(colnames(Liver_cnv) == Liver_meta$cell))

library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

common_genes <- Reduce(
  intersect,
  list(rownames(Primary_cnv), rownames(Lung_cnv), rownames(Liver_cnv))
)

gene_order <- Primary@gene_order
gene_order$gene <- rownames(gene_order)

colnames(gene_order)[1:3] <- c("chr", "start", "stop")

gene_order <- gene_order %>%
  filter(gene %in% common_genes)

chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY",
                as.character(1:22), "X", "Y")

gene_order$chr <- as.character(gene_order$chr)
gene_order$chr <- factor(
  gene_order$chr,
  levels = unique(c(chr_levels, sort(unique(gene_order$chr))))
)

gene_order <- gene_order %>%
  arrange(chr, start)

genes_per_bin <- 100

gene_order <- gene_order %>%
  group_by(chr) %>%
  mutate(
    bin_number = ceiling(row_number() / genes_per_bin),
    bin_id = paste0(chr, "_bin", bin_number)
  ) %>%
  ungroup()

bin_cnv_matrix <- function(mat, gene_order) {
  mat <- mat[gene_order$gene, , drop = FALSE]
  
  bin_sum <- rowsum(mat, group = gene_order$bin_id, reorder = FALSE)
  bin_n <- table(factor(gene_order$bin_id, levels = rownames(bin_sum)))
  
  bin_mean <- sweep(bin_sum, 1, as.numeric(bin_n), "/")
  return(bin_mean)
}

Primary_bin <- bin_cnv_matrix(Primary_cnv, gene_order)
Lung_bin <- bin_cnv_matrix(Lung_cnv, gene_order)
Liver_bin <- bin_cnv_matrix(Liver_cnv, gene_order)

aggregate_by_barcode <- function(mat, meta, min_cells = 3) {
  meta <- meta[match(colnames(mat), meta$cell), ]
  stopifnot(all(meta$cell == colnames(mat)))
  
  keep <- !is.na(meta$barcode) & meta$barcode != ""
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  
  barcode_counts <- table(meta$barcode)
  keep_barcodes <- names(barcode_counts)[barcode_counts >= min_cells]
  
  keep <- meta$barcode %in% keep_barcodes
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  
  grp <- meta$barcode
  
  summed <- rowsum(t(mat), group = grp, reorder = FALSE)
  counts <- table(factor(grp, levels = rownames(summed)))
  
  averaged <- t(sweep(summed, 1, as.numeric(counts), "/"))
  
  return(averaged)
}

Primary_barcode_cnv <- aggregate_by_barcode(Primary_bin, Primary_meta, min_cells = 3)
Lung_barcode_cnv <- aggregate_by_barcode(Lung_bin, Lung_meta, min_cells = 3)
Liver_barcode_cnv <- aggregate_by_barcode(Liver_bin, Liver_meta, min_cells = 3)

make_barcode_group_lookup <- function(...) {
  meta_all <- bind_rows(...)
  
  meta_all %>%
    filter(!is.na(barcode), !is.na(group)) %>%
    count(barcode, group, sort = TRUE) %>%
    group_by(barcode) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(barcode, group)
}

barcode_group_lookup <- make_barcode_group_lookup(
  Primary_meta %>% mutate(sample = "Primary"),
  Lung_meta %>% mutate(sample = "Lung"),
  Liver_meta %>% mutate(sample = "Liver")
)

primary_barcodes <- colnames(Primary_barcode_cnv)
lung_barcodes <- colnames(Lung_barcode_cnv)
liver_barcodes <- colnames(Liver_barcode_cnv)

matched_barcodes <- primary_barcodes[
  primary_barcodes %in% union(lung_barcodes, liver_barcodes)
]

matched_barcodes <- intersect(
  matched_barcodes,
  barcode_group_lookup$barcode[
    barcode_group_lookup$group %in% c("barcode group 1", "barcode group 2")
  ]
)

mat_list <- list()
anno_list <- list()

for (bc in matched_barcodes) {
  if (bc %in% colnames(Primary_barcode_cnv)) {
    mat_list[[paste(bc, "Primary", sep = "|")]] <- Primary_barcode_cnv[, bc]
  }
  if (bc %in% colnames(Lung_barcode_cnv)) {
    mat_list[[paste(bc, "Lung", sep = "|")]] <- Lung_barcode_cnv[, bc]
  }
  if (bc %in% colnames(Liver_barcode_cnv)) {
    mat_list[[paste(bc, "Liver", sep = "|")]] <- Liver_barcode_cnv[, bc]
  }
}

heat_mat <- do.call(cbind, mat_list)

anno_col <- data.frame(
  barcode = sub("\\|.*", "", colnames(heat_mat)),
  sample = sub(".*\\|", "", colnames(heat_mat))
)

anno_col <- anno_col %>%
  left_join(barcode_group_lookup, by = "barcode") %>%
  as.data.frame()

rownames(anno_col) <- colnames(heat_mat)

plot_mat <- log2(pmax(heat_mat, 1e-6))

pheatmap(
  plot_mat,
  cluster_rows = T,
  cluster_cols = T,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = anno_col,
  color = colorRampPalette(c("blue", "white", "red"))(101),
  breaks = seq(-0.5, 0.5, length.out = 102),
  main = "Barcode-level CNV profiles across primary and metastatic tumors"
)

make_delta_matrix <- function(primary_mat, met_mat, comparison_name) {
  common_bc <- intersect(colnames(primary_mat), colnames(met_mat))
  
  delta_list <- lapply(common_bc, function(bc) {
    log2(pmax(met_mat[, bc], 1e-6)) - log2(pmax(primary_mat[, bc], 1e-6))
  })
  
  delta_mat <- do.call(cbind, delta_list)
  colnames(delta_mat) <- paste(common_bc, comparison_name, sep = "|")
  
  return(delta_mat)
}

delta_lung <- make_delta_matrix(
  Primary_barcode_cnv,
  Lung_barcode_cnv,
  "Lung_vs_Primary"
)

delta_liver <- make_delta_matrix(
  Primary_barcode_cnv,
  Liver_barcode_cnv,
  "Liver_vs_Primary"
)

delta_mat <- cbind(delta_lung, delta_liver)

anno_delta <- data.frame(
  barcode = sub("\\|.*", "", colnames(delta_mat)),
  comparison = sub(".*\\|", "", colnames(delta_mat))
)

anno_delta <- anno_delta %>%
  left_join(barcode_group_lookup, by = "barcode") %>%
  as.data.frame()

rownames(anno_delta) <- colnames(delta_mat)

pheatmap(
  delta_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = anno_delta,
  color = colorRampPalette(c("blue", "white", "red"))(101),
  breaks = seq(-0.3, 0.3, length.out = 102),
  main = "CNV changes in metastases relative to matched primary barcodes"
)

compute_barcode_similarity <- function(primary_mat, met_mat, comparison_name) {
  common_bc <- intersect(colnames(primary_mat), colnames(met_mat))
  
  res <- lapply(common_bc, function(bc) {
    x <- log2(pmax(primary_mat[, bc], 1e-6))
    y <- log2(pmax(met_mat[, bc], 1e-6))
    
    data.frame(
      barcode = bc,
      comparison = comparison_name,
      spearman_cor = suppressWarnings(cor(x, y, method = "spearman")),
      pearson_cor = suppressWarnings(cor(x, y, method = "pearson")),
      mean_abs_delta = mean(abs(y - x), na.rm = TRUE)
    )
  })
  
  bind_rows(res)
}

sim_lung <- compute_barcode_similarity(
  Primary_barcode_cnv,
  Lung_barcode_cnv,
  "Primary vs Lung"
)

sim_liver <- compute_barcode_similarity(
  Primary_barcode_cnv,
  Liver_barcode_cnv,
  "Primary vs Liver"
)

sim_df <- bind_rows(sim_lung, sim_liver) %>%
  left_join(barcode_group_lookup, by = "barcode")

ggplot(sim_df, aes(x = group, y = spearman_cor, fill = comparison)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    aes(color = comparison),
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Spearman correlation of CNV profile",
    title = "Matched barcode CNV similarity between primary and metastases"
  )

ggplot(sim_df, aes(x = group, y = mean_abs_delta, fill = comparison)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    aes(color = comparison),
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean absolute CNV difference",
    title = "Magnitude of CNV changes in matched metastatic barcodes"
  )


