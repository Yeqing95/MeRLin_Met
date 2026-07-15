setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/")
library(Seurat)
library(monocle3)
library(ggplot2)

data <- readRDS(file = "Integrated.rds")
DimPlot(data, group.by = "Group")

group1 <- read.table("group1_bc.txt")$V1
group2 <- read.table("group2_bc.txt")$V1
data$barcode_group <- "barcode group 3"
data$barcode_group[which(data$Barcode %in% group1)] <- "barcode group 1"
data$barcode_group[which(data$Barcode %in% group2)] <- "barcode group 2"
DimPlot(data, group.by = "barcode_group")

#data <- subset(data, barcode_group == "barcode group 1")
#data <- SCTransform(object = data, vst.flavor = "v2", vars.to.regress = c("percent_mito","S.Score","G2M.Score"))
#data <- RunPCA(data)
#ElbowPlot(data)
#data <- RunUMAP(data, reduction = "pca", dims = 1:10)
#DimPlot(data, group.by = "Group")

expr_matrix <- GetAssayData(data, slot = "count", assay = "SCT")
cell_metadata <- data@meta.data
gene_annotation <- as.data.frame(rownames(data), row.names = rownames(data))
colnames(gene_annotation) <- "gene_short_name"

cds <- new_cell_data_set(expr_matrix, 
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30, method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", 
                        umap.min_dist = 0.3, umap.n_neighbors = 30)

#plot_cells(cds,color_cells_by = "Group")
#plot_cells(cds,color_cells_by = "barcode_group")

reducedDims(cds)$UMAP <- data@reductions$umap@cell.embeddings
cds <- cluster_cells(cds, resolution = 0.005)

cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           cell_size = 0.5, 
           trajectory_graph_color = "grey30")

fit_res <- fit_models(cds, model_formula_str = "~ pseudotime")
coef_res <- coefficient_table(fit_res)
coef_res <- subset(coef_res, term == "pseudotime")
coef_res <- coef_res[,c("gene_id", "num_cells_expressed", "status", "term", "estimate", "std_err", "test_val", "p_value", "normalized_effect", "q_value")]

plot_genes_in_pseudotime(cds["USP13",], color_cells_by = "Group")
plot_genes_in_pseudotime(cds["ROPN1B",], color_cells_by = "pseudotime")

