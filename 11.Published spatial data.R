library(Seurat)
library(ggplot2)
library(AUCell)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("~/Documents/Wistar/Haiyin/GSE245582_RAW/WM4237/T0_S1/")

data <- Read10X(data.dir = "matrix")
data <- CreateSeuratObject(data, min.cells = 1, min.features = 1)
img <- Read10X_Image(image.dir = "spatial")
img <- img[Cells(data)]
DefaultAssay(img) <- DefaultAssay(data)
data[["slice"]] <- img

data <- subset(data, cells = rownames(GetTissueCoordinates(data)))

mouse_genes <- data[1:15410,]
human_genes <- data[15411:36952,]
mouse_counts <- Matrix::colSums(mouse_genes@assays$RNA$counts)
human_counts <- Matrix::colSums(human_genes@assays$RNA$counts)
mouse_ratio <- mouse_counts / (mouse_counts + human_counts)
data[["percent_mouse"]] <- mouse_ratio

data <- subset(data, percent_mouse <= 0.1)
data <- data[15411:36952,]

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data$celltype <- "Melanoma"


setwd("~/Documents/Wistar/Haiwei/ADI project/published spatial/")
data <- Read10X_h5("CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5")
data <- CreateSeuratObject(data, min.cells = 1, min.features = 1)
img <- Read10X_Image(image.dir = "spatial")
img <- img[Cells(data)]
DefaultAssay(img) <- DefaultAssay(data)
data[["slice"]] <- img

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)

melanoma_markers <- c("MLANA", "PMEL", "TYR", "MITF", "DCT", "SOX10")
stromal_markers  <- c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1")

melanoma_markers <- intersect(melanoma_markers, rownames(data))
stromal_markers  <- intersect(stromal_markers, rownames(data))

data <- AddModuleScore(data, features = list(melanoma_markers), name = "MelanomaScore")
data <- AddModuleScore(data, features = list(stromal_markers), name = "StromalScore")

data$TumorIndex = data$MelanomaScore1 - data$StromalScore1
data$celltype <- "Other"
data$celltype[data$TumorIndex > 1] <- "Melanoma"
SpatialDimPlot(data, group.by = "celltype", cols = c("Melanoma" = "#F8766D", "Other" = "grey80"))


setwd("~/Downloads/GSM8376641/")
data <- Read10X_h5(filename = "filtered_feature_bc_matrix.h5")
data <- CreateSeuratObject(data, min.cells = 1, min.features = 1)
img <- Read10X_Image(image.dir = "spatial")
img <- img[Cells(data)]
DefaultAssay(img) <- DefaultAssay(data)
data[["slice"]] <- img

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)

melanoma_markers <- c("MLANA", "PMEL", "TYR", "MITF", "DCT")
data <- AddModuleScore(data, features = list(melanoma_markers), name = "MelanomaScore")

data$celltype <- "Other"
data$celltype[data$MelanomaScore1 > 0.5] <- "Melanoma"
SpatialDimPlot(data, group.by = "celltype", cols = c("Melanoma" = "#F8766D", "Other" = "grey80"), image.alpha = 0.3)

mel <- subset(data, celltype == "Melanoma")
SpatialFeaturePlot(mel, features = "OLFML3", image.alpha = 0.3)

lipid <- c("ASAH1","ATP6AP1","CCN1","CD44","CEBPD","COX5B","COX7B","CPT1C","EMP1","FUNDC2",
           "HPGD","LGALS3BP","LPCAT2","LXN","MARCKS","MITF","NDUFA8","NDUFS7","PAFAH1B3",
           "PEBP1","PLA1A","PLAAT3","PLTP","PTGR1","RDH5","SMS","TFAP2A","UGCG")
NC <- c("ALDH1A3","AXL","DOCK4","EBF1","EFNA5","EMILIN1","FGF13","ITGA10","MEF2C","NCALD",
        "NEDD9","NOX4","OLFML3","OLIG1","PXN","SCRG1","SEMA3A","SERPINA3","SFRP1","SPARC",
        "TSPAN7","ZEB1","ZMYND8","ZNF521")
NM <- c("S100A3","S100A4","LINC02241","VCX","SPON2","ELMO1","VCX3B","CNIH3","PI15")

counts <- GetAssayData(object = mel, assay = "RNA", layer = "data")

mel$`Lipid signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = lipid)))
mel$`NC signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NC)))
mel$`Not met signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NM)))

SpatialFeaturePlot(mel, features = "Lipid signature", image.alpha = 0.3)
SpatialFeaturePlot(mel, features = "NC signature", image.alpha = 0.3)
SpatialFeaturePlot(mel, features = "Not met signature", image.alpha = 0.3)

library(dplyr)
library(tibble)
library(dbscan)
library(FNN)
library(igraph)
library(purrr)
library(broom)
library(tidyr)

# =========================
# User settings
# =========================
celltype_col <- "celltype"
lipid_col <- "Lipid signature"
nc_col <- "NC signature"

min_component_size <- 15

# =========================
# 1. Extract metadata + coordinates
# =========================
coords <- as.data.frame(Seurat::GetTissueCoordinates(mel))

if ("cell" %in% colnames(coords)) {
  coords <- coords %>% rename(barcode = cell)
} else {
  coords <- coords %>% rownames_to_column("barcode")
}

x_col <- intersect(c("x", "imagecol", "col", "pxl_col_in_fullres"), colnames(coords))[1]
y_col <- intersect(c("y", "imagerow", "row", "pxl_row_in_fullres"), colnames(coords))[1]

df <- mel@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(
    coords %>% select(barcode, x = all_of(x_col), y = all_of(y_col)),
    by = "barcode"
  ) %>%
  mutate(
    is_mel = .data[[celltype_col]] == "Melanoma",
    Lipid = .data[[lipid_col]],
    NC = .data[[nc_col]]
  )

# =========================
# 2. Build spatial neighbors
# =========================
xy_all <- as.matrix(df[, c("x", "y")])

# Estimate Visium spot-neighbor distance
spot_dist <- median(FNN::get.knn(xy_all, k = 1)$nn.dist[, 1], na.rm = TRUE)
eps <- spot_dist * 1.55

nb <- dbscan::frNN(xy_all, eps = eps)$id
nb <- lapply(seq_along(nb), function(i) setdiff(nb[[i]], i))

df$n_all_nb <- lengths(nb)
df$n_mel_nb <- purrr::map_int(nb, ~ sum(df$is_mel[.x], na.rm = TRUE))

# Boundary melanoma spot:
# 1) has non-melanoma neighbor, or
# 2) has too few melanoma neighbors, which catches tissue-edge / isolated edges
df <- df %>%
  mutate(
    is_boundary = is_mel & (n_mel_nb < n_all_nb | n_mel_nb < 5)
  )

# =========================
# 3. Connected components among melanoma spots
# =========================
mel_barcodes <- df$barcode[df$is_mel]

edge_tbl <- purrr::map_dfr(which(df$is_mel), function(i) {
  j <- nb[[i]]
  j <- j[df$is_mel[j]]
  
  if (length(j) == 0) {
    return(tibble(from = character(), to = character()))
  }
  
  tibble(
    from = df$barcode[i],
    to = df$barcode[j]
  )
}) %>%
  mutate(
    from2 = pmin(from, to),
    to2 = pmax(from, to)
  ) %>%
  distinct(from2, to2) %>%
  transmute(from = from2, to = to2)

g_mel <- igraph::graph_from_data_frame(
  edge_tbl,
  directed = FALSE,
  vertices = tibble(name = mel_barcodes)
)

comp <- igraph::components(g_mel)$membership

df$component <- NA_character_
df$component[match(names(comp), df$barcode)] <- paste0("tumor_", comp)

comp_size <- df %>%
  filter(is_mel) %>%
  count(component, name = "component_n")

df <- df %>%
  left_join(comp_size, by = "component")

# =========================
# 4. Compute tumor depth:
# shortest graph distance to nearest boundary melanoma spot
# =========================

df$edge_dist <- NA_real_
df$depth_norm <- NA_real_
df$is_depth_boundary <- FALSE

valid_components <- comp_size$component[comp_size$component_n >= min_component_size]

for (cc in valid_components) {
  
  # melanoma spots in this connected tumor component
  verts <- df$barcode[df$component == cc & df$is_mel]
  gc <- igraph::induced_subgraph(g_mel, vids = verts)
  
  if (igraph::vcount(gc) < 3) {
    next
  }
  
  # use pre-defined boundary melanoma spots
  boundary_verts <- df$barcode[df$component == cc & df$is_boundary]
  boundary_verts <- intersect(boundary_verts, igraph::V(gc)$name)
  
  # Important fallback:
  # If no boundary is detected, or if every spot is called boundary,
  # use low-degree melanoma spots as the graph boundary.
  if (length(boundary_verts) == 0 || length(boundary_verts) == igraph::vcount(gc)) {
    
    deg <- igraph::degree(gc)
    
    # Take the lowest-degree 20% spots as boundary
    n_boundary <- max(1, ceiling(0.20 * length(deg)))
    boundary_verts <- names(sort(deg, decreasing = FALSE))[seq_len(n_boundary)]
  }
  
  # shortest graph distance from each melanoma spot to nearest boundary spot
  D <- igraph::distances(
    gc,
    v = boundary_verts,
    to = igraph::V(gc),
    weights = NA
  )
  
  dmin <- apply(D, 2, min)
  dmin[!is.finite(dmin)] <- NA_real_
  names(dmin) <- igraph::V(gc)$name
  
  # store edge distance only
  df$edge_dist[match(names(dmin), df$barcode)] <- dmin
  df$is_depth_boundary[match(boundary_verts, df$barcode)] <- TRUE
}

# =========================
# 5. Normalize edge distance within each tumor component
# edge = 0, deepest/core-like spots = 1
# =========================

df <- df %>%
  group_by(component) %>%
  mutate(
    max_edge_dist = if (all(is.na(edge_dist))) {
      NA_real_
    } else {
      max(edge_dist, na.rm = TRUE)
    },
    depth_norm = case_when(
      is_mel & !is.na(edge_dist) & max_edge_dist > 0 ~ edge_dist / max_edge_dist,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()


# =========================
# 4. Compute tumor depth:
# shortest graph distance to nearest boundary melanoma spot
# =========================

df$edge_dist <- NA_real_
df$depth_norm <- NA_real_
df$is_depth_boundary <- FALSE

valid_components <- comp_size$component[comp_size$component_n >= min_component_size]

for (cc in valid_components) {
  
  # melanoma spots in this connected tumor component
  verts <- df$barcode[df$component == cc & df$is_mel]
  gc <- igraph::induced_subgraph(g_mel, vids = verts)
  
  if (igraph::vcount(gc) < 3) {
    next
  }
  
  # use pre-defined boundary melanoma spots
  boundary_verts <- df$barcode[df$component == cc & df$is_boundary]
  boundary_verts <- intersect(boundary_verts, igraph::V(gc)$name)
  
  # Important fallback:
  # If no boundary is detected, or if every spot is called boundary,
  # use low-degree melanoma spots as the graph boundary.
  if (length(boundary_verts) == 0 || length(boundary_verts) == igraph::vcount(gc)) {
    
    deg <- igraph::degree(gc)
    
    # Take the lowest-degree 20% spots as boundary
    n_boundary <- max(1, ceiling(0.20 * length(deg)))
    boundary_verts <- names(sort(deg, decreasing = FALSE))[seq_len(n_boundary)]
  }
  
  # shortest graph distance from each melanoma spot to nearest boundary spot
  D <- igraph::distances(
    gc,
    v = boundary_verts,
    to = igraph::V(gc),
    weights = NA
  )
  
  dmin <- apply(D, 2, min)
  dmin[!is.finite(dmin)] <- NA_real_
  names(dmin) <- igraph::V(gc)$name
  
  # store edge distance only
  df$edge_dist[match(names(dmin), df$barcode)] <- dmin
  df$is_depth_boundary[match(boundary_verts, df$barcode)] <- TRUE
}

# =========================
# 5. Normalize edge distance within each tumor component
# edge = 0, deepest/core-like spots = 1
# =========================

df <- df %>%
  group_by(component) %>%
  mutate(
    max_edge_dist = if (all(is.na(edge_dist))) {
      NA_real_
    } else {
      max(edge_dist, na.rm = TRUE)
    },
    depth_norm = case_when(
      is_mel & !is.na(edge_dist) & max_edge_dist > 0 ~ edge_dist / max_edge_dist,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

ana <- df %>%
  filter(
    is_mel,
    component_n >= min_component_size,
    !is.na(depth_norm),
    !is.na(Lipid),
    !is.na(NC)
  ) %>%
  mutate(
    component = droplevels(factor(component)),
    #Lipid_z = as.numeric(scale(Lipid)),
    #NC_z = as.numeric(scale(NC))
    Lipid_z = as.numeric(Lipid),
    NC_z = as.numeric(NC)
  )

plot_df <- ana %>%
  select(barcode, component, x, y, depth_norm, Lipid_z, NC_z) %>%
  pivot_longer(
    cols = c(Lipid_z, NC_z),
    names_to = "signature",
    values_to = "zscore"
  )

ggplot(ana, aes(x = x, y = y, color = depth_norm)) +
  geom_point(size = 1.6) +
  coord_fixed() +
  scale_y_reverse() +
  theme_classic() +
  labs(color = "Tumor depth")

plot_df <- ana %>%
  select(barcode, component, depth_norm, Lipid_z, NC_z) %>%
  mutate(
    depth_group = case_when(
      depth_norm == 0 ~ "Edge",
      depth_norm > 0.65 ~ "Core",
      TRUE ~ "Middle"
    ),
    depth_group = factor(depth_group, levels = c("Edge", "Middle", "Core"))
  ) %>%
  pivot_longer(
    cols = c(Lipid_z, NC_z),
    names_to = "signature",
    values_to = "zscore"
  ) %>%
  mutate(
    signature = recode(
      signature,
      Lipid_z = "Lipid signature",
      NC_z = "NC signature"
    )
  )

ggplot(plot_df, aes(x = depth_group, y = zscore, fill = depth_group)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.75,
    color = "black"
  ) +
  geom_jitter(
    aes(color = depth_group),
    width = 0.15,
    size = 1.8,
    alpha = 0.7
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 0.8,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  facet_wrap(~ signature, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Edge" = "#4DBBD5",
      "Middle" = "#00A087",
      "Core" = "#E64B35"
    )
  ) +
  scale_color_manual(
    values = c(
      "Edge" = "#4DBBD5",
      "Middle" = "#00A087",
      "Core" = "#E64B35"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(
    x = "Tumor depth group",
    y = "Signature score"
  )

wilcox.test(x = plot_df$zscore[which(plot_df$signature == "Lipid signature" & 
                                       plot_df$depth_group == "Edge")],
            y = plot_df$zscore[which(plot_df$signature == "Lipid signature" & 
                                       plot_df$depth_group == "Middle")])


overlap_test <- function(q = 0.80) {
  
  tmp <- ana %>%
    group_by(component) %>%
    mutate(
      Lipid_high = Lipid >= quantile(Lipid, q, na.rm = TRUE),
      NC_high = NC >= quantile(NC, q, na.rm = TRUE)
    ) %>%
    ungroup()
  
  a <- sum(tmp$Lipid_high & tmp$NC_high)
  b <- sum(tmp$Lipid_high & !tmp$NC_high)
  c <- sum(!tmp$Lipid_high & tmp$NC_high)
  d <- sum(!tmp$Lipid_high & !tmp$NC_high)
  
  tab <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Lipid_high = c("TRUE", "FALSE"),
      NC_high = c("TRUE", "FALSE")
    )
  )
  
  ft <- fisher.test(tab, alternative = "less")
  
  tibble(
    quantile_cutoff = q,
    observed_overlap = a,
    expected_overlap = sum(tmp$Lipid_high) * sum(tmp$NC_high) / nrow(tmp),
    odds_ratio = unname(ft$estimate),
    fisher_p_less_overlap = ft$p.value,
    jaccard = a / sum(tmp$Lipid_high | tmp$NC_high)
  )
}

overlap_results <- map_dfr(c(0.75, 0.80, 0.85, 0.90), overlap_test)
overlap_results

overlap_perm_test <- function(q = 0.80, B = 1000) {
  
  tmp <- ana %>%
    group_by(component) %>%
    mutate(
      Lipid_high = Lipid >= quantile(Lipid, q, na.rm = TRUE),
      NC_high = NC >= quantile(NC, q, na.rm = TRUE)
    ) %>%
    ungroup()
  
  obs <- sum(tmp$Lipid_high & tmp$NC_high)
  
  perm <- replicate(B, {
    tmp2 <- tmp %>%
      group_by(component) %>%
      mutate(NC_high_perm = sample(NC_high)) %>%
      ungroup()
    
    sum(tmp2$Lipid_high & tmp2$NC_high_perm)
  })
  
  tibble(
    quantile_cutoff = q,
    observed_overlap = obs,
    perm_expected_overlap = mean(perm),
    p_less_overlap = (sum(perm <= obs) + 1) / (B + 1)
  )
}

perm_overlap_results <- map_dfr(c(0.75, 0.80, 0.85, 0.90), overlap_perm_test)
perm_overlap_results

