# test_isometric.R
# Script to verify plot_isometric_spatial function

source("isometric_spatial.R")
library(Seurat)
library(SpatialCellChat)

message("=== STEP 1: Loading Datasets ===")
seurat_path <- "/Users/sameet/Data/brain_data_analysis/spatial-data-analysis/data-board-spatial/all-4-brain-data/20251030T012536Z-88f14/all-4-brain-data.rds"
cellchat_path <- "/Users/sameet/Data/brain_data_analysis/spatial-data-analysis/spatial-analysis-jun042026/spatial_cellchat_WT_2026-06-15/spatial_cellchat_object.rds"

if (!file.exists(seurat_path)) {
  stop("Seurat RDS file not found!")
}
if (!file.exists(cellchat_path)) {
  stop("SpatialCellChat RDS file not found!")
}

message("Loading Seurat object...")
x <- readRDS(seurat_path)
print(x)

message("Loading SpatialCellChat object...")
chat_wt <- readRDS(cellchat_path)
print(chat_wt)

message("\n=== STEP 2: Running 2-Layer Plot (Seurat Genes) ===")
# Layer 1: Gfap (Astrocyte marker)
# Layer 2: Lgals3 (Microglia marker)
layers_2 <- list(
  list(
    type = "seurat",
    seurat_obj = x,
    feature = "Gfap",
    alpha = 0.45,
    legend_title = "Astrocyte Marker (Gfap)",
    palette = c("#eeeeee", "#3b4cc0", "#b40426") # Custom light-blue-red
  ),
  list(
    type = "seurat",
    seurat_obj = x,
    feature = "Lgals3",
    alpha = 0.5,
    legend_title = "Microglia Marker (Lgals3)",
    palette = "magma" # Predefined Viridis-style palette
  )
)

p2 <- plot_isometric_spatial(
  layers = layers_2,
  image = "WT",
  z_step = 0.7,
  point_size = 0.55,
  output_dir = "test_plots",
  filename_base = "test_2layer_seurat"
)
message("2-Layer Plot completed successfully.")


message("\n=== STEP 3: Running 3-Layer Plot (Mixed Seurat & SpatialCellChat) ===")
# Layer 1: seurat_clusters (Discrete cell groupings)
# Layer 2: SpatialCellChat incoming COLLAGEN communication field (Arrows + discrete sender/receiver labels)
# Layer 3: Gpnmb gene expression (Continuous microglia state marker)
layers_3 <- list(
  list(
    type = "seurat",
    seurat_obj = x,
    feature = "seurat_clusters",
    alpha = 0.35,
    legend_title = "Seurat Clusters",
    palette = "Set1" # RColorBrewer palette name
  ),
  list(
    type = "cellchat",
    spatial_cellchat_obj = chat_wt,
    pathway = "COLLAGEN",
    pattern = "incoming",
    alpha = 0.4,
    legend_title = "Collagen Targets",
    arrow_color = "black", # Draw clear black arrows for flow
    arrow_scale_multiplier = 1.2
  ),
  list(
    type = "seurat",
    seurat_obj = x,
    feature = "Gpnmb",
    alpha = 0.5,
    legend_title = "Gpnmb Expression",
    palette = "viridis"
  )
)

p3 <- plot_isometric_spatial(
  layers = layers_3,
  image = "WT",
  z_step = 0.75,
  point_size = 0.55,
  output_dir = "test_plots",
  filename_base = "test_3layer_mixed"
)
message("3-Layer Plot completed successfully.")

message("\n=== STEP 3b: Running 3-Layer Plot (seurat_clusters, SPP1 incoming, COMPLEMENT incoming) ===")
layers_spp1_comp <- list(
  list(
    type = "seurat",
    seurat_obj = x,
    feature = "seurat_clusters",
    alpha = 0.35,
    legend_title = "Seurat Clusters",
    palette = "Set1"
  ),
  list(
    type = "cellchat",
    spatial_cellchat_obj = chat_wt,
    pathway = "SPP1",
    pattern = "incoming",
    alpha = 0.4,
    legend_title = "SPP1 Targets",
    arrow_color = "blue",
    arrow_scale_multiplier = 1.2
  ),
  list(
    type = "cellchat",
    spatial_cellchat_obj = chat_wt,
    pathway = "COMPLEMENT",
    pattern = "incoming",
    alpha = 0.4,
    legend_title = "Complement Targets",
    arrow_color = "red",
    arrow_scale_multiplier = 1.2
  )
)

p3b <- plot_isometric_spatial(
  layers = layers_spp1_comp,
  image = "WT",
  z_step = 0.75,
  point_size = 0.55,
  output_dir = "test_plots",
  filename_base = "test_3layer_spp1_complement"
)
message("SPP1 & Complement 3-Layer Plot completed successfully.")

message("\n=== STEP 4: Running Input Limit Validation ===")
# Attempt to plot 4 layers to verify the limit validation works
layers_4 <- list(
  list(type = "seurat", seurat_obj = x, feature = "Gfap"),
  list(type = "seurat", seurat_obj = x, feature = "Lgals3"),
  list(type = "seurat", seurat_obj = x, feature = "Gpnmb"),
  list(type = "seurat", seurat_obj = x, feature = "Trem2")
)

tryCatch({
  p4 <- plot_isometric_spatial(layers = layers_4, image = "WT")
  message("WARNING: 4 layers plotting was allowed! Validation failed.")
}, error = function(e) {
  message("Success! 4 layers plot blocked with error: ", e$message)
})

message("\n=== VERIFICATION RUN COMPLETE ===")
