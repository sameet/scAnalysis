# This file will hold generic functions for processing samples and visualizations.
# The visualizations are a bit customized and the custom code will be held here.
# I will add this to a Git repo so we can use the same code on all of my analysis platforms.

# suppressPackageStartupMessages(library(Seurat))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(patchwork))
# suppressPackageStartupMessages(library(pheatmap))
# suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
  library(patchwork)
  library(glmGamPoi)
  library(harmony)
})

process_h5 <- function(fn, p_name) {
  of_dir <- make_opdir(p_name)
  ofn <- file.path(
    of_dir,
    paste(paste(p_name, Sys.Date(), sep = "-"), ".RDS", sep = "")
  )

  sce <- CreateSeuratObject(Read10X_h5(fn), project = p_name)

  sce[["percent.mito"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

  # plot QC
  # p <- VlnPlot(
  #   sce,
  #   features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  #   ncol = 3
  # )
  # ggsave(
  #   file.path(of_dir, "voilin_plot_QC.pdf"),
  #   plot = p,
  #   width = 14,
  #   height = 7,
  #   units = "in"
  # )

  # p2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mito")
  # p3 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # p4 <- p2 + p3 + plot_layout(guides = "collect")
  # ggsave(
  #   file.path(of_dir, "qc-plot-2.pdf"),
  #   plot = p4,
  #   width = 14,
  #   height = 7,
  #   units = "in"
  # )

  sce <- SCTransform(
    sce,
    # method = "glmGamPoi",
    vars.to.regress = "percent.mito",
    verbose = FALSE
  )
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- RunPCA(sce, verbose = FALSE) |>
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)
  saveRDS(sce, ofn)
  sce
}

make_opdir <- function(project_name) {
  # Create output directory to hold all the graphs generated for each sample
  of_dir <- paste(project_name, "_analysis", sep = "")
  message(paste("Saving to ", of_dir))
  command <- paste("mkdir ", of_dir, sep = "")
  system(command)
  of_dir
}

process_sample <- function(dir_name, project_name = "scRNA") {
  of_dir <- make_opdir(project_name)

  # create an output path for the sample
  ofn <- file.path(
    of_dir,
    paste(paste(project_name, Sys.Date(), sep = ""), ".RDS", sep = "")
  )
  # ofn <- paste(paste(project_name, Sys.Date(), sep = "-"), ".RDS", sep = "")
  # ofn <- file.path(of_dir, ofn)
  print(ofn)

  # Read in the count matrix from cellranger output
  exp_mat <- Read10X(dir_name)
  sce <- CreateSeuratObject(exp_mat, project = project_name)
  sce[["percent.mito"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

  # plot QC
  p <- VlnPlot(
    sce,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
    ncol = 3
  )
  ggsave(
    file.path(of_dir, "voilin_plot_QC.pdf"),
    plot = p,
    width = 14,
    height = 7,
    units = "in"
  )

  p2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mito")
  p3 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p4 <- p2 + p3 + plot_layout(guides = "collect")
  ggsave(
    file.path(of_dir, "qc-plot-2.pdf"),
    plot = p4,
    width = 14,
    height = 7,
    units = "in"
  )

  sce <- SCTransform(
    sce,
    method = "glmGamPoi",
    vars.to.regress = "percent.mito",
    verbose = F
  )
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- RunPCA(sce, verbose = F) |>
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>%
    FindClusters(verbose = FALSE)
  saveRDS(sce, ofn)
  return(sce)
}

integrate_sce <- function(sce_list, ...) {
  use_features <- SelectIntegrationFeatures(object.list = sce_list)
  # use_sce_list <- lapply(sce_list, function(x) {
  #   x <- ScaleData(x, features = use_features, verbose = FALSE)
  #   x <- RunPCA(x, features = use_features, verbose = FALSE)
  # })
  sce_list <- PrepSCTIntegration(sce_list, anchor.features = use_features)
  use_anchors <- FindIntegrationAnchors(
    object.list = use_sce_list,
    anchor.features = use_features
  )
  combined <- IntegrateData(anchorset = use_anchors, ...)

  DefaultAssay(combined) <- "integrated"

  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
  combined <- FindClusters(combined, resolution = 0.5)
  return(combined)
}

integrate_sce_2 <- function(x, y) {
  use_sce_l <- list(x, y)
  use_features <- SelectIntegrationFeatures(
    object.list = use_sce_l,
    verbose = F
  )
  use_sce_list <- lapply(use_sce_l, function(sc) {
    sc <- ScaleData(sc, features = use_features, verbose = FALSE)
    sc <- RunPCA(sc, features = use_features, verbose = FALSE)
  })
  use_anchors <- FindIntegrationAnchors(
    object.list = use_sce_list,
    anchor.features = use_features,
    reduction = "rpca"
  )
  combined <- IntegrateData(anchorset = use_anchors)

  DefaultAssay(combined) <- "integrated"

  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
  combined <- FindClusters(combined, resolution = 0.5)
  return(combined)
}

get_cluster_colors <- function(sce, clust_name = "seurat_clusters") {
  sce@meta.data |>
    pull(get(clust_name)) |>
    unique() -> clusters

  max_cols <- brewer.pal.info["Paired", ]$maxcolors

  req_colors <- length(clusters)

  if (req_colors > max_colors) {
    use_cols <- colorRampPalette(rev(brewer.pal(max_cols, "Paired")))(
      req_colors
    )
  }

  if (req_colors < max_colors) {
    use_cols <- rev(brewer.pal(req_colors, "Paired"))
  }

  names(use_cols) <- clusters
  use_cols
}

#' Map ENSEMBL IDs to Gene Symbols in a Seurat Object
#'
#' @param seurat_obj A Seurat object with ENSEMBL IDs as rownames.
#' @param species_db An AnnotationDb object (e.g., org.Hs.eg.db).
#' @return A new Seurat object with Gene Symbols as rownames, summed counts for duplicates.
map_ensembl_to_symbol <- function(seurat_obj, species_db = org.Hs.eg.db) {
  suppressPackageStartupMessages({
    library(Matrix)
    library(AnnotationDbi)
  })

  message("Mapping ENSEMBL IDs to Gene Symbols...")
  
  # Get counts (works for V5 and earlier)
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  ensembl_ids <- rownames(counts)

  # Map IDs
  mapping <- AnnotationDbi::select(species_db, keys = ensembl_ids, columns = "SYMBOL", keytype = "ENSEMBL")
  
  # Remove NAs
  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  
  # Filter to ensembl IDs that are in the counts
  mapping <- mapping[mapping$ENSEMBL %in% ensembl_ids, ]
  
  if (nrow(mapping) == 0) {
    stop("No ENSEMBL IDs could be mapped to Symbols.")
  }

  # If one ENSEMBL maps to many Symbols (rare), pick the first one to avoid duplicating counts
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  
  message(paste0("Mapped ", nrow(mapping), " ENSEMBL IDs to ", length(unique(mapping$SYMBOL)), " unique Symbols."))

  # Subset counts to those that were mapped
  counts_subset <- counts[mapping$ENSEMBL, ]
  
  # Aggregate by Symbol using sparse matrix multiplication for efficiency
  unique_symbols <- unique(mapping$SYMBOL)
  mapping$SYMBOL <- factor(mapping$SYMBOL, levels = unique_symbols)
  
  # Create a grouping matrix G: Rows = Symbols, Cols = Ensembl IDs
  i_idx <- as.numeric(mapping$SYMBOL)
  j_idx <- 1:nrow(mapping)
  
  G <- sparseMatrix(i = i_idx, j = j_idx, x = 1, 
                    dims = c(length(unique_symbols), nrow(mapping)))
  
  summed_counts <- G %*% counts_subset
  rownames(summed_counts) <- unique_symbols
  colnames(summed_counts) <- colnames(counts_subset)
  
  # Create new Seurat object
  # Transfer metadata
  new_seurat <- CreateSeuratObject(counts = summed_counts, meta.data = seurat_obj[[]], project = Project(seurat_obj))
  
  return(new_seurat)
}
# make_phate <-

merge_harmony <- function(sce_list) {
  sce_merge <- merge(
    sce_list[[1]],
    y = sce_list[2:length(sce_list)],
    add.cell.ids = TRUE
  )
  sce_merge <- FindVariableFeatures(
    sce_merge,
    selection.method = "vst",
    nfeatures = 2000
  )
  sce_merge <- ScaleData(sce_merge)
  sce_merge <- RunPCA(sce_merge, npcs = 50)
  # Run Harmony to merge the samples

  sce_merge <- RunHarmony(
    sce_merge,
    group.by.vars = "orig.ident",
    # reduction = "pca",
    assay.use = "SCT",
    project.dims = TRUE,
    ncores = 8
  )
  sce_merge
}
