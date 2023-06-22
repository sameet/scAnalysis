# This file will hold generic functions for processing samples and visualizations.
# The visualizations are a bit customized and the custom code will be held here.
# I will add this to a Git repo so we can use the same code on all of my analysis platforms.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(glmGamPoi))

make_opdir <- function(project_name) {
  # Create output directory to hold all the graphs generated for each sample
  of_dir <- paste(project_name, "_analysis", sep = "")
  command <- paste("mkdir ", of_dir, sep = "")
  system(command)
  return(of_dir)
}

process_sample <- function(dir_name, project_name = "scRNA") {
  
  
  of_dir <- make_opdir(dir_name)
  
  
  # create an output path for the sample
  ofn <- paste(paste(project_name, Sys.Date(), sep = "-"), ".RDS", sep = "")
  ofn <- file.path(of_dir, ofn)
  print(ofn)
  
  # Read in the count matrix from cellranger output
  exp_mat <- Read10X(dir_name)
  sce <- CreateSeuratObject(exp_mat, project = project_name)
  sce[["percent.mito"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  
  # plot QC
  p <- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  ggsave(file.path(of_dir, "voilin_plot_QC.pdf"), plot = p, width = 14, height = 7, units = "in")
  
  p2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mito")
  p3 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p4 <- p2 + p3 + plot_layout(guides = "collect")
  ggsave(file.path(of_dir, "qc-plot-2.pdf"), plot = p4, width = 14, height = 7, units = "in")
  
  sce <- SCTransform(sce, method = "glmGamPoi", vars.to.regress = "percent.mito", verbose = F)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- RunPCA(sce) |>
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters(verbose = FALSE)
  saveRDS(sce, ofn)
  return(sce)
}

integrate_sce <- function(sce_list) {
  use_features <- SelectIntegrationFeatures(object.list = sce_list)
  use_sce_list <- lapply(sce_list, function(x) {
    x <- ScaleData(x, features = use_features, verbose = FALSE)
    x <- RunPCA(x, features = use_features, verbose = FALSE)
  })
  use_anchors <- FindIntegrationAnchors(object.list = use_sce_list, anchor.features = use_features)
  combined <- IntegrateData(anchorset = use_anchors)
  
  DefaultAssay(combined) <- "integrated"
  
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
  combined <- FindClusters(combined, resolution = 0.5)
  return(combined)
}