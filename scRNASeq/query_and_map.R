#!/usr/bin/env Rscript

# Script to query CellxGene using metadata from a Parquet file
# and map ENSEMBL IDs to Gene Symbols with count aggregation.

suppressPackageStartupMessages({
  library(arrow)
  library(cellxgene.census)
  library(Seurat)
  library(org.Hs.eg.db)
})

# Source generic functions from the workspace
if (file.exists("sc_functions.R")) {
  source("sc_functions.R")
} else {
  stop("sc_functions.R not found in the current directory.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  message("Usage: Rscript query_and_map.R <metadata.parquet> [output.RDS] [species]")
  message("Example: Rscript query_and_map.R cell_metadata.parquet seurat_mapped.RDS 'Homo sapiens'")
  quit(status = 1)
}

metadata_path <- args[1]
output_path <- if (length(args) >= 2) args[2] else "seurat_mapped.RDS"
species <- if (length(args) >= 3) args[3] else "Homo sapiens"

if (!file.exists(metadata_path)) {
  stop(paste("Metadata file not found:", metadata_path))
}

message(paste("--- Starting Workflow ---"))
message(paste("Metadata:", metadata_path))
message(paste("Output:", output_path))
message(paste("Species:", species))

# 1. Read metadata Parquet file
message("Reading metadata...")
metadata <- arrow::read_parquet(metadata_path)

if (!"soma_joinid" %in% colnames(metadata)) {
  stop("Metadata parquet must contain a 'soma_joinid' column.")
}

use_ids <- metadata$soma_joinid
message(paste("Found", length(use_ids), "cells in metadata."))

# 2. Query CellxGene Census
message("Connecting to CellxGene Census...")
census <- cellxgene.census::open_soma()

message(paste("Fetching Seurat object for", species, "..."))
seurat_obj <- cellxgene.census::get_seurat(
  census = census,
  organism = species,
  obs_coords = use_ids
)

cellxgene.census::close_soma()
message("Successfully fetched data from CellxGene.")

# 3. Map ENSEMBL to Symbol and sum counts
# This function is defined in sc_functions.R
new_seurat <- map_ensembl_to_symbol(seurat_obj)

# 4. Save results
message(paste("Saving processed Seurat object to", output_path, "..."))
saveRDS(new_seurat, output_path)

message("--- Workflow Complete ---")
