# Visium Isometric Spatial Stacked Visualization - Documentation

This module provides a tool for visualizing Visium spatial transcriptomics and SpatialCellChat communication fields in a 3D-like isometric stacked layout. By stacking multiple layers of biological information (e.g., gene expression, cluster assignments, and cell-cell communication flow fields), researchers can compare spatial structures and signaling patterns.

---

## Prerequisites

Ensure the following R packages are installed and loaded in your environment:
* `Seurat` (v5.0 or later recommended)
* `ggplot2`
* `ggnewscale` (for independent scales and legends per layer)
* `dplyr`
* `scatterpie` (required for scatter pie layers)
* `SpatialCellChat` (required for plotting communication fields)

---

## Function: `plot_isometric_spatial`

### Signature
```R
plot_isometric_spatial(
  seurat_obj = NULL,
  layers = list(),
  image = NULL,
  z_step = 0.8,
  theta = 30,
  phi = 40,
  pad = 0.05,
  point_size = 1.0,
  normalize_independently = FALSE,
  reverse_y = TRUE,
  output_dir = ".",
  filename_base = "isometric_spatial_plot",
  width = 10,
  height = 8
)
```

### Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `seurat_obj` | `Seurat` | `NULL` | Default Seurat object to use if layers don't specify their own. |
| `layers` | `list` | `list()` | List of layer definitions (1 to 3 layers). See [Layer Schema](#layer-schema) below. |
| `image` | `character` | `NULL` | Image FOV to use (e.g., `"WT"`, `"SRT"`). If `NULL`, uses the first image or active FOV. |
| `z_step` | `numeric` | `0.8` | Vertical spacing between stacked layers. Reduce this to bring layers closer together. |
| `theta` | `numeric` | `30` | Rotation angle around the Z-axis in degrees. |
| `phi` | `numeric` | `40` | Tilt/elevation angle in degrees. |
| `pad` | `numeric` | `0.05` | Margin size around the coordinates when drawing the background platforms. |
| `point_size` | `numeric` | `1.0` | Size of Visium spots. Reduce this for high-resolution rendering. |
| `normalize_independently` | `logical` | `FALSE` | If `TRUE`, normalizes each layer's coordinates to `[0, 1]` independently. If `FALSE`, normalizes globally across layers to preserve visual alignment. |
| `reverse_y` | `logical` | `TRUE` | If `TRUE` (default), reverses the Y coordinate axis so the top of the tissue section is plotted at the top of the layer visual. |
| `output_dir` | `character` | `"."` | Directory where output files will be written. |
| `filename_base` | `character` | `"isometric_spatial_plot"` | Base filename for output. |
| `width` | `numeric` | `10` | Width of the output figure in inches. |
| `height` | `numeric` | `8` | Height of the output figure in inches. |

---

## Layer Schema

The `layers` argument is a list containing up to **3** layer lists. The layers are drawn in the order they are added (Layer 1 is at the bottom, Layer 3 is at the top).

### 1. Seurat Feature Layer (Gene Expression or Metadata)
```R
list(
  type = "seurat",                  # optional, default
  seurat_obj = x,                   # optional, overrides global seurat_obj
  feature = "Gfap",                 # character. Gene name or metadata column
  alpha = 0.5,                      # numeric. Spot opacity (0 to 1)
  legend_title = "Astrocyte Marker",# character. Custom legend title
  palette = "viridis",              # character or vector. Predefined palette (e.g. "magma", "Set1") or vector of custom colors
  assay = "SCT",                    # optional. Assay name
  layer = "data"                    # optional. Layer name
)
```

### 2. Scatter Pie Layer (Cell Type Proportions or Multiple Genes)
```R
list(
  type = "scatterpie",              # required
  seurat_obj = x,                   # optional
  features = c("Astro", "Micro"),   # required. Vector of features
  alpha = 0.6,                      # numeric. Pie opacity
  legend_title = "Cell Types",      # character. Custom legend title
  palette = "Set2",                 # character or vector. Palette for pies
  pie_radius = 0.01                 # optional. Manual radius for pies
)
```

### 3. SpatialCellChat Field Layer (Communication Flow)
```R
list(
  type = "cellchat",                # required
  spatial_cellchat_obj = chat_wt,   # required. SpatialCellChat object
  pathway = "SPP1",                 # required. Signaling pathway name
  pattern = "incoming",             # required. "incoming" or "outgoing"
  alpha = 0.4,                      # numeric. Spot opacity (0 to 1)
  legend_title = "SPP1 Targets",    # character. Custom legend title
  arrow_color = "blue",             # character. Arrow color ("map" to match cell colors, or static name like "black")
  arrow_scale_multiplier = 1.2      # numeric. Scaling multiplier for arrow lengths
)
```

> [!IMPORTANT]
> **SpatialCellChat Prerequisite & Robust Handling**: Before plotting a `cellchat` layer, standard `SpatialCellChat` cell-cell communication probability calculation must have been run on the object.
> - **On-the-Fly Field Computation**: If the spatial communication field for the requested `pathway` has not yet been computed on the object, the function will automatically attempt to compute it on the fly using `SpatialCellChat::computeCommunField(chat, signaling.name = pathway)`, as long as the package is available.
> - **Error Skipping**: If the calculation fails (e.g., if the pathway is not present in the probability calculations) or the package is missing, a warning will be displayed and the invalid layer will be skipped, while the remaining valid layers will still compile and plot successfully.

---

## Styling & Saturation Tips

1. **Avoid Visual Overcrowding**: Visium datasets contain thousands of spots. For crisp rendering, use a smaller `point_size` (e.g., `0.5` to `0.8`).
2. **Reducing Color Saturation**: Use softer opacity (`alpha`) values on the layer spots (e.g., `0.35` for clusters, `0.4` for cellchat, `0.5` for genes). This yields a clean, muted look on-screen. 
3. **Full Saturation Legends**: The legends will automatically render at full saturation (`alpha = 1.0`) and larger point sizes even if the layers themselves are semi-transparent, ensuring the legends remain legible.
4. **Bringing Layers Closer**: If there is too much whitespace in the stack, decrease `z_step` (e.g., to `0.7` or `0.6`).

---

## Code Examples

### Example 1: 2-Layer Gene Stack
Compare two marker genes side-by-side:
```R
source("isometric_spatial.R")
library(Seurat)

# Load Seurat object
x <- readRDS("seurat_object.rds")

# Plot Gfap on bottom and Lgals3 on top
layers_2 <- list(
  list(
    feature = "Gfap",
    alpha = 0.45,
    legend_title = "Astrocyte Marker (Gfap)",
    palette = c("#eeeeee", "#3b4cc0", "#b40426") # custom white-blue-red
  ),
  list(
    feature = "Lgals3",
    alpha = 0.5,
    legend_title = "Microglia Marker (Lgals3)",
    palette = "magma"
  )
)

plot_isometric_spatial(
  seurat_obj = x,
  layers = layers_2,
  image = "WT",
  z_step = 0.7,
  point_size = 0.55,
  output_dir = "plots",
  filename_base = "2layer_gene_stack"
)
```

### Example 2: 3-Layer Mixed Stack (Clusters + CellChat Field + Gene Expression)
Compare cell cluster boundaries, ligand-receptor flow fields, and activation state expression:
```R
source("isometric_spatial.R")
library(Seurat)
library(SpatialCellChat)

x <- readRDS("seurat_object.rds")
chat_wt <- readRDS("cellchat_object.rds")

layers_3 <- list(
  list(
    type = "seurat",
    feature = "seurat_clusters",
    alpha = 0.35,
    legend_title = "Seurat Clusters",
    palette = "Set1"
  ),
  list(
    type = "cellchat",
    spatial_cellchat_obj = chat_wt,
    pathway = "COLLAGEN",
    pattern = "incoming",
    alpha = 0.4,
    legend_title = "Collagen Incoming Targets",
    arrow_color = "black",
    arrow_scale_multiplier = 1.2
  ),
  list(
    type = "seurat",
    feature = "Gpnmb",
    alpha = 0.5,
    legend_title = "Microglia Activation (Gpnmb)",
    palette = "viridis"
  )
)

plot_isometric_spatial(
  seurat_obj = x,
  layers = layers_3,
  image = "WT",
  z_step = 0.75,
  point_size = 0.55,
  output_dir = "plots",
  filename_base = "3layer_mixed_stack"
)
```
