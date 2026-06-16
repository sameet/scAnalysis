# isometric_spatial.R
# Module for Visium spatial transcriptomics isometric stacked visualization

library(Seurat)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(grid)

#' Plot Multiple Spatial Transcriptomics Layers in Isometric 3D Stack
#'
#' @param seurat_obj Default Seurat object to use if layers don't specify their own.
#' @param layers A list of layer definitions (maximum 3). Each layer is a list with:
#'   - \code{type}: "seurat" (default) or "cellchat".
#'   - \code{seurat_obj}: Seurat object (optional, defaults to main \code{seurat_obj}).
#'   - \code{spatial_cellchat_obj}: SpatialCellChat object (required if type is "cellchat").
#'   - \code{feature}: Feature to plot (gene name or metadata column for "seurat" type).
#'   - \code{pathway}: Pathway name (required for "cellchat" type, e.g., "TGFb").
#'   - \code{pattern}: "incoming" or "outgoing" (required for "cellchat" type).
#'   - \code{assay}: Assay name (optional, defaults to DefaultAssay).
#'   - \code{layer}: Layer name (optional, defaults to "data").
#'   - \code{alpha}: Alpha transparency for the spots (default 0.8).
#'   - \code{legend_title}: Custom title for the legend (defaults to Layer 1, Layer 2, etc.).
#'   - \code{palette}: Custom color palette vector or predefined palette name.
#'   - \code{arrow_color}: Color for cellchat arrows. If NULL or "map", colors arrows by cell type (default). Otherwise, a static color string (e.g. "red", "black").
#'   - \code{arrow_scale_multiplier}: Multiplier to adjust the auto-scaled length of communication field arrows (default 1.0).
#' @param image Optional image FOV name to extract coordinates from (e.g., "WT"). If NULL, uses the first image or active FOV.
#' @param z_step Vertical spacing between layers (default 1.2).
#' @param theta Angle of rotation around the Z-axis in degrees (default 30).
#' @param phi Angle of tilt/elevation in degrees (default 40).
#' @param pad Padding to add around the platforms (default 0.05).
#' @param point_size Size of the spots (default 1.5).
#' @param normalize_independently Logical. If TRUE, normalizes each layer's coordinates to [0,1] independently. If FALSE (default), normalizes globally across all layers.
#' @param reverse_y Logical. If TRUE (default), reverses the Y coordinate axis so the top of the tissue section is plotted at the top of the layer visual.
#' @param output_dir Directory to save the output files (default ".").
#' @param filename_base Base filename for the saved PNG and PDF files (default "isometric_spatial_plot").
#' @param width Width of the output figure in inches (default 10).
#' @param height Height of the output figure in inches (default 8).
#' @return A ggplot object representing the isometric stack.
#' @export
plot_isometric_spatial <- function(
  seurat_obj = NULL,
  layers = list(),
  image = NULL,
  z_step = 1.2,
  theta = 30,
  phi = 40,
  pad = 0.05,
  point_size = 1.5,
  normalize_independently = FALSE,
  reverse_y = TRUE,
  output_dir = ".",
  filename_base = "isometric_spatial_plot",
  width = 10,
  height = 8
) {
  # 1. Validation
  if (length(layers) == 0) {
    stop("Please provide at least one layer definition in 'layers'.")
  }
  if (length(layers) > 3) {
    stop("The number of layers is limited to 3 for now.")
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Convert angles to radians
  theta_rad <- theta * pi / 180
  phi_rad <- phi * pi / 180
  
  # 2. Extract Data for Each Layer
  layer_dfs <- list()
  valid_layers <- list()
  
  for (i in seq_along(layers)) {
    lyr <- layers[[i]]
    
    df <- tryCatch({
      lyr_type <- ifelse(is.null(lyr$type), "seurat", lyr$type)
      
      if (lyr_type == "seurat") {
        # Resolve Seurat object
        obj <- lyr$seurat_obj
        if (is.null(obj)) obj <- seurat_obj
        if (is.null(obj)) {
          stop(paste("Layer", i, "does not have a Seurat object specified and no default was provided."))
        }
        
        # Extract coordinates
        if (!is.null(image)) {
          coords <- Seurat::GetTissueCoordinates(obj, image = image)
        } else {
          coords <- Seurat::GetTissueCoordinates(obj)
        }
        coords <- as.data.frame(coords)
        
        # Identify coordinate columns
        col_names <- colnames(coords)
        x_col <- if ("x" %in% col_names) "x" else if ("imagerow" %in% col_names) "imagerow" else "imagecol"
        y_col <- if ("y" %in% col_names) "y" else if ("imagecol" %in% col_names) "imagecol" else "imagerow"
        
        # Extract feature value
        feature <- lyr$feature
        if (is.null(feature)) {
          stop(paste("Layer", i, "must specify a 'feature' to plot."))
        }
        
        if (feature %in% colnames(obj@meta.data)) {
          value <- obj@meta.data[[feature]]
          names(value) <- rownames(obj@meta.data)
        } else {
          assay_name <- ifelse(is.null(lyr$assay), DefaultAssay(obj), lyr$assay)
          layer_name <- ifelse(is.null(lyr$layer), "data", lyr$layer)
          
          # In Seurat v5, we get assay data using GetAssayData
          exp_data <- GetAssayData(obj, assay = assay_name, layer = layer_name)
          if (!feature %in% rownames(exp_data)) {
            stop(paste("Feature", feature, "not found in metadata or assay", assay_name))
          }
          value <- exp_data[feature, ]
        }
        
        # Align barcodes
        cells <- rownames(coords)
        value <- value[cells]
        
        data.frame(
          cell = cells,
          x = coords[[x_col]],
          y = coords[[y_col]],
          value = value,
          stringsAsFactors = FALSE
        )
        
      } else if (lyr_type == "cellchat") {
        chat <- lyr$spatial_cellchat_obj
        if (is.null(chat)) {
          stop(paste("Layer", i, "is of type 'cellchat' but no 'spatial_cellchat_obj' was provided."))
        }
        
        pathway <- lyr$pathway
        pattern <- lyr$pattern
        if (is.null(pathway) || is.null(pattern)) {
          stop(paste("Layer", i, "is of type 'cellchat' and must specify both 'pathway' and 'pattern'."))
        }
        
        # Extract field data
        net0 <- methods::slot(chat, "netP")
        field_array <- net0$field[[pattern]]
        signaling_names <- if (!is.null(field_array)) dimnames(field_array)[[3]] else NULL
        
        # If the pathway has not been computed, try to compute it on the fly
        if (is.null(signaling_names) || !(pathway %in% signaling_names)) {
          if (requireNamespace("SpatialCellChat", quietly = TRUE)) {
            message(paste("Communication field for pathway", pathway, "not found in SpatialCellChat object. Attempting to compute it on the fly..."))
            
            chat_updated <- tryCatch({
              SpatialCellChat::computeCommunField(chat, signaling.name = pathway)
            }, error = function(e) {
              stop(paste("Could not compute communication field for pathway:", pathway, "-", e$message))
            })
            
            if (!is.null(chat_updated)) {
              chat <- chat_updated
              # Update local references
              net0 <- methods::slot(chat, "netP")
              field_array <- net0$field[[pattern]]
              signaling_names <- if (!is.null(field_array)) dimnames(field_array)[[3]] else NULL
            }
          } else {
            stop("SpatialCellChat package is not available to calculate communication fields on the fly.")
          }
        }
        
        if (is.null(signaling_names) || !(pathway %in% signaling_names)) {
          stop(paste("Pathway", pathway, "not computed in the SpatialCellChat object's", pattern, "fields."))
        }
        
        # Extract communication field vectors
        field_use <- data.frame(field_array[, , pathway])
        # Swapping matches SpatialCellChat::netVisual_CommunField swapping
        dx <- field_use[, 2]
        dy <- field_use[, 1]
        
        # Extract coordinates
        coordinates_df <- as.data.frame(chat@images$coordinates)
        col_names <- colnames(coordinates_df)
        
        if ("x_cent" %in% col_names && "y_cent" %in% col_names) {
          x_vals <- coordinates_df[["x_cent"]]
          y_vals <- coordinates_df[["y_cent"]]
        } else if ("x" %in% col_names && "y" %in% col_names) {
          x_vals <- coordinates_df[["x"]]
          y_vals <- coordinates_df[["y"]]
        } else {
          # Fallback to swapping if columns are in raw (y, x) order
          x_vals <- coordinates_df[, 2]
          y_vals <- coordinates_df[, 1]
        }
        cells <- rownames(coordinates_df)
        labels <- chat@idents
        
        data.frame(
          cell = cells,
          x = x_vals,
          y = y_vals,
          dx = dx,
          dy = dy,
          mag = sqrt(dx^2 + dy^2),
          value = labels,
          stringsAsFactors = FALSE
        )
        
      } else {
        stop(paste("Unsupported layer type:", lyr_type))
      }
    }, error = function(e) {
      warning(paste("Layer", i, "failed to compile and will be skipped:", e$message))
      NULL
    })
    
    if (!is.null(df)) {
      layer_dfs[[length(layer_dfs) + 1]] <- df
      valid_layers[[length(valid_layers) + 1]] <- lyr
    }
  }
  
  # Update layers list to match the valid layers plotted
  layers <- valid_layers
  
  if (length(layer_dfs) == 0) {
    stop("All specified layers failed to compile.")
  }
  
  # 3. Coordinate Normalization
  if (normalize_independently) {
    for (i in seq_along(layer_dfs)) {
      x_min <- min(layer_dfs[[i]]$x, na.rm = TRUE)
      x_max <- max(layer_dfs[[i]]$x, na.rm = TRUE)
      y_min <- min(layer_dfs[[i]]$y, na.rm = TRUE)
      y_max <- max(layer_dfs[[i]]$y, na.rm = TRUE)
      range_max <- max(x_max - x_min, y_max - y_min, na.rm = TRUE)
      
      layer_dfs[[i]]$x_norm <- (layer_dfs[[i]]$x - x_min) / range_max
      layer_dfs[[i]]$y_norm <- (layer_dfs[[i]]$y - y_min) / range_max
      
      if ("dx" %in% colnames(layer_dfs[[i]])) {
        layer_dfs[[i]]$dx_norm <- layer_dfs[[i]]$dx / range_max
        layer_dfs[[i]]$dy_norm <- layer_dfs[[i]]$dy / range_max
      }
    }
  } else {
    # Global normalization across all layers to preserve relative alignment
    all_x <- unlist(lapply(layer_dfs, function(d) d$x))
    all_y <- unlist(lapply(layer_dfs, function(d) d$y))
    x_min <- min(all_x, na.rm = TRUE)
    x_max <- max(all_x, na.rm = TRUE)
    y_min <- min(all_y, na.rm = TRUE)
    y_max <- max(all_y, na.rm = TRUE)
    range_max <- max(x_max - x_min, y_max - y_min, na.rm = TRUE)
    
    for (i in seq_along(layer_dfs)) {
      layer_dfs[[i]]$x_norm <- (layer_dfs[[i]]$x - x_min) / range_max
      layer_dfs[[i]]$y_norm <- (layer_dfs[[i]]$y - y_min) / range_max
      
      if ("dx" %in% colnames(layer_dfs[[i]])) {
        layer_dfs[[i]]$dx_norm <- layer_dfs[[i]]$dx / range_max
        layer_dfs[[i]]$dy_norm <- layer_dfs[[i]]$dy / range_max
      }
    }
  }
  
  # Apply Y reversion if requested
  if (reverse_y) {
    for (i in seq_along(layer_dfs)) {
      layer_dfs[[i]]$y_norm <- 1 - layer_dfs[[i]]$y_norm
    }
  }
  
  # 4. Perform Isometric 3D Projection
  # Order added: Layer 1 at bottom (z=0), Layer 2 at z=z_step, Layer 3 at z=2*z_step
  for (i in seq_along(layer_dfs)) {
    z_offset <- (i - 1) * z_step
    df <- layer_dfs[[i]]
    
    df$x_proj <- (df$x_norm - df$y_norm) * cos(theta_rad)
    df$y_proj <- (df$x_norm + df$y_norm) * sin(theta_rad) * sin(phi_rad) + z_offset
    
    if ("dx_norm" %in% colnames(df)) {
      # If y was reversed, the direction of dy is negated
      dy_val <- if (reverse_y) -df$dy_norm else df$dy_norm
      dx_val <- df$dx_norm
      
      # Project vectors
      df$dx_proj <- (dx_val - dy_val) * cos(theta_rad)
      df$dy_proj <- (dx_val + dy_val) * sin(theta_rad) * sin(phi_rad)
    }
    
    layer_dfs[[i]] <- df
  }
  
  # 5. Build ggplot
  p <- ggplot()
  
  # Prepare background platforms (using same pad range [0, 1])
  x_min_n <- -pad
  x_max_n <- 1 + pad
  y_min_n <- -pad
  y_max_n <- 1 + pad
  
  corners_base <- data.frame(
    x_n = c(x_min_n, x_max_n, x_max_n, x_min_n),
    y_n = c(y_min_n, y_min_n, y_max_n, y_max_n)
  )
  
  # Draw platforms and compute corner coordinates for structural lines
  platforms_corners <- list()
  for (i in seq_along(layer_dfs)) {
    z_offset <- (i - 1) * z_step
    
    corners_proj <- corners_base
    corners_proj$x_proj <- (corners_proj$x_n - corners_proj$y_n) * cos(theta_rad)
    corners_proj$y_proj <- (corners_proj$x_n + corners_proj$y_n) * sin(theta_rad) * sin(phi_rad) + z_offset
    
    platforms_corners[[i]] <- corners_proj
    
    # Render background platform
    p <- p + geom_polygon(
      data = corners_proj,
      aes(x = x_proj, y = y_proj),
      fill = "grey96",
      color = "grey80",
      linewidth = 0.5,
      alpha = 0.4
    )
  }
  
  # Draw vertical connecting cage lines if there are multiple layers
  if (length(layer_dfs) > 1) {
    bottom_corners <- platforms_corners[[1]]
    top_corners <- platforms_corners[[length(layer_dfs)]]
    
    for (j in 1:4) {
      cage_df <- data.frame(
        x = c(bottom_corners$x_proj[j], top_corners$x_proj[j]),
        y = c(bottom_corners$y_proj[j], top_corners$y_proj[j])
      )
      p <- p + geom_path(
        data = cage_df,
        aes(x = x, y = y),
        linetype = "dashed",
        color = "grey60",
        linewidth = 0.4
      )
    }
  }
  
  # Draw each layer's spots/vectors with its own legend scale
  for (i in seq_along(layer_dfs)) {
    df <- layer_dfs[[i]]
    lyr <- layers[[i]]
    lyr_type <- ifelse(is.null(lyr$type), "seurat", lyr$type)
    lyr_alpha <- ifelse(is.null(lyr$alpha), 0.8, lyr$alpha)
    
    # Legend title: User provided, or default "Layer I"
    leg_title <- lyr$legend_title
    if (is.null(leg_title)) {
      leg_title <- paste("Layer", i)
    }
    
    # Trigger ggnewscale for subsequent layers
    if (i > 1) {
      p <- p + ggnewscale::new_scale_color()
    }
    
    # Determine color scale mapping
    is_continuous <- is.numeric(df$value)
    
    # Draw spots
    p <- p + geom_point(
      data = df,
      aes(x = x_proj, y = y_proj, color = value),
      alpha = lyr_alpha,
      size = point_size
    )
    
    # Apply color scale
    if (is_continuous) {
      # Setup continuous color scale
      pal <- lyr$palette
      if (is.null(pal)) {
        # Default distinct continuous palettes per layer
        pal_options <- list(
          c("#440154", "#21908C", "#FDE725"), # Viridis
          c("#0B0010", "#7C1D6F", "#FCA50A"), # Inferno / Magma
          c("#03051A", "#3B4CC0", "#B40426")  # Blue-to-Red
        )
        pal <- pal_options[[((i - 1) %% 3) + 1]]
      }
      
      # Define colorbar guide (no override.aes needed, colorbar is opaque by default)
      guide_opt <- guide_colorbar(title = leg_title)
      
      if (length(pal) == 1 && is.character(pal)) {
        p <- p + scale_color_viridis_c(option = pal, name = leg_title, guide = guide_opt)
      } else {
        p <- p + scale_color_gradientn(colors = pal, name = leg_title, guide = guide_opt)
      }
      
    } else {
      # Setup discrete/categorical color scale
      df$value <- as.factor(df$value)
      levels_value <- levels(df$value)
      
      pal <- lyr$palette
      if (is.null(pal)) {
        # Auto distinct qualitative colors
        pal <- scales::hue_pal()(length(levels_value))
      } else if (length(pal) == 1 && is.character(pal)) {
        # Preset RColorBrewer palette name
        pal <- colorRampPalette(RColorBrewer::brewer.pal(9, pal))(length(levels_value))
      }
      
      # Define legend guide with override.aes to force full saturation and size 3
      guide_opt <- guide_legend(title = leg_title, override.aes = list(alpha = 1.0, size = 3.0))
      
      p <- p + scale_color_manual(values = pal, name = leg_title, guide = guide_opt)
    }
    
    # Draw communication field arrows if CellChat type
    if (lyr_type == "cellchat") {
      # Filter magnitude 0 displacements
      df_arrow <- df %>% filter(mag > 0)
      
      if (nrow(df_arrow) > 0) {
        # Scale arrows so the max magnitude spans ~0.05 normalized space
        max_mag <- max(df_arrow$mag, na.rm = TRUE)
        base_scale <- 0.05
        
        arrow_mult <- ifelse(is.null(lyr$arrow_scale_multiplier), 1.0, lyr$arrow_scale_multiplier)
        scale_factor <- (base_scale / max_mag) * arrow_mult
        
        df_arrow$x_end <- df_arrow$x_proj + df_arrow$dx_proj * scale_factor
        df_arrow$y_end <- df_arrow$y_proj + df_arrow$dy_proj * scale_factor
        
        arrow_col <- lyr$arrow_color
        if (is.null(arrow_col) || arrow_col == "map") {
          # Map arrow colors to cell type colors (shared scale)
          p <- p + geom_segment(
            data = df_arrow,
            aes(x = x_proj, y = y_proj, xend = x_end, yend = y_end, color = value),
            arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
            alpha = lyr_alpha,
            linewidth = 0.35,
            show.legend = FALSE # Spots provide the legend
          )
        } else {
          # Use static arrow color
          p <- p + geom_segment(
            data = df_arrow,
            aes(x = x_proj, y = y_proj, xend = x_end, yend = y_end),
            color = arrow_col,
            arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
            alpha = lyr_alpha,
            linewidth = 0.35
          )
        }
      }
    }
  }
  
  # 6. Apply Theme Styling
  p <- p +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
      legend.text = element_text(size = 8, family = "sans"),
      legend.title = element_text(face = "bold", size = 9, family = "sans"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_fixed()
  
  # 7. Save Files
  png_path <- file.path(output_dir, paste0(filename_base, ".png"))
  pdf_path <- file.path(output_dir, paste0(filename_base, ".pdf"))
  
  message("Saving PNG plot to: ", png_path)
  ggsave(png_path, plot = p, width = width, height = height, dpi = 300)
  
  message("Saving PDF plot to: ", pdf_path)
  ggsave(pdf_path, plot = p, width = width, height = height)
  
  return(p)
}
