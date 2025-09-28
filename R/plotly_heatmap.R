# Custom plotly heatmap functions for SeuratExplorer

#' Create an interactive plotly heatmap from Seurat data
#'
#' @param obj Seurat object
#' @param features Features to plot
#' @param group.by Grouping variable
#' @param assay Assay to use
#' @param slot Data slot to use
#' @param downsample Number of cells to downsample to
#' @param color_palette Color palette option
#' @return A plotly heatmap object
#' @export
plotly_heatmap <- function(obj, features, group.by = "seurat_clusters", 
                          assay = "RNA", slot = "data", 
                          downsample = 1000, color_palette = "default") {
  requireNamespace("plotly")
  requireNamespace("viridis")
  
  # Subset to selected features
  if (!all(features %in% rownames(obj[[assay]]))) {
    available_features <- features[features %in% rownames(obj[[assay]])]
    if (length(available_features) == 0) {
      stop("None of the requested features are available in the assay")
    }
    features <- available_features
  }
  
  # Downsample cells if needed
  if (ncol(obj) > downsample) {
    cells_to_use <- sample(colnames(obj), size = downsample)
    obj_subset <- obj[, cells_to_use]
  } else {
    obj_subset <- obj
  }
  
  # Get expression data
  if (slot == "data") {
    expr_data <- GetAssayData(obj_subset, assay = assay, layer = "data")
  } else if (slot == "counts") {
    expr_data <- GetAssayData(obj_subset, assay = assay, layer = "counts")
  } else {
    expr_data <- GetAssayData(obj_subset, assay = assay, layer = slot)
  }
  
  # Subset to features
  expr_data <- expr_data[features, , drop = FALSE]
  
  # Get cell metadata - fix the column access issue
  if (group.by %in% colnames(obj_subset@meta.data)) {
    cell_meta <- obj_subset@meta.data[, group.by, drop = FALSE]
  } else {
    # Use Seurat::Idents if group.by column doesn't exist
    cell_meta <- data.frame(cluster = Seurat::Idents(obj_subset))
    colnames(cell_meta) <- group.by
  }
  
  # Order cells by group
  cell_order <- order(cell_meta[, group.by])
  expr_data <- expr_data[, cell_order, drop = FALSE]
  cell_meta <- cell_meta[cell_order, , drop = FALSE]
  
  # Create plotly heatmap
  if (color_palette == "default") {
    colorscale <- list(c(0, "blue"), c(0.5, "white"), c(1, "red"))
  } else {
    # Use viridis palettes
    viridis_colors <- viridis::viridis_pal(option = color_palette)(256)
    colorscale <- lapply(seq(0, 1, length.out = 256), function(i) {
      list(i, viridis_colors[round(i * 255) + 1])
    })
  }
  
  p <- plotly::plot_ly(
    z = as.matrix(expr_data),
    type = "heatmap",
    colorscale = colorscale,
    showscale = TRUE,
    hovertemplate = paste(
      "Feature: %{y}<br>",
      "Cell: %{x}<br>",
      "Expression: %{z}<br>",
      "<extra></extra>"
    )
  ) %>%
    plotly::layout(
      title = paste("Expression Heatmap -", length(features), "features,", ncol(expr_data), "cells"),
      xaxis = list(title = "Cells", showticklabels = FALSE),
      yaxis = list(title = "Features", tickvals = seq_along(features) - 1, ticktext = features),
      margin = list(l = 150, r = 50, t = 50, b = 50)
    ) %>%
    plotly::config(
      displayModeBar = TRUE,
      modeBarButtonsToAdd = list("pan2d", "zoom2d", "resetScale2d"),
      scrollZoom = TRUE
    )
  
  return(p)
}

#' Create an interactive plotly heatmap from averaged expression
#'
#' @param obj Seurat object
#' @param features Features to plot
#' @param group.by Grouping variable
#' @param assay Assay to use
#' @param slot Data slot to use
#' @param color_palette Color palette option
#' @return A plotly heatmap object
#' @export
plotly_averaged_heatmap <- function(obj, features, group.by = "seurat_clusters", 
                                   assay = "RNA", slot = "data", 
                                   color_palette = "default") {
  requireNamespace("plotly")
  requireNamespace("viridis")
  
  # Get averaged expression
  avg_expr <- AverageExpression(
    obj, 
    features = features,
    group.by = group.by,
    assays = assay,
    layer = slot
  )
  
  # Extract the matrix
  if (length(avg_expr) == 1) {
    expr_matrix <- avg_expr[[1]]
  } else {
    expr_matrix <- avg_expr[[assay]]
  }
  
  # Create plotly heatmap
  if (color_palette == "default") {
    colorscale <- list(c(0, "blue"), c(0.5, "white"), c(1, "red"))
  } else {
    # Use viridis palettes
    viridis_colors <- viridis::viridis_pal(option = color_palette)(256)
    colorscale <- lapply(seq(0, 1, length.out = 256), function(i) {
      list(i, viridis_colors[round(i * 255) + 1])
    })
  }
  
  p <- plotly::plot_ly(
    z = as.matrix(expr_matrix),
    type = "heatmap",
    colorscale = colorscale,
    showscale = TRUE,
    hovertemplate = paste(
      "Feature: %{y}<br>",
      "Group: %{x}<br>",
      "Avg Expression: %{z}<br>",
      "<extra></extra>"
    )
  ) %>%
    plotly::layout(
      title = paste("Averaged Expression Heatmap -", nrow(expr_matrix), "features,", ncol(expr_matrix), "groups"),
      xaxis = list(title = "Groups"),
      yaxis = list(title = "Features"),
      margin = list(l = 150, r = 50, t = 50, b = 50)
    ) %>%
    plotly::config(
      displayModeBar = TRUE,
      modeBarButtonsToAdd = list("pan2d", "zoom2d", "resetScale2d"),
      scrollZoom = TRUE
    )
  
  return(p)
}