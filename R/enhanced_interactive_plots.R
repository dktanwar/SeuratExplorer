# Enhanced Interactive Plotting Functions for SeuratExplorer
# These functions provide significant improvements over the current implementation

#' Enhanced interactive DimPlot with custom tooltips and selection capabilities
#' @param obj Seurat object
#' @param reduction Dimensional reduction to use
#' @param group.by Metadata column for grouping
#' @param metadata_cols Additional metadata columns to show in tooltip
#' @param enable_selection Enable cell selection tools
#' @param downsample Downsample large datasets for performance
interactive_dimplot_enhanced <- function(obj,
                                         reduction = "umap",
                                         group.by = "seurat_clusters",
                                         metadata_cols = c("nFeature_RNA", "nCount_RNA"),
                                         enable_selection = TRUE,
                                         downsample = 10000,
                                         height = 500) {
  requireNamespace("plotly")
  requireNamespace("dplyr")

  # Downsample for performance if needed
  if (ncol(obj) > downsample) {
    set.seed(123)
    cells_sample <- sample(colnames(obj), downsample)
    obj_plot <- obj[, cells_sample]
    message(paste("Downsampling to", downsample, "cells for interactive display"))
  } else {
    obj_plot <- obj
  }

  # Get coordinates and metadata
  coords <- Embeddings(obj_plot, reduction = reduction)
  meta <- obj_plot@meta.data

  # Prepare tooltip data
  tooltip_cols <- c(group.by, metadata_cols)
  tooltip_data <- meta[, tooltip_cols, drop = FALSE]

  # Create hover text
  hover_text <- apply(tooltip_data, 1, function(x) {
    paste(names(x), ":", x, collapse = "<br>")
  })

  # Create the plot data
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    group = meta[[group.by]],
    hover_text = hover_text,
    cell_id = rownames(meta),
    stringsAsFactors = FALSE
  )

  # Create plotly object
  p <- plotly::plot_ly(
    data = plot_data,
    x = ~x,
    y = ~y,
    color = ~group,
    text = ~hover_text,
    hovertemplate = "%{text}<extra></extra>",
    type = "scatter",
    mode = "markers",
    marker = list(size = 3, opacity = 0.7),
    source = "dimplot" # Important for event handling
  ) %>%
    plotly::layout(
      title = list(
        text = paste("Interactive", toupper(reduction), "Plot"),
        font = list(size = 16)
      ),
      xaxis = list(
        title = paste0(toupper(reduction), "_1"),
        showgrid = TRUE,
        gridcolor = "rgba(128,128,128,0.2)"
      ),
      yaxis = list(
        title = paste0(toupper(reduction), "_2"),
        showgrid = TRUE,
        gridcolor = "rgba(128,128,128,0.2)"
      ),
      hovermode = "closest",
      dragmode = if (enable_selection) "select" else "pan",
      height = height
    ) %>%
    plotly::config(
      displayModeBar = TRUE,
      modeBarButtonsToAdd = if (enable_selection) list("select2d", "lasso2d") else NULL,
      displaylogo = FALSE,
      toImageButtonOptions = list(
        format = "png",
        filename = "interactive_dimplot",
        height = 800,
        width = 1000,
        scale = 2
      )
    )

  return(p)
}

#' Enhanced interactive FeaturePlot with expression-aware tooltips
#' @param obj Seurat object
#' @param features Features to plot
#' @param reduction Dimensional reduction
#' @param slot Data slot to use
#' @param assay Assay to use
interactive_featureplot_enhanced <- function(obj,
                                             features,
                                             reduction = "umap",
                                             slot = "data",
                                             assay = "RNA",
                                             downsample = 10000,
                                             height = 500) {
  requireNamespace("plotly")
  requireNamespace("viridis")

  # Downsample if needed
  if (ncol(obj) > downsample) {
    set.seed(123)
    cells_sample <- sample(colnames(obj), downsample)
    obj_plot <- obj[, cells_sample]
  } else {
    obj_plot <- obj
  }

  # Get coordinates
  coords <- Embeddings(obj_plot, reduction = reduction)

  # Get expression data
  DefaultAssay(obj_plot) <- assay
  expr_data <- FetchData(obj_plot, vars = features, slot = slot)

  # Handle single vs multiple features
  if (length(features) == 1) {
    plot_data <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      expression = expr_data[[features]],
      cell_id = rownames(expr_data),
      stringsAsFactors = FALSE
    )

    # Create hover text with expression info
    plot_data$hover_text <- paste0(
      "Cell: ", plot_data$cell_id, "<br>",
      features, ": ", round(plot_data$expression, 3), "<br>",
      "Coordinates: (", round(plot_data$x, 2), ", ", round(plot_data$y, 2), ")"
    )

    p <- plotly::plot_ly(
      data = plot_data,
      x = ~x,
      y = ~y,
      color = ~expression,
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      type = "scatter",
      mode = "markers",
      marker = list(size = 3),
      colors = viridis::viridis(256),
      source = "featureplot"
    ) %>%
      plotly::layout(
        title = list(text = paste("Interactive Feature Plot:", features)),
        xaxis = list(title = paste0(toupper(reduction), "_1")),
        yaxis = list(title = paste0(toupper(reduction), "_2")),
        height = height
      ) %>%
      plotly::colorbar(title = "Expression")
  } else {
    # Multiple features - create subplot
    plot_list <- list()

    for (i in seq_along(features)) {
      feature <- features[i]
      plot_data <- data.frame(
        x = coords[, 1],
        y = coords[, 2],
        expression = expr_data[[feature]],
        stringsAsFactors = FALSE
      )

      plot_data$hover_text <- paste0(
        feature, ": ", round(plot_data$expression, 3)
      )

      p_sub <- plotly::plot_ly(
        data = plot_data,
        x = ~x,
        y = ~y,
        color = ~expression,
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        type = "scatter",
        mode = "markers",
        marker = list(size = 2),
        colors = viridis::viridis(256),
        showscale = TRUE
      ) %>%
        plotly::layout(
          title = feature,
          xaxis = list(title = ""),
          yaxis = list(title = "")
        )

      plot_list[[i]] <- p_sub
    }

    # Create subplot
    n_features <- length(features)
    ncol <- ceiling(sqrt(n_features))
    nrow <- ceiling(n_features / ncol)

    p <- plotly::subplot(
      plot_list,
      nrows = nrow,
      ncols = ncol,
      titleX = TRUE,
      titleY = TRUE,
      margin = 0.05
    ) %>%
      plotly::layout(
        title = list(text = "Interactive Feature Plot"),
        height = height * nrow / 2
      )
  }

  return(p)
}

#' Create selection event handlers for interactive plots
#' @param input Shiny input object
#' @param output Shiny output object
#' @param obj Seurat object
create_selection_handlers <- function(input, output, obj) {
  # Store selected cells
  selected_cells <- reactiveValues(dimplot = NULL, featureplot = NULL)

  # Handle DimPlot selection
  observeEvent(plotly::event_data("plotly_selected", source = "dimplot"), {
    event <- plotly::event_data("plotly_selected", source = "dimplot")
    if (!is.null(event)) {
      # Get selected point indices
      selected_indices <- event$pointNumber + 1 # R indexing
      selected_cells$dimplot <- colnames(obj)[selected_indices]

      # Update UI to show selection info
      output$selection_info <- renderText({
        paste("Selected", length(selected_cells$dimplot), "cells")
      })
    }
  })

  # Handle FeaturePlot selection
  observeEvent(plotly::event_data("plotly_selected", source = "featureplot"), {
    event <- plotly::event_data("plotly_selected", source = "featureplot")
    if (!is.null(event)) {
      selected_indices <- event$pointNumber + 1
      selected_cells$featureplot <- colnames(obj)[selected_indices]
    }
  })

  # Download handler for selected cells
  output$download_selected_cells <- downloadHandler(
    filename = function() {
      "selected_cells.csv"
    },
    content = function(file) {
      if (!is.null(selected_cells$dimplot)) {
        cells_df <- data.frame(
          cell_barcode = selected_cells$dimplot,
          stringsAsFactors = FALSE
        )
        write.csv(cells_df, file, row.names = FALSE)
      }
    }
  )

  return(selected_cells)
}

#' Performance-optimized plotting for large datasets
#' @param obj Seurat object
#' @param max_cells Maximum cells to display
#' @param sampling_method Method for downsampling ('random', 'stratified')
optimize_for_performance <- function(obj,
                                     max_cells = 10000,
                                     sampling_method = "stratified") {
  if (ncol(obj) <= max_cells) {
    return(obj)
  }

  if (sampling_method == "random") {
    # Simple random sampling
    set.seed(123)
    cells_keep <- sample(colnames(obj), max_cells)
  } else if (sampling_method == "stratified") {
    # Stratified sampling by cluster to maintain representation
    clusters <- obj$seurat_clusters
    cells_per_cluster <- floor(max_cells / length(levels(clusters)))

    cells_keep <- c()
    for (cluster in levels(clusters)) {
      cluster_cells <- colnames(obj)[clusters == cluster]
      if (length(cluster_cells) > cells_per_cluster) {
        sampled_cells <- sample(
          cluster_cells,
          min(cells_per_cluster, length(cluster_cells))
        )
      } else {
        sampled_cells <- cluster_cells
      }
      cells_keep <- c(cells_keep, sampled_cells)
    }
  }

  return(obj[, cells_keep])
}

#' Add cross-plot linking functionality
#' @param plot_outputs List of plotly outputs to link
create_linked_plots <- function(plot_outputs) {
  # This would be implemented in the server logic
  # When selection changes in one plot, highlight same cells in others

  observe({
    # Get selection from any source
    selection <- plotly::event_data("plotly_selected")

    if (!is.null(selection)) {
      # Update all linked plots to highlight selected cells
      for (plot_id in names(plot_outputs)) {
        plotly::plotlyProxy(plot_id) %>%
          plotly::plotlyProxyInvoke(
            "restyle",
            list(opacity = c(0.3, 1.0)),
            list(selection$pointNumber)
          )
      }
    }
  })
}

# Example usage in server.R:
#
# # Enhanced interactive DimPlot
# output$dimplot_interactive <- plotly::renderPlotly({
#   interactive_dimplot_enhanced(
#     obj = data$obj,
#     reduction = input$DimDimensionReduction,
#     group.by = input$DimClusterResolution,
#     metadata_cols = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#     enable_selection = TRUE
#   )
# })
#
# # Selection handlers
# selected_cells <- create_selection_handlers(input, output, data$obj)


#' Fixed interactive DimPlot with proper tooltips and height control
#' @param p ggplot object from Seurat::DimPlot
#' @param obj Seurat object
#' @param reduction Reduction type
#' @param group.by Grouping column
#' @param hw_ratio Height/width ratio
#' @param width_px Width in pixels
#' @param tooltip_vars Variables to show in tooltip
#' @param theme_choice Theme for the plot
interactive_dimplot_fixed <- function(p, obj, reduction = "umap", group.by = "seurat_clusters",
                                      hw_ratio = 0.9, width_px = 600,
                                      tooltip_vars = c("nFeature_RNA", "nCount_RNA"),
                                      theme_choice = "plotly_white") {
  requireNamespace("plotly")

  # Calculate height from width and ratio
  height_px <- width_px * hw_ratio

  # Get coordinates and metadata for custom tooltips
  coords <- Embeddings(obj, reduction = reduction)
  meta <- obj@meta.data

  # Create meaningful tooltip text
  if (!is.null(tooltip_vars) && all(tooltip_vars %in% colnames(meta))) {
    hover_text <- apply(meta[, c(group.by, tooltip_vars), drop = FALSE], 1, function(x) {
      paste(paste(names(x), x, sep = ": "), collapse = "<br>")
    })
  } else {
    hover_text <- paste("Cluster:", meta[[group.by]])
  }

  # Extract plot data from ggplot
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]

  # Create plotly from scratch with width and height in plot_ly
  plotly_obj <- plotly::plot_ly(
    x = coords[, 1],
    y = coords[, 2],
    color = meta[[group.by]],
    colors = unique(plot_data$colour),
    text = hover_text,
    hovertemplate = "%{text}<extra></extra>",
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.7),
    width = width_px,  # Moved here
    height = height_px # Moved here
  ) %>%
    plotly::layout(
      title = list(text = paste("Interactive", toupper(reduction)), font = list(size = 16)),
      xaxis = list(
        title = paste0(toupper(reduction), "_1"),
        showgrid = TRUE,
        gridcolor = "rgba(128,128,128,0.2)"
      ),
      yaxis = list(
        title = paste0(toupper(reduction), "_2"),
        showgrid = TRUE,
        gridcolor = "rgba(128,128,128,0.2)"
      ),
      template = theme_choice,
      legend = list(title = list(text = group.by))
    ) %>%
    plotly::config(
      displayModeBar = TRUE,
      toImageButtonOptions = list(
        format = "png",
        filename = "interactive_dimplot",
        height = height_px,
        width = width_px,
        scale = 2
      )
    )

  return(plotly_obj)
}

