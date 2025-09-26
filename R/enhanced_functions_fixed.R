# Enhanced Interactive Plotting Functions - FIXES for SeuratExplorer Issues
# This addresses all the issues you mentioned

# =============================================================================
# ISSUE 1: Interactive plots don't respect H/W ratio and proper height
# =============================================================================

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

  # Create plotly from scratch instead of converting ggplot
  plotly_obj <- plotly::plot_ly(
    x = coords[, 1],
    y = coords[, 2],
    color = meta[[group.by]],
    colors = unique(plot_data$colour),
    text = hover_text,
    hovertemplate = "%{text}<extra></extra>",
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.7)
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
      width = width_px,
      height = height_px,
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

#' Fixed interactive FeaturePlot with proper height control and meaningful tooltips
interactive_featureplot_fixed <- function(p, obj, features, reduction = "umap", slot = "data",
                                          hw_ratio = 0.9, width_px = 600,
                                          theme_choice = "plotly_white") {
  requireNamespace("plotly")

  height_px <- width_px * hw_ratio
  coords <- Embeddings(obj, reduction = reduction)

  # Get expression data
  expr_data <- FetchData(obj, vars = features, slot = slot)

  if (length(features) == 1) {
    feature <- features[1]

    # Create meaningful hover text with expression values
    hover_text <- paste0(
      "Cell: ", rownames(coords), "<br>",
      feature, ": ", round(expr_data[[feature]], 3), "<br>",
      "Coordinates: (", round(coords[, 1], 2), ", ", round(coords[, 2], 2), ")"
    )

    plotly_obj <- plotly::plot_ly(
      x = coords[, 1],
      y = coords[, 2],
      color = expr_data[[feature]],
      text = hover_text,
      hovertemplate = "%{text}<extra></extra>",
      type = "scatter",
      mode = "markers",
      marker = list(size = 4, opacity = 0.8),
      colorscale = "Viridis"
    ) %>%
      plotly::layout(
        title = list(text = paste("Interactive Feature Plot:", feature)),
        xaxis = list(title = paste0(toupper(reduction), "_1")),
        yaxis = list(title = paste0(toupper(reduction), "_2")),
        width = width_px,
        height = height_px,
        template = theme_choice
      ) %>%
      plotly::colorbar(title = "Expression")
  } else {
    # Multiple features - create subplots
    plot_list <- list()
    n_features <- length(features)
    ncol_subplot <- ceiling(sqrt(n_features))
    nrow_subplot <- ceiling(n_features / ncol_subplot)

    for (i in seq_along(features)) {
      feature <- features[i]
      hover_text <- paste0(feature, ": ", round(expr_data[[feature]], 3))

      p_sub <- plotly::plot_ly(
        x = coords[, 1],
        y = coords[, 2],
        color = expr_data[[feature]],
        text = hover_text,
        hovertemplate = "%{text}<extra></extra>",
        type = "scatter",
        mode = "markers",
        marker = list(size = 3),
        colorscale = "Viridis",
        showscale = (i == 1)
      )
      plot_list[[i]] <- p_sub
    }

    plotly_obj <- plotly::subplot(
      plot_list,
      nrows = nrow_subplot,
      ncols = ncol_subplot,
      subplot_titles = features
    ) %>%
      plotly::layout(
        title = list(text = "Interactive Feature Plot"),
        width = width_px,
        height = height_px * nrow_subplot,
        template = theme_choice
      )
  }

  return(plotly_obj)
}

#' Interactive VlnPlot with proper height control
interactive_vlnplot_fixed <- function(p, obj, features, group.by = "seurat_clusters",
                                      hw_ratio = 0.9, width_px = 600,
                                      theme_choice = "plotly_white") {
  requireNamespace("plotly")

  height_px <- width_px * hw_ratio

  # Convert ggplot to plotly with custom height
  plotly_obj <- plotly::ggplotly(p, tooltip = "all") %>%
    plotly::layout(
      title = list(text = paste("Interactive Violin Plot:", paste(features, collapse = ", "))),
      width = width_px,
      height = height_px,
      template = theme_choice
    )

  return(plotly_obj)
}

# =============================================================================
# ISSUE 2: Interactive Ridge Plot
# =============================================================================

#' Interactive Ridge Plot
interactive_ridgeplot <- function(p, obj, features, group.by = "seurat_clusters",
                                  hw_ratio = 0.9, width_px = 600,
                                  theme_choice = "plotly_white") {
  requireNamespace("plotly")

  height_px <- width_px * hw_ratio

  plotly_obj <- plotly::ggplotly(p, tooltip = "all") %>%
    plotly::layout(
      title = list(text = paste("Interactive Ridge Plot:", paste(features, collapse = ", "))),
      width = width_px,
      height = height_px,
      template = theme_choice
    )

  return(plotly_obj)
}

# =============================================================================
# ISSUE 3: Interactive Complex Heatmap using InteractiveComplexHeatmap
# =============================================================================

#' Interactive Complex Heatmap with viridis colors
#' @param obj Seurat object
#' @param features Features to plot
#' @param group.by Grouping variable
#' @param color_palette Color palette choice
interactive_complex_heatmap <- function(obj, features, group.by = "seurat_clusters",
                                        color_palette = "viridis",
                                        cluster_columns = TRUE, cluster_rows = TRUE) {
  requireNamespace("ComplexHeatmap")
  requireNamespace("InteractiveComplexHeatmap")
  requireNamespace("viridis")

  # Get average expression data
  avg_expr <- AverageExpression(obj,
    features = features, group.by = group.by,
    assays = DefaultAssay(obj), slot = "data"
  )

  expr_matrix <- as.matrix(avg_expr[[DefaultAssay(obj)]])

  # Scale the data
  expr_scaled <- t(scale(t(expr_matrix)))

  # Choose color palette
  color_fun <- switch(color_palette,
    "viridis" = circlize::colorRamp2(c(-2, 0, 2), viridis::viridis(3)),
    "plasma" = circlize::colorRamp2(c(-2, 0, 2), viridis::plasma(3)),
    "inferno" = circlize::colorRamp2(c(-2, 0, 2), viridis::inferno(3)),
    "magma" = circlize::colorRamp2(c(-2, 0, 2), viridis::magma(3)),
    "cividis" = circlize::colorRamp2(c(-2, 0, 2), viridis::cividis(3)),
    "default" = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )

  # Create the heatmap
  ht <- ComplexHeatmap::Heatmap(
    expr_scaled,
    name = "Z-score",
    col = color_fun,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 12),
    heatmap_legend_param = list(title = "Expression Z-score")
  )

  # Make it interactive
  InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput("heatmap")

  return(ht)
}

# =============================================================================
# ISSUE 4: Enhanced Violin Plot Layout (like your example)
# =============================================================================

#' Enhanced VlnPlot with custom layout
create_enhanced_vlnplot <- function(obj, features, group.by = "samples",
                                    pt.size = 0, ncol = 3, theme_choice = "classic") {
  requireNamespace("Seurat")
  requireNamespace("ggplot2")

  p <- Seurat::VlnPlot(
    obj,
    features = features,
    group.by = group.by,
    pt.size = pt.size,
    ncol = ncol
  )

  # Apply theme
  theme_fun <- switch(theme_choice,
    "classic" = ggplot2::theme_classic(),
    "minimal" = ggplot2::theme_minimal(),
    "bw" = ggplot2::theme_bw(),
    "void" = ggplot2::theme_void(),
    "dark" = ggplot2::theme_dark(),
    ggplot2::theme_classic()
  )

  p <- p & theme_fun

  return(p)
}

# =============================================================================
# ISSUE 5: Feature Selection and Import System
# =============================================================================

#' Feature Selection UI Module
feature_selection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Selected Features"),
    verbatimTextOutput(ns("selected_features")),
    br(),
    fluidRow(
      column(4, actionButton(ns("import_to_violin"), "Import to Violin Plot",
        class = "btn-primary btn-sm"
      )),
      column(4, actionButton(ns("import_to_dot"), "Import to Dot Plot",
        class = "btn-primary btn-sm"
      )),
      column(4, actionButton(ns("import_to_feature"), "Import to Feature Plot",
        class = "btn-primary btn-sm"
      ))
    )
  )
}

#' Feature Selection Server Module
feature_selection_server <- function(id, selected_features_reactive) {
  moduleServer(id, function(input, output, session) {
    output$selected_features <- renderText({
      features <- selected_features_reactive()
      if (is.null(features) || length(features) == 0) {
        "No features selected"
      } else {
        paste(features, collapse = ", ")
      }
    })

    # Return reactive values for button clicks
    list(
      import_to_violin = reactive(input$import_to_violin),
      import_to_dot = reactive(input$import_to_dot),
      import_to_feature = reactive(input$import_to_feature),
      selected_features = selected_features_reactive
    )
  })
}

# =============================================================================
# ISSUE 6: Theme Selection System
# =============================================================================

#' Theme Selection UI
theme_selection_ui <- function(id, for_interactive = TRUE) {
  ns <- NS(id)

  if (for_interactive) {
    # Plotly themes
    choices <- c(
      "plotly_white" = "plotly_white",
      "plotly_dark" = "plotly_dark",
      "ggplot2" = "ggplot2",
      "seaborn" = "seaborn",
      "simple_white" = "simple_white",
      "none" = "none"
    )
  } else {
    # ggplot2 themes
    choices <- c(
      "classic" = "classic",
      "minimal" = "minimal",
      "bw" = "bw",
      "void" = "void",
      "dark" = "dark",
      "light" = "light"
    )
  }

  selectInput(ns("theme"), "Plot Theme:", choices = choices, selected = choices[1])
}

# =============================================================================
# ISSUE 7: Tooltip Configuration System
# =============================================================================

#' Tooltip Configuration UI
tooltip_config_ui <- function(id, available_vars) {
  ns <- NS(id)

  tagList(
    h5("Hover Information"),
    checkboxGroupInput(
      ns("tooltip_vars"),
      "Variables to show in hover:",
      choices = available_vars,
      selected = available_vars[1:min(3, length(available_vars))]
    )
  )
}

# =============================================================================
# EXAMPLE USAGE IN SERVER.R - HOW TO INTEGRATE THESE FIXES
# =============================================================================

# In your server.R, replace the interactive plot outputs like this:

# # Fixed Interactive DimPlot with proper height control
# output$dimplot_interactive <- plotly::renderPlotly({
#   req(input$DimDimensionReduction)
#
#   cds <- data$obj
#   Seurat::Idents(cds) <- input$DimClusterResolution
#
#   # Generate base plot (same logic as static)
#   if(is.null(input$DimHighlightedClusters)) {
#     dim_cells_highlighted <- NULL
#   } else {
#     dim_cells_highlighted <- colnames(cds)[cds@meta.data[,input$DimClusterResolution] %in% input$DimHighlightedClusters]
#   }
#
#   cds@meta.data[,input$DimClusterResolution] <- factor(cds@meta.data[,input$DimClusterResolution],
#                                                        levels = input$DimClusterOrder)
#
#   p <- Seurat::DimPlot(cds,
#                        reduction = input$DimDimensionReduction,
#                        label = input$DimShowLabel,
#                        pt.size = input$DimPointSize,
#                        label.size = input$DimLabelSize,
#                        group.by = input$DimClusterResolution,
#                        cells.highlight = dim_cells_highlighted)
#
#   if(!input$DimShowLegend) {
#     p <- p & NoLegend()
#   }
#
#   # Use fixed interactive function with proper height
#   interactive_dimplot_fixed(
#     p = p,
#     obj = cds,
#     reduction = input$DimDimensionReduction,
#     group.by = input$DimClusterResolution,
#     hw_ratio = input$DimPlotHWRatio,
#     width_px = session$clientData$output_dimplot_width,
#     tooltip_vars = input$tooltip_vars,  # from tooltip config UI
#     theme_choice = input$plot_theme     # from theme selection UI
#   )
# })

# =============================================================================
# ADDITIONAL HELPER FUNCTIONS
# =============================================================================

#' Update gene input from feature selection
update_gene_input <- function(session, input_id, new_genes) {
  updateTextAreaInput(session, input_id, value = paste(new_genes, collapse = "\n"))
}

#' Color palette choices for different plot types
get_color_palette_choices <- function(plot_type = "heatmap") {
  switch(plot_type,
    "heatmap" = c(
      "Viridis" = "viridis",
      "Plasma" = "plasma",
      "Inferno" = "inferno",
      "Magma" = "magma",
      "Cividis" = "cividis",
      "Default" = "default"
    ),
    "dimplot" = c(
      "Default" = "default",
      "Set1" = "Set1",
      "Set2" = "Set2",
      "Dark2" = "Dark2",
      "Paired" = "Paired"
    )
  )
}
