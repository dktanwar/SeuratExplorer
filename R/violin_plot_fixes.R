# Violin Plot Fixes for SeuratExplorer
# 1. Fix interactive violin plot to properly handle multiple genes (separate subplots)
# 2. Add option to create separate plots by condition instead of just side-by-side comparison

#' Enhanced interactive violin plot that properly handles multiple genes
#' @param p ggplot object from Seurat::VlnPlot
#' @param obj Seurat object
#' @param features Features being plotted
#' @param group.by Grouping column
#' @param slot Data slot
#' @param split_by Split by variable
#' @param separate_conditions Create separate plots by condition
#' @param hw_ratio Height/width ratio
#' @param width_px Width in pixels
interactive_vlnplot_enhanced <- function(p, obj, features, group.by = "seurat_clusters",
                                         slot = "data", split_by = NULL,
                                         separate_conditions = FALSE,
                                         hw_ratio = 0.9, width_px = 600,
                                         theme_choice = "plotly_white") {
  requireNamespace("plotly")
  requireNamespace("dplyr")

  height_px <- width_px * hw_ratio

  # Handle multiple genes case
  if (length(features) > 1 && !separate_conditions) {
    # Create separate subplot for each gene (like static version)
    plot_list <- list()

    for (i in seq_along(features)) {
      feature <- features[i]

      # Create individual violin plot for this feature
      p_single <- Seurat::VlnPlot(obj,
        features = feature,
        group.by = group.by,
        split.by = split_by,
        pt.size = 0
      ) # Set pt.size to 0 for cleaner look

      # Convert to plotly
      p_plotly <- plotly::ggplotly(p_single) %>%
        plotly::layout(
          title = feature,
          showlegend = (i == 1), # Only show legend for first plot
          template = theme_choice
        )

      plot_list[[i]] <- p_plotly
    }

    # Create subplot layout
    n_features <- length(features)
    ncol_subplot <- ceiling(sqrt(n_features))
    nrow_subplot <- ceiling(n_features / ncol_subplot)

    combined_plot <- plotly::subplot(
      plot_list,
      nrows = nrow_subplot,
      ncols = ncol_subplot,
      subplot_titles = features,
      shareY = FALSE,
      titleX = TRUE,
      titleY = TRUE,
      margin = 0.08
    ) %>%
      plotly::layout(
        title = list(text = paste("Interactive Violin Plot:", paste(features, collapse = ", "))),
        width = width_px,
        height = height_px * max(1, nrow_subplot * 0.8),
        template = theme_choice
      )

    return(combined_plot)
  } else if (separate_conditions && !is.null(split_by)) {
    # Create separate plots for each condition
    conditions <- unique(obj@meta.data[[split_by]])
    plot_list <- list()

    for (i in seq_along(conditions)) {
      condition <- conditions[i]

      # Subset object for this condition
      obj_subset <- subset(obj, subset = obj@meta.data[[split_by]] == condition)

      # Create violin plot for this condition
      if (length(features) == 1) {
        p_condition <- Seurat::VlnPlot(obj_subset,
          features = features,
          group.by = group.by,
          pt.size = 0
        )
      } else {
        # Multiple features, single condition - create combined plot
        p_condition <- Seurat::VlnPlot(obj_subset,
          features = features,
          group.by = group.by,
          pt.size = 0,
          ncol = ceiling(sqrt(length(features)))
        )
      }

      # Convert to plotly
      p_plotly <- plotly::ggplotly(p_condition) %>%
        plotly::layout(
          title = paste("Condition:", condition),
          template = theme_choice
        )

      plot_list[[i]] <- p_plotly
    }

    # Arrange condition plots
    n_conditions <- length(conditions)
    if (n_conditions == 2) {
      # Side by side for 2 conditions
      combined_plot <- plotly::subplot(
        plot_list,
        nrows = 1,
        ncols = 2,
        subplot_titles = paste("Condition:", conditions),
        shareY = TRUE,
        titleX = TRUE,
        margin = 0.08
      )
    } else {
      # Grid layout for more conditions
      ncol_subplot <- ceiling(sqrt(n_conditions))
      nrow_subplot <- ceiling(n_conditions / ncol_subplot)

      combined_plot <- plotly::subplot(
        plot_list,
        nrows = nrow_subplot,
        ncols = ncol_subplot,
        subplot_titles = paste("Condition:", conditions),
        shareY = FALSE,
        titleX = TRUE,
        margin = 0.08
      )
    }

    combined_plot <- combined_plot %>%
      plotly::layout(
        title = list(text = "Interactive Violin Plot - Separate Conditions"),
        width = width_px,
        height = height_px * max(1, length(conditions) * 0.6),
        template = theme_choice
      )

    return(combined_plot)
  } else {
    # Single gene or standard split plot - use original approach
    plotly_obj <- plotly::ggplotly(p, tooltip = "all") %>%
      plotly::layout(
        title = list(text = paste("Interactive Violin Plot:", paste(features, collapse = ", "))),
        width = width_px,
        height = height_px,
        template = theme_choice
      )

    return(plotly_obj)
  }
}

# =============================================================================
# UI MODIFICATIONS - Add to violin plot settings box
# =============================================================================

#' UI additions for enhanced violin plot controls
ui_violin_plot_enhancements <- function() {
  tagList(
    # Add option for separate condition plots
    conditionalPanel(
      condition = "input.VlnSplitBy != 'None'",
      checkboxInput("VlnSeparateConditions",
        "Create separate plots by condition",
        value = FALSE
      ),
      helpText("When checked, creates separate violin plots for each condition instead of side-by-side comparison",
        style = "font-size: 11px; color: #666;"
      )
    ),

    # Add column layout control for multiple genes
    conditionalPanel(
      condition = "output.Vlnplot_multiple_genes",
      sliderInput("VlnPlotNCols", "Subplot Columns:",
        min = 1, max = 4, value = 2, step = 1
      ),
      helpText("Number of columns for multiple gene subplots",
        style = "font-size: 11px; color: #666;"
      )
    ),

    # Interactive plot specific controls
    conditionalPanel(
      condition = "input.vlnplot_mode == 'interactive'",
      selectInput("VlnPlotTheme", "Interactive Theme:",
        choices = c(
          "White" = "plotly_white",
          "Dark" = "plotly_dark",
          "ggplot2" = "ggplot2",
          "Simple" = "simple_white"
        ),
        selected = "plotly_white"
      )
    )
  )
}

# =============================================================================
# SERVER MODIFICATIONS - Replace existing violin plot interactive output
# =============================================================================

#' Fixed server function for interactive violin plots
server_vlnplot_enhanced <- function() {
  # Enhanced Interactive VlnPlot output
  output$vlnplot_interactive <- plotly::renderPlotly({
    if (verbose) {
      message("SeuratExplorer: preparing enhanced vlnplot_interactive...")
    }

    if (any(is.na(features_vlnplot$features_current))) {
      interactive_empty_plot()
    } else {
      cds <- data$obj
      cds@meta.data[, input$VlnClusterResolution] <- factor(cds@meta.data[, input$VlnClusterResolution],
        levels = input$VlnClusterOrder
      )
      SeuratObject::Idents(cds) <- input$VlnClusterResolution

      if (!any(features_vlnplot$features_current %in% c(rownames(cds[[input$VlnAssay]]), data$extra_qc_options))) {
        interactive_empty_plot()
      } else {
        # Get split by condition (revised)
        split_by_var <- if (is.null(VlnSplit.Revised()) || VlnSplit.Revised() == "None") {
          NULL
        } else {
          VlnSplit.Revised()
        }

        # Determine separation mode
        separate_conditions <- !is.null(input$VlnSeparateConditions) && input$VlnSeparateConditions

        if (length(features_vlnplot$features_current) == 1) {
          # Single Gene
          p <- Seurat::VlnPlot(cds,
            features = features_vlnplot$features_current,
            assay = input$VlnAssay,
            layer = input$VlnSlot,
            split.by = split_by_var,
            split.plot = input$VlnSplitPlot,
            pt.size = input$VlnPointSize,
            alpha = input$VlnPointAlpha,
            idents = input$VlnIdentsSelected
          ) &
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
              axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize)
            )
        } else {
          # Multiple genes - different handling for static vs interactive
          if (!separate_conditions) {
            # For multiple genes, create individual plots that will be combined as subplots
            p <- Seurat::VlnPlot(cds,
              features = features_vlnplot$features_current[1], # Just use first gene for base plot
              assay = input$VlnAssay,
              layer = input$VlnSlot,
              split.by = split_by_var,
              split.plot = input$VlnSplitPlot,
              pt.size = input$VlnPointSize,
              alpha = input$VlnPointAlpha,
              idents = input$VlnIdentsSelected
            ) &
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize)
              )
          } else {
            # For separate conditions, create base plot
            p <- Seurat::VlnPlot(cds,
              features = features_vlnplot$features_current,
              assay = input$VlnAssay,
              layer = input$VlnSlot,
              split.by = NULL, # No split for separate conditions
              pt.size = input$VlnPointSize,
              alpha = input$VlnPointAlpha,
              idents = input$VlnIdentsSelected,
              ncol = input$VlnPlotNCols %||% 2
            ) &
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize)
              )
          }
        }

        # Apply color palette
        if (!is.null(input$Vlnfillcolorplatte) && input$Vlnfillcolorplatte != "default" &&
          (is.null(split_by_var) || !separate_conditions)) {
          fill.colors <- getColors(
            color.platte = color_list,
            choice = input$Vlnfillcolorplatte,
            n = length(levels(Idents(cds)))
          )
          names(fill.colors) <- levels(Idents(cds))
          p <- p & scale_fill_manual(values = fill.colors)
        }

        # Convert to interactive plot with enhanced features
        interactive_vlnplot_enhanced(
          p = p,
          obj = cds,
          features = features_vlnplot$features_current,
          group.by = input$VlnClusterResolution,
          slot = input$VlnSlot,
          split_by = split_by_var,
          separate_conditions = separate_conditions,
          hw_ratio = input$VlnPlotHWRatio,
          width_px = session$clientData$output_vlnplot_width %||% 600,
          theme_choice = input$VlnPlotTheme %||% "plotly_white"
        )
      }
    }
  })
}

# =============================================================================
# INTEGRATION INSTRUCTIONS
# =============================================================================

# To integrate these fixes:

# 1. Add the ui_violin_plot_enhancements() content to the violin plot settings box in ui.R
#    Insert after line ~207 in the existing violin plot box

# 2. Replace the existing vlnplot_interactive output in server.R with server_vlnplot_enhanced()

# 3. Add the interactive_vlnplot_enhanced function to functions.R

# Example integration for ui.R:
# In the violin plot tab_list[["vlnplot"]] section, after the existing controls, add:

# # Enhanced controls
# ui_violin_plot_enhancements()

# Example integration for server.R:
# Replace the existing output$vlnplot_interactive <- ... block with:
# server_vlnplot_enhanced()

# =============================================================================
# ADDITIONAL HELPER FUNCTIONS
# =============================================================================

#' Create violin plot with proper subplot layout for static plots
create_static_vlnplot_with_layout <- function(obj, features, group.by, split.by = NULL,
                                              separate_conditions = FALSE, ncol = NULL, ...) {
  if (separate_conditions && !is.null(split.by)) {
    # Create separate plots for each condition
    conditions <- unique(obj@meta.data[[split.by]])
    plot_list <- list()

    for (condition in conditions) {
      obj_subset <- subset(obj, subset = obj@meta.data[[split.by]] == condition)

      p_condition <- Seurat::VlnPlot(obj_subset,
        features = features,
        group.by = group.by,
        ncol = ncol,
        ...
      ) +
        ggtitle(paste("Condition:", condition))

      plot_list[[length(plot_list) + 1]] <- p_condition
    }

    # Combine plots
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = length(conditions))
    return(combined_plot)
  } else {
    # Standard violin plot
    return(Seurat::VlnPlot(obj,
      features = features, group.by = group.by,
      split.by = split.by, ncol = ncol, ...
    ))
  }
}
