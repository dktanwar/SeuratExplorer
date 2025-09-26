# UI and Server Integration Examples
# How to modify existing UI and server code to implement the fixes

# =============================================================================
# 1. MODIFY UI.R - Add theme selection and tooltip configuration
# =============================================================================

# Add these UI elements to the existing boxes in ui.R:

# For DimPlot box - add after existing controls:
ui_dimplot_additions <- function() {
  tagList(
    # Theme selection for interactive plots
    conditionalPanel(
      condition = "input.dimplot_mode == 'interactive'",
      selectInput("DimPlotTheme", "Interactive Theme:",
                  choices = c("White" = "plotly_white",
                             "Dark" = "plotly_dark",
                             "ggplot2" = "ggplot2",
                             "Seaborn" = "seaborn",
                             "Simple" = "simple_white"),
                  selected = "plotly_white")
    ),

    # Theme selection for static plots
    conditionalPanel(
      condition = "input.dimplot_mode == 'static'",
      selectInput("DimPlotThemeStatic", "Static Theme:",
                  choices = c("Classic" = "classic",
                             "Minimal" = "minimal",
                             "BW" = "bw",
                             "Void" = "void",
                             "Dark" = "dark"),
                  selected = "classic")
    ),

    # Tooltip configuration for interactive plots
    conditionalPanel(
      condition = "input.dimplot_mode == 'interactive'",
      tags$hr(),
      h5("Hover Information"),
      checkboxGroupInput("DimPlotTooltipVars",
                        "Show in hover:",
                        choices = c("Cluster" = "cluster",
                                   "Feature Count" = "nFeature_RNA",
                                   "UMI Count" = "nCount_RNA",
                                   "Mitochondrial %" = "percent.mt"),
                        selected = c("cluster", "nFeature_RNA", "nCount_RNA"))
    )
  )
}

# For FeaturePlot box - similar additions
ui_featureplot_additions <- function() {
  tagList(
    conditionalPanel(
      condition = "input.featureplot_mode == 'interactive'",
      selectInput("FeaturePlotTheme", "Interactive Theme:",
                  choices = c("White" = "plotly_white",
                             "Dark" = "plotly_dark",
                             "ggplot2" = "ggplot2"),
                  selected = "plotly_white")
    ),
    conditionalPanel(
      condition = "input.featureplot_mode == 'static'",
      selectInput("FeaturePlotThemeStatic", "Static Theme:",
                  choices = c("Classic" = "classic",
                             "Minimal" = "minimal",
                             "BW" = "bw"),
                  selected = "classic")
    )
  )
}

# For Heatmap - add viridis color palette options
ui_heatmap_additions <- function() {
  tagList(
    # Add after existing heatmap controls
    selectInput("HeatmapColorPalette", "Color Palette:",
                choices = c("Viridis" = "viridis",
                           "Plasma" = "plasma",
                           "Inferno" = "inferno",
                           "Magma" = "magma",
                           "Cividis" = "cividis",
                           "Default" = "default"),
                selected = "viridis"),

    # Interactive heatmap option
    div(style = "margin-bottom: 10px;",
        shinyWidgets::radioGroupButtons(
          inputId = "heatmap_mode",
          label = "Heatmap Type:",
          choices = c("Static" = "static", "Interactive" = "interactive"),
          selected = "static",
          status = "primary",
          size = "sm"
        )
    ),

    conditionalPanel(
      condition = "input.heatmap_mode == 'interactive'",
      checkboxInput("HeatmapClusterColumns", "Cluster Columns", TRUE),
      checkboxInput("HeatmapClusterRows", "Cluster Rows", TRUE)
    )
  )
}

# For VlnPlot - add ncol option like your example
ui_vlnplot_additions <- function() {
  tagList(
    # Add multi-column layout option
    conditionalPanel(
      condition = "output.Vlnplot_multiple_genes",
      sliderInput("VlnPlotNCols", "Number of Columns:",
                  min = 1, max = 5, value = 3, step = 1)
    )
  )
}

# For Search Features - add import buttons
ui_search_features_additions <- function() {
  tagList(
    # Add at the bottom of the search features box
    conditionalPanel(
      condition = "typeof output.dataset_features !== 'undefined' && output.dataset_features !== null",
      tags$hr(),
      h5("Import Selected Features"),
      fluidRow(
        column(4, actionButton("ImportToViolinPlot", "To Violin Plot",
                              class = "btn-primary btn-sm", width = "100%")),
        column(4, actionButton("ImportToDotPlot", "To Dot Plot",
                              class = "btn-primary btn-sm", width = "100%")),
        column(4, actionButton("ImportToFeaturePlot", "To Feature Plot",
                              class = "btn-primary btn-sm", width = "100%"))
      ),
      textOutput("ImportMessage")
    )
  )
}

# =============================================================================
# 2. MODIFY SERVER.R - Implement the fixes
# =============================================================================

# Replace existing interactive plot outputs with these fixed versions:

server_dimplot_fixed <- function() {
  # Fixed Interactive DimPlot with proper height and tooltips
  output$dimplot_interactive <- plotly::renderPlotly({
    req(input$DimDimensionReduction)
    if(verbose){message("SeuratExplorer: preparing dimplot_interactive_fixed...")}

    cds <- data$obj
    Seurat::Idents(cds) <- input$DimClusterResolution

    # Get highlighted cells (same logic as before)
    if(is.null(input$DimHighlightedClusters)){
      dim_cells_highlighted <- NULL
    } else {
      dim_cells_highlighted <- colnames(cds)[cds@meta.data[,input$DimClusterResolution] %in% input$DimHighlightedClusters]
    }

    cds@meta.data[,input$DimClusterResolution] <- factor(cds@meta.data[,input$DimClusterResolution],
                                                         levels = input$DimClusterOrder)

    # Generate base plot (same as static version)
    if (is.null(DimSplit.Revised())) {
      p <- Seurat::DimPlot(cds,
                           reduction = input$DimDimensionReduction,
                           label = input$DimShowLabel,
                           pt.size = input$DimPointSize,
                           label.size = input$DimLabelSize,
                           group.by = input$DimClusterResolution,
                           cells.highlight = dim_cells_highlighted)
    } else {
      plot_numbers <- length(levels(cds@meta.data[,DimSplit.Revised()]))
      p <- Seurat::DimPlot(cds, reduction = input$DimDimensionReduction,
                           label = input$DimShowLabel, pt.size = input$DimPointSize,
                           label.size = input$DimLabelSize,
                           group.by = input$DimClusterResolution,
                           split.by = DimSplit.Revised(),
                           ncol = ceiling(sqrt(plot_numbers)),
                           cells.highlight = dim_cells_highlighted)
    }

    if(!input$DimShowLegend){
      p <- p & NoLegend()
    }

    # Apply static theme if selected
    if (!is.null(input$DimPlotThemeStatic)) {
      theme_fun <- switch(input$DimPlotThemeStatic,
        "classic" = ggplot2::theme_classic(),
        "minimal" = ggplot2::theme_minimal(),
        "bw" = ggplot2::theme_bw(),
        "void" = ggplot2::theme_void(),
        "dark" = ggplot2::theme_dark(),
        ggplot2::theme_classic()
      )
      p <- p & theme_fun
    }

    # Use fixed interactive function
    interactive_dimplot_fixed(
      p = p,
      obj = cds,
      reduction = input$DimDimensionReduction,
      group.by = input$DimClusterResolution,
      hw_ratio = input$DimPlotHWRatio,  # This will now work!
      width_px = session$clientData$output_dimplot_width %||% 600,
      tooltip_vars = input$DimPlotTooltipVars %||% c("nFeature_RNA", "nCount_RNA"),
      theme_choice = input$DimPlotTheme %||% "plotly_white"
    )
  })
}

server_featureplot_fixed <- function() {
  # Fixed Interactive FeaturePlot
  output$featureplot_interactive <- plotly::renderPlotly({
    req(input$FeatureSlot)
    if(verbose){message("SeuratExplorer: preparing featureplot_interactive_fixed...")}

    # Same validation logic as before
    if(input$FeatureMinCutoff == 0){
      expr_min_cutoff <- NA
    } else {
      expr_min_cutoff <- paste0('q', round(input$FeatureMinCutoff))
    }
    if(input$FeatureMaxCutoff == 100){
      expr_max_cutoff <- NA
    } else {
      expr_max_cutoff <- paste0('q', round(input$FeatureMaxCutoff))
    }

    if (any(is.na(features_dimplot$features_current))) {
      interactive_empty_plot()
    } else {
      cds <- data$obj
      Seurat::Idents(cds) <- input$FeatureClusterResolution
      Seurat::DefaultAssay(cds) <- input$FeatureAssay

      if(!any(features_dimplot$features_current %in% c(rownames(cds[[input$FeatureAssay]]),data$extra_qc_options))){
        interactive_empty_plot()
      } else {
        # Generate base plot (same logic as static)
        if(is.null(FeatureSplit.Revised())) {
          p <- Seurat::FeaturePlot(cds,
                                   features = features_dimplot$features_current,
                                   pt.size = input$FeaturePointSize,
                                   reduction = input$FeatureDimensionReduction,
                                   slot = input$FeatureSlot,
                                   cols = c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor),
                                   label = input$FeatureShowLabel,
                                   label.size = input$FeatureLabelSize,
                                   alpha = input$FeaturePointAlpha,
                                   min.cutoff = expr_min_cutoff,
                                   max.cutoff = expr_max_cutoff)
        } else {
          p <- Seurat::FeaturePlot(cds,
                                   features = features_dimplot$features_current,
                                   pt.size = input$FeaturePointSize,
                                   reduction = input$FeatureDimensionReduction,
                                   slot = input$FeatureSlot,
                                   cols = c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor),
                                   split.by = FeatureSplit.Revised(),
                                   label = input$FeatureShowLabel,
                                   label.size = input$FeatureLabelSize,
                                   alpha = input$FeaturePointAlpha,
                                   min.cutoff = expr_min_cutoff,
                                   max.cutoff = expr_max_cutoff)
          if (length(features_dimplot$features_current) == 1) {
            plot_numbers <- length(levels(cds@meta.data[,FeatureSplit.Revised()]))
            p <- p + patchwork::plot_layout(ncol = ceiling(sqrt(plot_numbers)),
                                           nrow = ceiling(plot_numbers/ceiling(sqrt(plot_numbers))))
          }
        }

        # Apply static theme
        if (!is.null(input$FeaturePlotThemeStatic)) {
          theme_fun <- switch(input$FeaturePlotThemeStatic,
            "classic" = ggplot2::theme_classic(),
            "minimal" = ggplot2::theme_minimal(),
            "bw" = ggplot2::theme_bw(),
            ggplot2::theme_classic()
          )
          p <- p & theme_fun
        }

        # Use fixed interactive function
        interactive_featureplot_fixed(
          p = p,
          obj = cds,
          features = features_dimplot$features_current,
          reduction = input$FeatureDimensionReduction,
          slot = input$FeatureSlot,
          hw_ratio = input$FeaturePlotHWRatio,  # This will now work!
          width_px = session$clientData$output_featureplot_width %||% 600,
          theme_choice = input$FeaturePlotTheme %||% "plotly_white"
        )
      }
    }
  })
}

server_vlnplot_enhanced <- function() {
  # Enhanced VlnPlot with ncol support
  output$vlnplot <- renderPlot({
    if(verbose){message("SeuratExplorer: preparing enhanced vlnplot...")}

    if (any(is.na(features_vlnplot$features_current))) {
      p <- empty_plot
    } else {
      cds <- data$obj
      cds@meta.data[,input$VlnClusterResolution] <- factor(cds@meta.data[,input$VlnClusterResolution],
                                                           levels = input$VlnClusterOrder)
      SeuratObject::Idents(cds) <- input$VlnClusterResolution

      if(!any(features_vlnplot$features_current %in% c(rownames(cds[[input$VlnAssay]]),data$extra_qc_options))){
        p <- empty_plot
      } else {
        # Use ncol parameter for multiple genes
        ncol_value <- if(length(features_vlnplot$features_current) > 1 && !is.null(input$VlnPlotNCols)) {
          input$VlnPlotNCols
        } else {
          NULL
        }

        if(length(features_vlnplot$features_current) == 1) {
          p <- Seurat::VlnPlot(cds,
                               features = features_vlnplot$features_current,
                               assay = input$VlnAssay,
                               layer = input$VlnSlot,
                               split.by = VlnSplit.Revised(),
                               split.plot = input$VlnSplitPlot,
                               pt.size = input$VlnPointSize,
                               alpha = input$VlnPointAlpha,
                               idents = input$VlnIdentsSelected) &
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                           axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize))
        } else {
          p <- Seurat::VlnPlot(cds,
                               features = features_vlnplot$features_current,
                               assay = input$VlnAssay,
                               layer = input$VlnSlot,
                               split.by = VlnSplit.Revised(),
                               split.plot = input$VlnSplitPlot,
                               stack = input$VlnStackPlot,
                               flip = input$VlnFlipPlot,
                               fill.by = input$VlnFillBy,
                               idents = input$VlnIdentsSelected,
                               pt.size = input$VlnPointSize,
                               alpha = input$VlnPointAlpha,
                               ncol = ncol_value) &  # Add ncol parameter!
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                           axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize))
        }

        if (input$Vlnfillcolorplatte != 'default' & input$VlnSplitBy == 'None'){
          fill.colors <- getColors(color.platte = color_list,
                                   choice = input$Vlnfillcolorplatte,
                                   n = length(levels(Idents(cds))))
          names(fill.colors) <- levels(Idents(cds))
          p <- p & scale_fill_manual(values = fill.colors)
        }
      }
    }

    # Save plot
    ggplot2::ggsave(paste0(temp_dir,"/vlnplot.pdf"),
                    p,
                    width = vlnplot_width() * px2cm,
                    height = vlnplot_width() * input$VlnPlotHWRatio * px2cm,
                    units = "cm",
                    limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_vlnplot_width * input$VlnPlotHWRatio})
}

# Feature import functionality
server_feature_import <- function() {
  # Get selected features from DT
  selected_features <- reactive({
    selected_rows <- input$dataset_features_rows_selected
    if (is.null(selected_rows) || length(selected_rows) == 0) {
      return(NULL)
    }

    # Get the features from selected rows
    features_df <- data$features_df  # Assuming this exists
    selected_genes <- features_df$FeatureName[selected_rows]
    return(selected_genes)
  })

  # Import to Violin Plot
  observeEvent(input$ImportToViolinPlot, {
    features <- selected_features()
    if (!is.null(features)) {
      updateTextAreaInput(session, "VlnGeneSymbol",
                         value = paste(features, collapse = "\n"))
      # Switch to violin plot tab
      updateTabItems(session, "tabs", selected = "vlnplot")

      output$ImportMessage <- renderText({
        paste("Imported", length(features), "features to Violin Plot")
      })
    }
  })

  # Import to Dot Plot
  observeEvent(input$ImportToDotPlot, {
    features <- selected_features()
    if (!is.null(features)) {
      updateTextAreaInput(session, "DotGeneSymbol",
                         value = paste(features, collapse = "\n"))
      updateTabItems(session, "tabs", selected = "dotplot")

      output$ImportMessage <- renderText({
        paste("Imported", length(features), "features to Dot Plot")
      })
    }
  })

  # Import to Feature Plot
  observeEvent(input$ImportToFeaturePlot, {
    features <- selected_features()
    if (!is.null(features)) {
      updateTextAreaInput(session, "FeatureGeneSymbol",
                         value = paste(features, collapse = "\n"))
      updateTabItems(session, "tabs", selected = "featureplot")

      output$ImportMessage <- renderText({
        paste("Imported", length(features), "features to Feature Plot")
      })
    }
  })
}
