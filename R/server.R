# server.R
## the server side

#' server side functions related to `explorer_sidebar_ui`
#'
#' @param input server input
#' @param output server output
#' @param session server session
#' @param verbose for debug use
#' @param data the Seurat object and related parameters
#'
#' @import Seurat SeuratObject
#' @importFrom grDevices dev.off pdf
#' @export
#' @return server side functions related to `explorer_sidebar_ui`
#'
explorer_server <- function(input, output, session, data, verbose=FALSE){
  temp_dir <- tempdir() # temporary directory, for save plots

  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }

  dir.create(temp_dir, showWarnings = FALSE)
  # to make shinyBS::updateCollapse() runs correctly, refer to: https://github.com/ebailey78/shinyBS/issues/92
  shiny::addResourcePath("sbs", system.file("www", package="shinyBS"))

  # Using an un-exported function from annother R package:
  # https://stackoverflow.com/questions/32535773/using-un-exported-function-from-another-r-package
  subset_Seurat <- utils::getFromNamespace('subset.Seurat', 'SeuratObject')

  ############################# Dimension Reduction Plot
  # define reductions choices UI
  output$DimReductions.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DimReductions.UI...")}
    selectInput("DimDimensionReduction", "Dimension Reduction:", choices = data$reduction_options, selected = data$reduction_default) # set default reduction
  })

  # define Cluster Annotation choice
  output$DimClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DimClusterResolution.UI...")}
    selectInput("DimClusterResolution","Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Cluster order
  output$DimClusterOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DimClusterOrder.UI...")}
    shinyjqui::orderInput(inputId = 'DimClusterOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$DimClusterResolution]),width = '100%')
  })

  # when change cluster resolution, open the shinyBS::bsCollapsePanel, otherwise will cause cluster order not update
  # a bad effect is: each time changing the resolution option, will collapse cluster order ui
  observeEvent(input$DimClusterResolution, ({
    if(verbose){message("SeuratExplorer: updateCollapse for collapseDimplot...")}
    shinyBS::updateCollapse(session, "collapseDimplot", open = "Change Cluster Order")
  }))

  # define Split Choice UI
  output$DimSplit.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DimSplit.UI...")}
    selectInput("DimSplit","Split by:", choices = c("None" = "None", data$split_options))
  })

  # Revise Split selection which will be appropriate for plot
  DimSplit.Revised <- reactive({
    req(input$DimSplit) # only run after split is ready
    if(verbose){message("SeuratExplorer: preparing DimSplit.Revised...")}
    # Revise the Split choice
    if(is.na(input$DimSplit) | input$DimSplit == "None") {
      return(NULL)
    }else{
      return(input$DimSplit)
    }
  })

  # define Cluster choice for highlight
  output$DimHighlightedClusters.UI <- renderUI({
    req(input$DimClusterResolution)
    if(verbose){message("SeuratExplorer: preparing DimHighlightedClusters.UI...")}
    shinyWidgets::pickerInput(inputId = "DimHighlightedClusters", label = "Highlight Clusters:",
                              choices = levels(data$obj@meta.data[,input$DimClusterResolution]),
                              selected = NULL,
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  dimplot_width  <- reactive({ session$clientData$output_dimplot_width })

  # Pixel (X) to Centimeter: 1 pixel (X)	= 0.0264583333 cm, if use this value, the picture is a little bit of small, unknown why.
  px2cm <- 0.03

  output$dimplot <- renderPlot({
    if(verbose){message("SeuratExplorer: preparing dimplot...")}
    isolate(cds <- data$obj) # not a memory saving way
    # for highlight cells
    if (any(is.null(input$DimHighlightedClusters))) {
      dim_cells_highlighted <- NULL
    }else{
      dim_cells_highlighted <- colnames(cds)[cds@meta.data[,input$DimClusterResolution] %in% input$DimHighlightedClusters]
    }

    cds@meta.data[,input$DimClusterResolution] <- factor(cds@meta.data[,input$DimClusterResolution], levels = input$DimClusterOrder)
    if (is.null(DimSplit.Revised())) { # not splited
      p <- Seurat::DimPlot(cds, reduction = input$DimDimensionReduction, label = input$DimShowLabel, pt.size = input$DimPointSize, label.size = input$DimLabelSize,
                      group.by = input$DimClusterResolution, cells.highlight = dim_cells_highlighted)
      }else{ # splited
      plot_numbers <- length(levels(cds@meta.data[,DimSplit.Revised()]))
      p <- Seurat::DimPlot(cds, reduction = input$DimDimensionReduction, label = input$DimShowLabel, pt.size = input$DimPointSize, label.size = input$DimLabelSize,
                      group.by = input$DimClusterResolution, split.by = DimSplit.Revised(), ncol = ceiling(sqrt(plot_numbers)), cells.highlight = dim_cells_highlighted)
      }
    ggplot2::ggsave(paste0(temp_dir,"/dimplot.pdf"), p, width = dimplot_width() * px2cm, height = dimplot_width() * input$DimPlotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_dimplot_width * input$DimPlotHWRatio}) # box plot: height = width default

  # refer to: https://stackoverflow.com/questions/14810409/how-to-save-plots-that-are-made-in-a-shiny-app
  output$downloaddimplot <- downloadHandler(
    filename = function(){'dimplot.pdf'},
    content = function(file) {
      file.copy(paste0(temp_dir,"/dimplot.pdf"), file, overwrite=TRUE)
    })

  ################################ Feature Plot
  # define reductions choices UI
  output$FeatureReductions.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing FeatureReductions.UI...")}
    selectInput("FeatureDimensionReduction", "Dimension Reduction:", choices = data$reduction_options, selected = data$reduction_default) # set default reduction
  })

  # define Cluster Annotation choice
  output$FeatureClusterResolution.UI <- renderUI({
    selectInput("FeatureClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  })

  # define Split Choice UI
  output$FeatureSplit.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing FeatureSplit.UI...")}
    selectInput("FeatureSplit","Split by:", choices = c("None" = "None", data$split_options))
  })

  # inform extra qc options for Gene symbol input
  output$Featurehints.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing Featurehints.UI...")}
    helpText(strong(paste("Also supports: ", paste(data$extra_qc_options, collapse = " "), ".",sep = "")),
             br(),
             strong("Tips: You can paste multiple genes from a column in excel."),style = "font-size:12px;")
  })


  # Revise Split selection which will be appropriate for DimPlot, FeaturePlot and Vlnplot functions.
  FeatureSplit.Revised <- reactive({
    req(input$FeatureSplit)
    if(verbose){message("SeuratExplorer: preparing FeatureSplit.Revised...")}
    # Revise the Split choice
    if(is.na(input$FeatureSplit) | input$FeatureSplit == "None") {
      return(NULL)
    }else{
      return(input$FeatureSplit)
    }
  })

  # Check the input gene
  Featureplot.Gene.Revised <- reactive({
    req(input$FeatureGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Featureplot.Gene.Revised...")}
    ifelse(is.na(input$FeatureGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$FeatureGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  featureplot_width  <- reactive({ session$clientData$output_featureplot_width })

  output$featureplot <- renderPlot({
    if(verbose){message("SeuratExplorer: preparing featureplot...")}
    if(input$FeatureMinCutoff == 0){expr_min_cutoff <- NA}else{expr_min_cutoff <- paste0('q', round(input$FeatureMinCutoff))}
    if(input$FeatureMaxCutoff == 100){expr_max_cutoff <- NA}else{expr_max_cutoff <- paste0('q', round(input$FeatureMaxCutoff))}
    if (any(is.na(Featureplot.Gene.Revised()))) { # when NA value
      p <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when all wrong input, show a blank pic.
    }else if (input$FeatureShowLabel) { # show label
      isolate(cds <- data$obj)
      Idents(cds) <- input$FeatureClusterResolution
      if(is.null(FeatureSplit.Revised())) { # not splited
        p <- Seurat::FeaturePlot(cds, features = Featureplot.Gene.Revised(), pt.size = input$FeaturePointSize, reduction = input$FeatureDimensionReduction,
                                 cols = c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor), label = TRUE, label.size = input$FeatureLabelSize,
                                 alpha = input$FeaturePointAlpha, min.cutoff = expr_min_cutoff, max.cutoff = expr_max_cutoff)
      }else{ # splited
        p <- Seurat::FeaturePlot(cds, features = Featureplot.Gene.Revised(), pt.size = input$FeaturePointSize, reduction = input$FeatureDimensionReduction,
                                 cols =  c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor), split.by = FeatureSplit.Revised(), label = TRUE,
                                 label.size = input$FeatureLabelSize, alpha = input$FeaturePointAlpha, min.cutoff = expr_min_cutoff, max.cutoff = expr_max_cutoff)
        if (length( Featureplot.Gene.Revised()) == 1) { # only one gene
          plot_numbers <- length(levels(cds@meta.data[,FeatureSplit.Revised()]))
          p <- p + patchwork::plot_layout(ncol = ceiling(sqrt(plot_numbers)),nrow = ceiling(plot_numbers/ceiling(sqrt(plot_numbers))))
        }
      }
    }else{ # not show label
      if(is.null(FeatureSplit.Revised())) { # not splited
        p <- Seurat::FeaturePlot(data$obj, features = Featureplot.Gene.Revised(), pt.size = input$FeaturePointSize, reduction = input$FeatureDimensionReduction,
                                 cols = c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor), alpha = input$FeaturePointAlpha,
                                 min.cutoff = expr_min_cutoff, max.cutoff = expr_max_cutoff)
      }else{ # splited
        p <- Seurat::FeaturePlot(data$obj, features = Featureplot.Gene.Revised(), pt.size = input$FeaturePointSize, reduction = input$FeatureDimensionReduction,
                                 cols =  c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor), split.by = FeatureSplit.Revised(),
                                 alpha = input$FeaturePointAlpha, min.cutoff = expr_min_cutoff, max.cutoff = expr_max_cutoff)
        if (length( Featureplot.Gene.Revised()) == 1) { # only one gene
          plot_numbers <- length(levels(data$obj@meta.data[,FeatureSplit.Revised()]))
          p <- p + patchwork::plot_layout(ncol = ceiling(sqrt(plot_numbers)),nrow = ceiling(plot_numbers/ceiling(sqrt(plot_numbers))))
        }
      }
    }
    ggplot2::ggsave(paste0(temp_dir,"/featureplot.pdf"), p, width = featureplot_width() * px2cm, height = featureplot_width() * input$FeaturePlotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_featureplot_width * input$FeaturePlotHWRatio}) # box plot: height = width default


  output$downloadfeatureplot <- downloadHandler(
    filename = function(){'featureplot.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/featureplot.pdf"))) { # problem: will throw an error when file not exists; or with a uncorrected input, will download the pic of previous corrected input.
        file.copy(paste0(temp_dir,"/featureplot.pdf"), file, overwrite=TRUE)
      }
    })

  ################################ Violin Plot
  # Check the input gene
  Vlnplot.Gene.Revised <- reactive({
    req(input$VlnGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Vlnplot.Gene.Revised...")}
    ifelse(is.na(input$VlnGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$VlnGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  output$Vlnhints.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing Vlnhints.UI...")}
    helpText(strong(paste("Also supports: ", paste(data$extra_qc_options, collapse = " "), ".",sep = "")),
             br(),
             strong("Tips: You can paste multiple genes from a column in excel."),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$VlnClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing VlnClusterResolution.UI...")}
    selectInput("VlnClusterResolution","Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Cluster order
  output$VlnClusterOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing VlnClusterOrder.UI...")}
    shinyjqui::orderInput(inputId = 'VlnClusterOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$VlnClusterResolution]),width = '100%')
  })

  # when change cluster resolution, open the shinyBS::bsCollapsePanel, otherwise will cause cluster order not update
  observeEvent(input$VlnClusterResolution, ({
    if(verbose){message("SeuratExplorer: updateCollapse for collapseVlnplot...")}
    shinyBS::updateCollapse(session, "collapseVlnplot", open = "0")
  }))

  # define the idents used
  output$VlnIdentsSelected.UI <- renderUI({
    req(input$VlnClusterResolution)
    if(verbose){message("SeuratExplorer: preparing VlnIdentsSelected.UI...")}
    shinyWidgets::pickerInput(inputId = "VlnIdentsSelected", label = "Clusters Used:",
                              choices = levels(data$obj@meta.data[,input$VlnClusterResolution]), selected = levels(data$obj@meta.data[,input$VlnClusterResolution]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  # define Split Choice UI
  output$VlnSplitBy.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing VlnSplitBy.UI...")}
    selectInput("VlnSplitBy","Split by:", choices = c("None" = "None", data$split_options))
  })

  # Conditional panel: show this panel when split.by is selected and the the level equals to 2
  output$Vlnplot_splitoption_twolevels = reactive({
    req(input$VlnSplitBy)
    if(verbose){message("SeuratExplorer: preparing Vlnplot_splitoption_twolevels...")}
    if (input$VlnSplitBy == "None"){
      return(FALSE)
    }else if(length(levels(data$obj@meta.data[,input$VlnSplitBy])) == 2) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  # Disable suspend for output$file_loaded,
  # When TRUE (the default), the output object will be suspended (not execute) when it is hidden on the web page.
  # When FALSE, the output object will not suspend when hidden, and if it was already hidden and suspended, then it will resume immediately.
  outputOptions(output, 'Vlnplot_splitoption_twolevels', suspendWhenHidden = FALSE)

  # Conditional panel: show this panel when input multiple gene symbols
  output$Vlnplot_multiple_genes = reactive({
    req(input$VlnGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Vlnplot_multiple_genes...")}
    if (length(Vlnplot.Gene.Revised()) > 1) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  outputOptions(output, 'Vlnplot_multiple_genes', suspendWhenHidden = FALSE)


  # Conditional panel: show this panel when input multiple genes and stack is set to TRUE
  output$Vlnplot_StackPlot = reactive({
    req(input$VlnStackPlot)
    req(input$VlnGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Vlnplot_StackPlot...")}
    if (length(Vlnplot.Gene.Revised()) > 1 & input$VlnStackPlot) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  outputOptions(output, 'Vlnplot_StackPlot', suspendWhenHidden = FALSE)

  # Revise Split selection which will be appropriate for DimPlot, FeaturePlot and Vlnplot functions.
  VlnSplit.Revised <- reactive({
    if(verbose){message("SeuratExplorer: preparing VlnSplit.Revised...")}
    req(input$VlnSplitBy)
    # Revise the Split choice
    if(is.na(input$VlnSplitBy) | input$VlnSplitBy == "None") {
      return(NULL)
    }else{
      return(input$VlnSplitBy)
    }
  })

  # reset VlnSplitPlot value to FALSE when change the split options
  observe({
    req(input$VlnSplitBy)
    if(verbose){message("SeuratExplorer: vlnplot update UI...")}
    updateCheckboxInput(session, "VlnSplitPlot", value = FALSE)
    updateCheckboxInput(session, "VlnStackPlot", value = FALSE)
    updateCheckboxInput(session, "VlnFlipPlot", value = FALSE)
    updateSelectInput(session, "VlnFillBy", selected = "feature")
  })

  # shiny related bug
  # debug in future! 2024.05.15
  # how to make sure renderPlot run after the observe(input$VlnSplitBy)[Warning: Error in SingleExIPlot: Unknown plot type: splitViolin,
  # for the VlnSplitPlot is not updated

  # seurat related bug
  # VlnPlot(cds,features = c("CD4","CD8A"),split.by = "orig.ident", stack = TRUE,group.by = "cca_clusters_res_0.2",flip = FALSE,split.plot = TRUE)
  # Error:
  # Error in `vln.geom()`:
  #   ! Problem while converting geom to grob.
  # Caused by error in `$<-.data.frame`:
  # Run `rlang::last_trace()` to see where the error occurred
  # not related to ggplot2, pathcwork, rlang versions

  vlnplot_width  <- reactive({ session$clientData$output_vlnplot_width })

  output$vlnplot <- renderPlot({
    req(Vlnplot.Gene.Revised())
    if(verbose){message("SeuratExplorer: preparing vlnplot...")}
    if (any(is.na(Vlnplot.Gene.Revised()))) { # when NA value
      p <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj)
      cds@meta.data[,input$VlnClusterResolution] <- factor(cds@meta.data[,input$VlnClusterResolution], levels = input$VlnClusterOrder)
      SeuratObject::Idents(cds) <- input$VlnClusterResolution
      if(length(Vlnplot.Gene.Revised()) == 1) { # only One Gene
        p <- Seurat::VlnPlot(cds, features = Vlnplot.Gene.Revised(), split.by = VlnSplit.Revised(), split.plot = input$VlnSplitPlot,
                             pt.size = input$VlnPointSize, alpha = input$VlnPointAlpha, idents = input$VlnIdentsSelected) &
          ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                         axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize))
      }else{ # multiple genes
        p <- Seurat::VlnPlot(cds, features = Vlnplot.Gene.Revised(), split.by = VlnSplit.Revised(), split.plot = input$VlnSplitPlot,
                             stack = input$VlnStackPlot,flip = input$VlnFlipPlot, fill.by = input$VlnFillBy, idents = input$VlnIdentsSelected,
                             pt.size = input$VlnPointSize, alpha = input$VlnPointAlpha) &
          ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                         axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize))
      }
      ggplot2::ggsave(paste0(temp_dir,"/vlnplot.pdf"), p, width = vlnplot_width() * px2cm, height = vlnplot_width() * input$VlnPlotHWRatio * px2cm, units = "cm", limitsize = FALSE)
      return(p)
    }
  }, height = function(){session$clientData$output_vlnplot_width * input$VlnPlotHWRatio}) # box plot: height = width default


  output$downloadvlnplot <- downloadHandler(
    filename = function(){'vlnplot.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/vlnplot.pdf"))) {
        file.copy(paste0(temp_dir,"/vlnplot.pdf"), file, overwrite=TRUE)
      }
    })


  ################################ Dot Plot

  # Check the input gene
  Dotplot.Gene.Revised <- reactive({
    req(input$DotGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Dotplot.Gene.Revised...")}
    ifelse(is.na(input$DotGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$DotGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  output$Dothints.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing Dothints.UI...")}
    helpText(strong(paste("Also supports: ", paste(data$extra_qc_options, collapse = " "), ".",sep = "")),
             br(),
             strong("Tips: You can paste multiple genes from a column in excel."),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$DotClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DotClusterResolution.UI...")}
    selectInput("DotClusterResolution","Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Cluster order
  output$DotClusterOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DotClusterOrder.UI...")}
    shinyjqui::orderInput(inputId = 'DotClusterOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$DotClusterResolution]),width = '100%')
  })

  # when change cluster resolution, open the shinyBS::bsCollapsePanel, otherwise will cause cluster order not update
  observeEvent(input$DotClusterResolution, ({
    if(verbose){message("SeuratExplorer: updateCollapse for collapseDotplot...")}
    shinyBS::updateCollapse(session, "collapseDotplot", open = "0")
  }))

  # define the idents used
  output$DotIdentsSelected.UI <- renderUI({
    req(input$DotClusterResolution)
    if(verbose){message("SeuratExplorer: preparing DotIdentsSelected.UI...")}
    shinyWidgets::pickerInput(inputId = "DotIdentsSelected", label = "Clusters Used:",
                              choices = levels(data$obj@meta.data[,input$DotClusterResolution]), selected = levels(data$obj@meta.data[,input$DotClusterResolution]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  # define Split Choice UI
  output$DotSplitBy.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing DotSplitBy.UI...")}
    selectInput("DotSplitBy","Split by:", choices = c("None" = "None", data$split_options))
  })

  # Revise Split selection which will be appropriate for DimPlot, FeaturePlot and Vlnplot functions.
  DotSplit.Revised <- reactive({
    req(input$DotSplitBy)
    if(verbose){message("SeuratExplorer: preparing DotSplit.Revised...")}
    # Revise the Split choice
    if(is.na(input$DotSplitBy) | input$DotSplitBy == "None") {
      return(NULL)
    }else{
      return(input$DotSplitBy)
    }
  })

  # Conditional panel: when split is NULL, You can set the corresponding color for highest and lowest value,
  # when split is not NULL, ggplot2 will generate colors for point.
  output$DotPlot_Split_isNone <- reactive({
    req(input$DotSplitBy)
    if(verbose){message("SeuratExplorer: preparing DotPlot_Split_isNone...")}
    if(is.na(input$DotSplitBy) | input$DotSplitBy == "None") {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  outputOptions(output, 'DotPlot_Split_isNone', suspendWhenHidden = FALSE)

  dotplot_width  <- reactive({ session$clientData$output_dotplot_width })

  output$dotplot <- renderPlot({
    req(Dotplot.Gene.Revised())
    if(verbose){message("SeuratExplorer: preparing dotplot...")}
    if (any(is.na(Dotplot.Gene.Revised()))) { # NA
      p <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj)
      Idents(cds) <- input$DotClusterResolution
      cds@meta.data[,input$DotClusterResolution] <- factor(cds@meta.data[,input$DotClusterResolution], levels = input$DotClusterOrder)
      if (is.null(DotSplit.Revised())) {
        p <- Seurat::DotPlot(cds, features = Dotplot.Gene.Revised(), group.by = input$DotClusterResolution,
                             idents = input$DotIdentsSelected,
                             split.by = DotSplit.Revised(), cluster.idents = input$DotClusterIdents, dot.scale = input$DotDotScale,
                             cols = c(input$DotPlotLowestExprColor, input$DotPlotHighestExprColor))
      }else{
        split.levels.length <- length(levels(cds@meta.data[,DotSplit.Revised()]))
        p <- Seurat::DotPlot(cds, features = Dotplot.Gene.Revised(), group.by = input$DotClusterResolution,
                             idents = input$DotIdentsSelected,
                             split.by = DotSplit.Revised(), cluster.idents = input$DotClusterIdents, dot.scale = input$DotDotScale,
                             cols = scales::hue_pal()(split.levels.length))
      }
      p <- p & ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$DotXlabelSize), axis.text.y = ggplot2::element_text(size = input$DotYlabelSize))
      if (input$DotRotateAxis) { p <- p + Seurat::RotatedAxis() }
      if (input$DotFlipCoordinate) { p <- p + ggplot2::coord_flip() }
    }
    ggplot2::ggsave(paste0(temp_dir,"/dotplot.pdf"), p, width = dotplot_width() * px2cm, height = dotplot_width() * input$DotPlotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_dotplot_width * input$DotPlotHWRatio}) # box plot: height = width default


  output$downloaddotplot <- downloadHandler(
    filename = function(){'dotplot.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/dotplot.pdf"))) {
        file.copy(paste0(temp_dir,"/dotplot.pdf"), file, overwrite=TRUE)
      }
    })

  ################################ Heatmap

  # Check the input gene
  Heatmap.Gene.Revised <- reactive({
    req(input$HeatmapGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Heatmap.Gene.Revised...")}
    ifelse(is.na(input$HeatmapGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$HeatmapGeneSymbol, GeneLibrary =  rownames(data$obj))))
  })

  output$Heatmaphints.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing Heatmaphints.UI...")}
    helpText(strong("Tips: You can paste multiple genes from a column in excel."),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$HeatmapClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing HeatmapClusterResolution.UI...")}
    selectInput("HeatmapClusterResolution","Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Cluster order
  output$HeatmapClusterOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing HeatmapClusterOrder.UI...")}
    shinyjqui::orderInput(inputId = 'HeatmapClusterOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$HeatmapClusterResolution]),width = '100%')
  })

  observeEvent(input$HeatmapClusterResolution, ({
    if(verbose){message("SeuratExplorer: updateCollapse for collapseHeatmap...")}
    shinyBS::updateCollapse(session, "collapseHeatmap", open = "0")
  }))

  heatmap_width  <- reactive({ session$clientData$output_heatmap_width })

  output$heatmap <- renderPlot({
    req(Heatmap.Gene.Revised())
    if(verbose){message("SeuratExplorer: preparing heatmap...")}
    if (any(is.na(Heatmap.Gene.Revised()))) { # NA
      p <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj)
      cds@meta.data[,input$HeatmapClusterResolution] <- factor(cds@meta.data[,input$HeatmapClusterResolution], levels = input$HeatmapClusterOrder)
      if (!all(Heatmap.Gene.Revised() %in% Seurat::VariableFeatures(cds))) {
        cds <- Seurat::ScaleData(object = cds, features = unique(c(Seurat::VariableFeatures(cds), Heatmap.Gene.Revised()))) # use only one gene to scaledata() will throw an error
      }

      p <- Seurat::DoHeatmap(object = cds, features = Heatmap.Gene.Revised(), group.by = input$HeatmapClusterResolution, size = input$HeatmapTextSize,
                        hjust = input$HeatmapTextHjust, vjust = input$HeatmapTextVjust, angle = input$HeatmapTextRatateAngle,
                        group.bar.height = input$HeatmapGroupBarHeight, lines.width = input$HeatmapLineWidth) &
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = input$HeatmapFeatureTextSize))
    }
    ggplot2::ggsave(paste0(temp_dir,"/heatmap.pdf"), p, width = heatmap_width() * px2cm, height = heatmap_width() * input$HeatmapPlotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_heatmap_width * input$HeatmapPlotHWRatio}) # box plot: height = width default


  output$downloadheatmap <- downloadHandler(
    filename = function(){'heatmap.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/heatmap.pdf"))) {
        file.copy(paste0(temp_dir,"/heatmap.pdf"), file, overwrite=TRUE)
      }
    })

  ################################ Averaged Heatmap
  # Check the input gene
  AveragedHeatmap.Gene.Revised <- reactive({
    req(input$AveragedHeatmapGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing AveragedHeatmap.Gene.Revised...")}
    ifelse(is.na(input$AveragedHeatmapGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$AveragedHeatmapGeneSymbol, GeneLibrary =  rownames(data$obj))))
  })

  output$AveragedHeatmaphints.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing AveragedHeatmaphints.UI...")}
    helpText(strong("Tips: You can paste multiple genes from a column in excel."),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$AveragedHeatmapClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing AveragedHeatmapClusterResolution.UI...")}
    selectInput("AveragedHeatmapClusterResolution","Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Cluster order
  output$AveragedHeatmapClusterOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing AveragedHeatmapClusterOrder.UI...")}
    shinyjqui::orderInput(inputId = 'AveragedHeatmapClusterOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$AveragedHeatmapClusterResolution]),width = '100%')
  })

  observeEvent(input$AveragedHeatmapClusterResolution, ({
    if(verbose){message("SeuratExplorer: updateCollapse for AveragedcollapseHeatmap...")}
    shinyBS::updateCollapse(session, "AveragedcollapseHeatmap", open = "0")
  }))

  # define the idents used
  output$AveragedHeatmapIdentsSelected.UI <- renderUI({
    req(input$AveragedHeatmapClusterResolution)
    if(verbose){message("SeuratExplorer: preparing AveragedHeatmapIdentsSelected.UI...")}
    shinyWidgets::pickerInput(inputId = "AveragedHeatmapIdentsSelected", label = "Clusters Used:",
                              choices = levels(data$obj@meta.data[,input$AveragedHeatmapClusterResolution]), selected = levels(data$obj@meta.data[,input$AveragedHeatmapClusterResolution]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  averagedheatmap_width  <- reactive({ session$clientData$output_averagedheatmap_width })

  output$averagedheatmap <- renderPlot({
    req(AveragedHeatmap.Gene.Revised())
    if(verbose){message("SeuratExplorer: preparing averagedheatmap...")}
    if (any(is.na(AveragedHeatmap.Gene.Revised()))) { # NA
      p <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj)
      cds@meta.data[,input$AveragedHeatmapClusterResolution] <- factor(cds@meta.data[,input$AveragedHeatmapClusterResolution], levels = input$AveragedHeatmapClusterOrder)
      Seurat::Idents(cds) <- input$AveragedHeatmapClusterResolution
      cds <- subset_Seurat(cds, idents = input$AveragedHeatmapIdentsSelected)
      p <- suppressMessages(AverageHeatmap(object = cds, markerGene = AveragedHeatmap.Gene.Revised(), group.by = input$AveragedHeatmapClusterResolution,
                          feature.fontsize = input$AveragedHeatmapFeatureTextSize,cluster.fontsize = input$AveragedHeatmapClusterTextSize,
                          column_names_rot = input$AveragedHeatmapClusterTextRatateAngle, cluster_columns = input$AveragedHeatmapClusterClusters,
                          cluster_rows = input$AveragedHeatmapClusterFeatures))
    }
    pdf(file = paste0(temp_dir,"/AveragedHeatmap.pdf"), width = (averagedheatmap_width() * px2cm)/2.54, height = (averagedheatmap_width() * input$AveragedHeatmapPlotHWRatio * px2cm)/2.54)
    print(p)
    dev.off()
    # ggplot2::ggsave(paste0(temp_dir,"/AveragedHeatmap.pdf"), p, width = averagedheatmap_width() * px2cm, height = averagedheatmap_width() * input$AveragedHeatmapPlotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_averagedheatmap_width * input$AveragedHeatmapPlotHWRatio}) # box plot: height = width default


  output$downloadaveragedheatmap <- downloadHandler(
    filename = function(){'AveragedHeatmap.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/AveragedHeatmap.pdf"))) {
        file.copy(paste0(temp_dir,"/AveragedHeatmap.pdf"), file, overwrite=TRUE)
      }
    })

  ################################ Ridge Plot

  # Check the input gene
  Ridgeplot.Gene.Revised <- reactive({
    req(input$RidgeplotGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Ridgeplot.Gene.Revised...")}
    ifelse(is.na(input$RidgeplotGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$RidgeplotGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  output$Ridgeplothints.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing Ridgeplothints.UI...")}
    helpText(strong(paste("Also supports: ", paste(data$extra_qc_options, collapse = " "), ".",sep = "")),
             br(),
             strong("Tips: You can paste multiple genes from a column in excel."),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$RidgeplotClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing RidgeplotClusterResolution.UI...")}
    selectInput("RidgeplotClusterResolution","Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Cluster order
  output$RidgeplotClusterOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing RidgeplotClusterOrder.UI...")}
    shinyjqui::orderInput(inputId = 'RidgeplotClusterOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$RidgeplotClusterResolution]),width = '100%')
  })

  observeEvent(input$RidgeplotClusterResolution, ({
    if(verbose){message("SeuratExplorer: updateCollapse for collapseRidgeplot...")}
    shinyBS::updateCollapse(session, "collapseRidgeplot", open = "0")
  }))

  # define the idents used
  output$RidgeplotIdentsSelected.UI <- renderUI({
    req(input$RidgeplotClusterResolution)
    if(verbose){message("SeuratExplorer: preparing RidgeplotIdentsSelected.UI...")}
    shinyWidgets::pickerInput(inputId = "RidgeplotIdentsSelected", label = "Clusters Used:",
                              choices = levels(data$obj@meta.data[,input$RidgeplotClusterResolution]), selected = levels(data$obj@meta.data[,input$RidgeplotClusterResolution]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  # Conditional panel: show this panel when input multiple genes and stack is set to TRUE
  output$Ridgeplot_stack_show = reactive({
    req(input$RidgeplotGeneSymbol)
    if(verbose){message("SeuratExplorer: preparing Ridgeplot_stack_show...")}
    if (length(Ridgeplot.Gene.Revised()) > 1) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  outputOptions(output, 'Ridgeplot_stack_show', suspendWhenHidden = FALSE)

  # Conditional panel: show this panel when input multiple genes and stack is set to TRUE
  output$Ridgeplot_stack_NotSelected = reactive({
    req(input$RidgeplotStackPlot)
    if(verbose){message("SeuratExplorer: preparing Ridgeplot_stack_NotSelected...")}
    !input$RidgeplotStackPlot
  })

  outputOptions(output, 'Ridgeplot_stack_NotSelected', suspendWhenHidden = FALSE)

  # reset VlnSplitPlot value to FALSE when change the input gene symbols
  observe({
    req(input$RidgeplotGeneSymbol)
    if(verbose){message("SeuratExplorer: update RidgeplotStackPlot...")}
    updateCheckboxInput(session, "RidgeplotStackPlot", value = FALSE)
  })

  ridgeplot_width  <- reactive({ session$clientData$output_ridgeplot_width })

  output$ridgeplot <- renderPlot({
    req(Ridgeplot.Gene.Revised())
    if(verbose){message("SeuratExplorer: preparing ridgeplot...")}
    if (any(is.na(Ridgeplot.Gene.Revised()))) { # NA
      p <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj)
      cds@meta.data[,input$RidgeplotClusterResolution] <- factor(cds@meta.data[,input$RidgeplotClusterResolution], levels = input$RidgeplotClusterOrder)
      Idents(cds) <- input$RidgeplotClusterResolution
      p <- Seurat::RidgePlot(object = cds, features = Ridgeplot.Gene.Revised(), ncol = input$RidgeplotNumberOfColumns,stack = input$RidgeplotStackPlot,
                             fill.by = input$RidgeplotFillBy, idents = input$RidgeplotIdentsSelected) &
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$RidgeplotXlabelSize),
                       axis.text.y = ggplot2::element_text(size = input$RidgeplotYlabelSize))
    }
    ggplot2::ggsave(paste0(temp_dir,"/ridgeplot.pdf"), p, width = ridgeplot_width() * px2cm, height = ridgeplot_width() * input$RidgeplotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_ridgeplot_width * input$RidgeplotHWRatio}) # box plot: height = width default


  output$downloadridgeplot <- downloadHandler(
    filename = function(){'ridgeplot.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/ridgeplot.pdf"))) {
        file.copy(paste0(temp_dir,"/ridgeplot.pdf"), file, overwrite=TRUE)
      }
    })

  ################################ Cell ratio Plot
  # define Fill choices
  output$CellratioFillChoice.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing CellratioFillChoice.UI...")}
    selectInput("CellratioFillChoice","Fill in choice:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define Fill order
  output$CellratioplotFillOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing CellratioplotFillOrder.UI...")}
    shinyjqui::orderInput(inputId = 'CellratioFillOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$CellratioFillChoice]),width = '100%')
  })


  # define X choices
  output$CellratioXChoice.UI <- renderUI({
    req(input$CellratioFillChoice)
    if(verbose){message("SeuratExplorer: preparing CellratioXChoice.UI...")}
    selectInput("CellratioXChoice","X axis choice:", choices = data$cluster_options[!data$cluster_options %in% input$CellratioFillChoice])
  })


  # define x choice order
  output$CellratioplotXOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing CellratioplotXOrder.UI...")}
    shinyjqui::orderInput(inputId = 'CellratioXOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$CellratioXChoice]),width = '100%')
  })

  # define Facet choices
  output$CellratioFacetChoice.UI <- renderUI({
    req(input$CellratioXChoice)
    if(verbose){message("SeuratExplorer: preparing CellratioFacetChoice.UI...")}
    selectInput("CellratioFacetChoice","Facet choice:", choices = c("None" = "None", data$cluster_options[!data$cluster_options %in% c(input$CellratioFillChoice, input$CellratioXChoice)]), selected = "None")
  })

  # Revise FacetChoice which will be appropriate for plot
  FacetChoice.Revised <- reactive({
    req(input$CellratioFacetChoice)
    if(verbose){message("SeuratExplorer: FacetChoice.Revised...")}
    # Revise the Split choice
    if(is.na(input$CellratioFacetChoice) | input$CellratioFacetChoice == "None") {
      return(NULL)
    }else{
      return(input$CellratioFacetChoice)
    }
  })

  # define Facet order
  output$CellratioplotFacetOrder.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing CellratioplotFacetOrder.UI...")}
    if (!is.null(FacetChoice.Revised())) {
      shinyjqui::orderInput(inputId = 'CellratioFacetOrder', label = 'Drag to order:', items = levels(data$obj@meta.data[,input$CellratioFacetChoice]),width = '100%')
    }else{

    }
  })

  cellratioplot_width  <- reactive({ session$clientData$output_cellratioplot_width})


  # plot
  output$cellratioplot <- renderPlot({
    req(input$CellratioXChoice)
    req(input$CellratioXOrder)
    req(input$CellratioFillChoice)
    req(input$CellratioFillOrder)
    if(verbose){message("SeuratExplorer: preparing cellratioplot...")}
    isolate(cds <- data$obj)
    if (is.null(FacetChoice.Revised())) { # not facet
      p <- cellRatioPlot(object = cds, sample.name = input$CellratioXChoice, sample.order = input$CellratioXOrder,
                         celltype.name = input$CellratioFillChoice, celltype.order = input$CellratioFillOrder,
                         facet.name = NULL, facet.order = NULL,
                         col.width = input$CellratioColumnWidth, flow.alpha = input$CellratioFlowAlpha,
                         flow.curve = input$CellratioFlowCurve, color.choice = input$fillcolorplatte)
    }else{
      p <- cellRatioPlot(object = cds, sample.name = input$CellratioXChoice, sample.order = input$CellratioXOrder,
                         celltype.name = input$CellratioFillChoice, celltype.order = input$CellratioFillOrder,
                         facet.name = FacetChoice.Revised(),
                         facet.order = input$CellratioFacetOrder,
                         col.width = input$CellratioColumnWidth, flow.alpha = input$CellratioFlowAlpha,
                         flow.curve = input$CellratioFlowCurve, color.choice = input$fillcolorplatte)
    }
    if (input$CellratioRotateAxis) { p <- p & ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) }
    ggplot2::ggsave(paste0(temp_dir,"/cellratioplot.pdf"), p, width = cellratioplot_width() * px2cm, height = cellratioplot_width() * input$CellratioplotHWRatio * px2cm, units = "cm", limitsize = FALSE)
    return(p)
  }, height = function(){session$clientData$output_cellratioplot_width * input$CellratioplotHWRatio}) # box plot: height = width default

  # download
  output$downloadcellratioplot <- downloadHandler(
    filename = function(){'cellratioplot.pdf'},
    content = function(file) {
      if (file.exists(paste0(temp_dir,"/cellratioplot.pdf"))) {
        file.copy(paste0(temp_dir,"/cellratioplot.pdf"), file, overwrite=TRUE)
      }
    })

  ################################ DEGs analysis
  # Warning
  output$degs_warning = renderText({
    paste0('Attention:
            This usually takes longer, please wait patiently. Make sure to save current results before a new analysis!')
  })
  # info
  output$degs_info = renderText({
    paste0('About:
            FindMarkers for All Clusters: calculate markers for all groups.
            Find DEGs for two groups: comparison between two groups, support subet cells before a comparison.')
  })

  DEGs <- reactiveValues(degs = NULL, degs_ready = FALSE)

  output$DEGs_ready <- reactive({
    return(DEGs$degs_ready)
  })

  outputOptions(output, 'DEGs_ready', suspendWhenHidden=FALSE)

  # Part-1: Cluster Markers

  # define Cluster Annotation choice
  output$ClusterMarkersClusterResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing ClusterMarkersClusterResolution.UI...")}
    selectInput("ClusterMarkersClusterResolution","Choose A Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  observeEvent(input$DEGsClusterMarkersAnalysis, {
    if(verbose){message("SeuratExplorer: preparing DEGsClusterMarkersAnalysis...")}
    isolate(cds <- data$obj)
    Seurat::Idents(cds) <- input$ClusterMarkersClusterResolution
    # if layers are splited

    if (length(as.character(Idents(cds))) < 2) {
      showModal(modalDialog(title = "Error...", "Please select a cluster resolution with more than one group!",easyClose = TRUE,footer = NULL, size = "l"))
    }else{
      showModal(modalDialog(title = "Calculating Cluster Markers...", "Please wait for a few minutes!", footer= NULL, size = "l"))

      # check_dependency(test = input$testuse) # Seurat will check dependency!
      cluster.markers <- Seurat::FindAllMarkers(cds, test.use = input$testuse, logfc.threshold = input$logfcthreshold,
                                                min.pct = input$minpct, min.diff.pct = ifelse(input$mindiffpct, input$mindiffpct, -Inf), only.pos = TRUE)
      removeModal()
      DEGs$degs <- cluster.markers
      DEGs$degs_ready <- TRUE
    }
  })

  # Part-2: Find DEGs for two groups
  # define Cluster Annotation choice
  output$IntraClusterDEGsCustomizedGroups.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing IntraClusterDEGsCustomizedGroups.UI...")}
    selectInput("IntraClusterDEGsCustomizedGroups","Group Cells By:", choices = data$cluster_options)
  })

  # define the idents used
  output$IntraClusterDEGsCustomizedGroupsCase.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    if(verbose){message("SeuratExplorer: preparing IntraClusterDEGsCustomizedGroupsCase.UI...")}
    selectInput("IntraClusterDEGsCustomizedGroupsCase","Choose Case Samples:", choices = levels(data$obj@meta.data[,input$IntraClusterDEGsCustomizedGroups]), multiple = TRUE)
  })

  # define the idents used
  output$IntraClusterDEGsCustomizedGroupsControl.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    req(input$IntraClusterDEGsCustomizedGroupsCase)
    if(verbose){message("SeuratExplorer: preparing IntraClusterDEGsCustomizedGroupsControl.UI...")}
    selectInput("IntraClusterDEGsCustomizedGroupsControl","Choose control Samples:", multiple = TRUE,
                choices = setdiff(levels(data$obj@meta.data[,input$IntraClusterDEGsCustomizedGroups]),input$IntraClusterDEGsCustomizedGroupsCase))
  })

  # define Cluster Annotation choice
  output$IntraClusterDEGsSubsetCells.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    if(verbose){message("SeuratExplorer: preparing IntraClusterDEGsSubsetCells.UI...")}
    selectInput("IntraClusterDEGsSubsetCells","Subset Cells By:", choices = setdiff(data$cluster_options, input$IntraClusterDEGsCustomizedGroups))
  })

  # define Cluster Annotation choice
  output$IntraClusterDEGsSubsetCellsSelectedClusters.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    req(input$IntraClusterDEGsSubsetCells)
    if(verbose){message("SeuratExplorer: preparing IntraClusterDEGsSubsetCellsSelectedClusters.UI...")}
    shinyWidgets::pickerInput(inputId = "IntraClusterDEGsSubsetCellsSelectedClusters", label = "Select Clusters:",
                              choices = levels(data$obj@meta.data[,input$IntraClusterDEGsSubsetCells]), selected = levels(data$obj@meta.data[,input$IntraClusterDEGsSubsetCells]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  # compare two groups, support subset clusters before comparison
  observeEvent(input$IntraClusterDEGssAnalysis, {
    if(verbose){message("SeuratExplorer: calculate DEGs...")}
    if (any(is.null(input$IntraClusterDEGsCustomizedGroupsCase), is.null(input$IntraClusterDEGsCustomizedGroupsControl), is.null(input$IntraClusterDEGsSubsetCellsSelectedClusters))) {
      showModal(modalDialog(title = "Error:","Please specify the case & control samples and clusters used. Press ESC to close.",easyClose = TRUE,footer = NULL))
    }else{
      showModal(modalDialog(title = "Calculating DEGs...", "Please wait for a few minutes!", footer= NULL, size = "l"))
      isolate(cds <- data$obj)
      Seurat::Idents(cds) <- input$IntraClusterDEGsSubsetCells
      cds <- subset_Seurat(cds, idents = input$IntraClusterDEGsSubsetCellsSelectedClusters)
      Seurat::Idents(cds) <- input$IntraClusterDEGsCustomizedGroups
      # check_dependency(test = input$testuse)
      cluster.markers <- Seurat::FindMarkers(cds, ident.1 = input$IntraClusterDEGsCustomizedGroupsCase, ident.2 = input$IntraClusterDEGsCustomizedGroupsControl,
                                             test.use = input$testuse, logfc.threshold = input$logfcthreshold,
                                             min.pct = input$minpct, min.diff.pct = ifelse(input$mindiffpct, input$mindiffpct, -Inf))
      removeModal()
      DEGs$degs <- cluster.markers
      DEGs$degs_ready <- TRUE
    }
  })

  # part-4: reset parameters
  observeEvent(input$SetDefault, {
    if(verbose){message("SeuratExplorer: reset DEGs parameters...")}
    updateSelectInput(session = session, inputId = "testuse", selected = "wilcox")
    updateSliderInput(session, "logfcthreshold", value = 0.1 )
    updateSliderInput(session, "minpct", value = 0.01 )
    updateSliderInput(session, "mindiffpct", value = 0 )
  })

  # part-5: output results
  output$dataset_degs <-  DT::renderDT(server=FALSE,{
    req(DEGs$degs)
    if(verbose){message("SeuratExplorer: preparing dataset_degs...")}
    # Show data
    if (nrow(DEGs$degs) == 0 | is.null(DEGs$degs)) {
      showModal(modalDialog(title = "Error", "None of DEGs found, contact technican for details!", footer= modalButton("Dismiss"), easyClose = TRUE, size = "l"))
      return(NULL)
    }else{
      data_res <- DT::datatable(DEGs$degs, extensions = 'Buttons',
                                options = list(scrollX=TRUE,
                                               # lengthMenu = c(5,10,15),
                                               paging = TRUE, searching = TRUE,
                                               fixedColumns = TRUE, autoWidth = TRUE,
                                               ordering = TRUE, dom = 'Bfrtip',
                                               buttons = list('copy',
                                                              list(extend = 'csv', title = "DEGs"),
                                                              list(extend = 'excel', title = "DEGs"))))
      for (acolumn in c("p_val","p_val_adj")) {
        if (acolumn %in% colnames(DEGs$degs)) {
          data_res <- DT::formatSignif(data_res, columns = acolumn, digits = 3)
        }
      }
      for (acolumn in c("avg_log2FC", "avg_diff", "avg_logFC")) {
        if (acolumn %in% colnames(DEGs$degs)) {
          data_res <- DT::formatRound(data_res, columns = acolumn, digits = 3)
        }
      }
      return(data_res)
    }
  })

  ################################ Top genes analysis
  # Warning
  output$topgenes_warning = renderText({
    paste0('Attention:
            This usually takes longer, please wait patiently. Save current results before a new analysis!')
  })
  # info
  output$topgenes_info = renderText({
    paste0('About:
            Find Top Genes by Cell: firstly, for each cell, find genes that has high UMI percentage, then summary those genes for each cluster, details see github page.
            Find Top Genes by mean UMI Counts: for each cluster, calculate the top n highly expressed genes by mean UMI counts.')
  })

  TopGenes <- reactiveValues(topgenes = NULL, topgenes_ready = FALSE)

  output$TopGenes_ready <- reactive({
    return(TopGenes$topgenes_ready)
  })

  outputOptions(output, 'TopGenes_ready', suspendWhenHidden=FALSE)

  # define Cluster Annotation choice
  output$TopGenesClusteResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing TopGenesClusteResolution.UI...")}
    selectInput("TopGenesClusterResolution","Choose A Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  observeEvent(input$TopGenesAnalysis, {
    if(verbose){message("SeuratExplorer: preparing TopGenesAnalysis...")}
    showModal(modalDialog(title = "Calculating Top Genes...", "Please wait for a few minutes!", footer= NULL, size = "l"))
    isolate(cds <- data$obj)
    TopGenes$topgenes <- top_genes(SeuratObj = cds, expr.cut = input$percentcut/100, group.by = input$TopGenesClusterResolution)
    removeModal()
    TopGenes$topgenes_ready <- TRUE
  })

  observeEvent(input$TopAccumulatedGenesAnalysis, {
    if(verbose){message("SeuratExplorer: preparing TopAccumulatedGenesAnalysis...")}
    showModal(modalDialog(title = "Calculating Top Genes...", "Please wait for a few minutes!", footer= NULL, size = "l"))
    isolate(cds <- data$obj)
    TopGenes$topgenes <- top_accumulated_genes(SeuratObj = cds, top = input$topcut, group.by = input$TopGenesClusterResolution)
    removeModal()
    TopGenes$topgenes_ready <- TRUE
  })

  output$dataset_topgenes <-  DT::renderDT(server=FALSE,{
    req(TopGenes$topgenes)
    if(verbose){message("SeuratExplorer: preparing topgenes...")}
    # Show data
    DT::datatable(TopGenes$topgenes, extensions = 'Buttons',
                  options = list(scrollX=TRUE,
                                 # lengthMenu = c(5,10,15),
                                 paging = TRUE, searching = TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip',
                                 buttons = list('copy',
                                                list(extend = 'csv', title = "top-features"),
                                                list(extend = 'excel', title = "top-features"))))
  })

  ################################ Feature Summary
  # info
  output$featuresummary_info = renderText({
    paste0('About:
            Summary interested features by cluster, such as the positive cell percentage and mean/median expression level.')
  })

  FeatureSummary <- reactiveValues(summary = NULL, summary_ready = FALSE)

  output$FeatureSummary_ready <- reactive({
    return(FeatureSummary$summary_ready)
  })

  outputOptions(output, 'FeatureSummary_ready', suspendWhenHidden=FALSE)

  # define Cluster Annotation choice
  output$FeatureSummaryClusteResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing FeatureSummaryClusteResolution.UI...")}
    selectInput("FeatureSummaryClusterResolution","Choose A Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  observeEvent(input$FeatureSummaryAnalysis, {
    if(verbose){message("SeuratExplorer: preparing FeatureSummaryAnalysis...")}
    if(is.na(input$FeatureSummarySymbol)){
      GeneRevised <- NA
    }else{
      GeneRevised <- CheckGene(InputGene = input$FeatureSummarySymbol, GeneLibrary =  rownames(data$obj))
    }
    if (any(is.na(GeneRevised))) {
      showModal(modalDialog(title = "Error", "None of the input genes can be found!", footer= modalButton("Dismiss"), easyClose = TRUE, size = "l"))
    }else{
      showModal(modalDialog(title = "Summarizing features...", "Please wait for a few minutes!", footer= NULL, size = "l"))
      isolate(cds <- data$obj)
      FeatureSummary$summary <- summary_features(SeuratObj = cds, features = GeneRevised, group.by = input$FeatureSummaryClusterResolution)
      removeModal()
      FeatureSummary$summary_ready <- TRUE
    }
  })

  output$dataset_featuresummary <-  DT::renderDT(server=FALSE,{
    req(FeatureSummary$summary)
    if(verbose){message("SeuratExplorer: preparing dataset_featuresummary...")}
    # Show data
    DT::datatable(FeatureSummary$summary, extensions = 'Buttons',
                  options = list(scrollX=TRUE,
                                 # lengthMenu = c(5,10,15),
                                 paging = TRUE, searching = TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip',
                                 buttons = list('copy',
                                                list(extend = 'csv', title = "feature-summary"),
                                                list(extend = 'excel', title = "feature-summary"))))
  })

  ################################ Feature Correlation
  # Warning
  output$featurecorrelation_warning = renderText({
    paste0('Attention:
            This usually takes longer, please wait patiently. Make sure to save current results before a new analysis!')
  })
  # info
  output$featurecorrelation_info = renderText({
    paste0('About:
            Find Top Correlated Gene Pairs: find top 1000 correlated gene pairs.
            Find Correlated Genes for A Gene: find the most correlated genes for input genes.
            Calculate Correlation for A Gene List: calculate the correlation value for each pair of the input genes.
     if nothing return, this is caused by the low expression of the input genes, very low expressed genes will be removed before analysis.')
  })

  FeatureCorrelation <- reactiveValues(summary = NULL, summary_ready = FALSE)

  output$FeatureCorrelation_ready <- reactive({
    return(FeatureCorrelation$summary_ready)
  })

  outputOptions(output, 'FeatureCorrelation_ready', suspendWhenHidden=FALSE)

  # define Cluster Annotation choice
  output$FeatureCorrelationClusteResolution.UI <- renderUI({
    if(verbose){message("SeuratExplorer: preparing FeatureCorrelationClusteResolution.UI...")}
    selectInput("FeatureCorrelationClusterResolution","Choose A Cluster Resolution:", choices = data$cluster_options, selected = data$cluster_default)
  })

  # define the idents used
  output$FeatureCorrelationIdentsSelected.UI <- renderUI({
    req(input$FeatureCorrelationClusterResolution)
    if(verbose){message("SeuratExplorer: preparing FeatureCorrelationIdentsSelected.UI...")}
    shinyWidgets::pickerInput(inputId = "FeatureCorrelationIdentsSelected", label = "Clusters Used:",
                              choices = levels(data$obj@meta.data[,input$FeatureCorrelationClusterResolution]), selected = levels(data$obj@meta.data[,input$FeatureCorrelationClusterResolution]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })



  observeEvent(input$TopCorrelationAnalysis, {
    if(verbose){message("SeuratExplorer: preparing TopCorrelationAnalysis...")}
    showModal(modalDialog(title = "Calculating", "Calculate top correlated gene pairs, which usually takes longer...", footer= NULL, size = "l"))
    isolate(cds <- data$obj)
    Seurat::Idents(cds) <- input$FeatureCorrelationClusterResolution
    cds <- subset_Seurat(cds, idents = input$FeatureCorrelationIdentsSelected)
    FeatureCorrelation$summary <- calculate_top_correlations(SeuratObj = cds, method = input$correlationmethod)
    removeModal()
    FeatureCorrelation$summary_ready <- TRUE
  })


  observeEvent(input$MostCorrelatedAnalysis, {
    if(verbose){message("SeuratExplorer: preparing MostCorrelatedAnalysis...")}
    feature.revised <- ReviseGene(Agene = trimws(input$MostCorrelatedAGene), GeneLibrary = rownames(data$obj))
    if(is.na(feature.revised)){
      showModal(modalDialog(title = "Error", "the input gene can not be found, please check...", footer= modalButton("Dismiss"), easyClose = TRUE, size = "l"))
    }else{
      showModal(modalDialog(title = "Calculating", "Calculate the most correlated genes for the input gene, which usually takes longer...", footer= NULL, size = "l"))
      isolate(cds <- data$obj)
      Seurat::Idents(cds) <- input$FeatureCorrelationClusterResolution
      cds <- subset_Seurat(cds, idents = input$FeatureCorrelationIdentsSelected)
      FeatureCorrelation$summary <- calculate_most_correlated(SeuratObj = cds, feature = feature.revised, method = input$correlationmethod)
      removeModal()
      FeatureCorrelation$summary_ready <- TRUE
    }
  })

  observeEvent(input$calculatecorrelation, {
    if(verbose){message("SeuratExplorer: preparing calculatecorrelation...")}
    if(is.na(input$CorrelationGeneList)){
      GeneRevised <- NA
    }else{
      GeneRevised <- CheckGene(InputGene = input$CorrelationGeneList, GeneLibrary =  rownames(data$obj))
    }
    if (any(is.na(GeneRevised))) {
      showModal(modalDialog(title = "Error", "None of the input genes can be found!", footer= modalButton("Dismiss"), easyClose = TRUE, size = "l"))
    }else if(length(GeneRevised) < 2){
      showModal(modalDialog(title = "Error", "Please input at least two genes!", footer= modalButton("Dismiss"), easyClose = TRUE, size = "l"))
    }else{
      showModal(modalDialog(title = "Calculating", "Calculate the correlation for the specified gene list...", footer= NULL, size = "l"))
      isolate(cds <- data$obj)
      Seurat::Idents(cds) <- input$FeatureCorrelationClusterResolution
      cds <- subset_Seurat(cds, idents = input$FeatureCorrelationIdentsSelected)
      FeatureCorrelation$summary <- calculate_correlation(SeuratObj = cds, features = GeneRevised, method = input$correlationmethod)
      removeModal()
      FeatureCorrelation$summary_ready <- TRUE
    }
  })

  output$dataset_correlation <-  DT::renderDT(server=FALSE,{
    req(FeatureCorrelation$summary)
    if(verbose){message("SeuratExplorer: preparing dataset_featuresummary...")}
    # Show data
    DT::datatable(FeatureCorrelation$summary, extensions = 'Buttons',
                  options = list(scrollX=TRUE,
                                 # lengthMenu = c(5,10,15),
                                 paging = TRUE, searching = TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip',
                                 buttons = list('copy',
                                                list(extend = 'csv', title = "feature-correlation"),
                                                list(extend = 'excel', title = "feature-correlation"))))
  })
}

#' Server
#' @import shiny shinydashboard shinyWidgets
#' @import ggplot2 Seurat SeuratObject
#' @importFrom utils write.csv
#'
#' @param input Input from the UI
#' @param output Output to send back to UI
#' @param session from shiny server function
#'
#' @export
#' @return the server functions of shiny app
#'
server <- function(input, output, session) {
  # set the limited upload file size
  options(shiny.maxRequestSize=5*1024^3)

  ## Dataset tab ----
  # reactiveValues: Create an object for storing reactive values,similar to a list,
  # but with special capabilities for reactive programming.
  data = reactiveValues(obj = NULL, loaded = FALSE, Name = NULL, Path = NULL, species = NULL,
                        reduction_options = NULL, reduction_default = NULL,
                        cluster_options = NULL, cluster_default = NULL,
                        split_maxlevel = 6, split_options = NULL,
                        extra_qc_options = NULL)

  # reductions_options: xy axis coordinate
  # cluster_options/split_options/extra_qc_options all are column name from seurat object meta.data, which will be used for later plot
  # load data after data selection
  observe({
    shiny::req(input$dataset_file) # req: Check for required values; dataset_file is a data.frame
    ext = tools::file_ext(input$dataset_file$datapath) # file_ext: returns the file (name) extensions
    # validate + need: check filename postfix, in not rds or qs2, will throw an error
    validate(need(expr = ext %in% c("rds","qs2","Rds"), message = "Please upload a .rds or a .qs2 file"))
    data$obj <- prepare_seurat_object(obj = readSeurat(path = input$dataset_file$datapath), verbose = getOption('SeuratExplorerVerbose'))
    data$reduction_options <- prepare_reduction_options(obj = data$obj, keywords = c("umap","tsne"), verbose = getOption('SeuratExplorerVerbose'))
    data$cluster_options <- prepare_cluster_options(df = data$obj@meta.data, verbose = getOption('SeuratExplorerVerbose'))
    data$split_options <- prepare_split_options(df = data$obj@meta.data, max.level = data$split_maxlevel, verbose = getOption('SeuratExplorerVerbose'))
    data$extra_qc_options <- prepare_qc_options(df = data$obj@meta.data, types = c("double","integer","numeric"), verbose = getOption('SeuratExplorerVerbose'))
  })

  # after data loaded,set loaded to TRUE
  observe({
    req(data$obj)
    data$loaded = !is.null(data$obj)
  })

  ############################### Render metadata table
  # Server set to TRUE: https://stackoverflow.com/questions/50039186/add-download-buttons-in-dtrenderdatatable
  # when sever is set to TRUE, to download the whole data in DT button extensions.https://github.com/rstudio/DT/issues/267
  output$dataset_meta <- DT::renderDT(server=TRUE,{
    req(data$obj)
    # Show data
    DT::datatable(data$obj@meta.data,
                  callback = DT::JS("$('div.dwnld').append($('#download_meta_data'));"),
                  extensions = 'Buttons',
                  options = list(scrollX=TRUE,
                                 # lengthMenu = c(5,10,15),
                                 #paging = TRUE,
                                 #searching = TRUE,
                                 #fixedColumns = TRUE,
                                 #autoWidth = TRUE,
                                 ordering = TRUE,
                                 dom = 'B<"dwnld">frtip',
                                 buttons = list('copy')
                  ))
  })

  output$download_meta_data <- downloadHandler(
    filename = function() {
      "cell-metadata.csv"
    },
    content = function(file) {
      write.csv(data$obj@meta.data, file)
    }
  )

  # Conditional panel control based on loaded obj, after loaded, show other UIs
  output$file_loaded = reactive({
    return(data$loaded)
  })

  outputOptions(output, 'file_loaded', suspendWhenHidden=FALSE)

  # Seurat Explorer functions
  explorer_server(input = input, output = output, session = session, data = data, verbose = getOption('SeuratExplorerVerbose'))

}
