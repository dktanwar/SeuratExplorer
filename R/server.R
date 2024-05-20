# server.R
## R shiny server side for SeuratExplorer

#' Server for SeuratExplorer shiny app
#' @import shiny
#' @import Seurat SeuratObject
#' @param input Input from the UI
#' @param output Output to send back to UI
#' @param session from shiny server function
#' @export
server <- function(input, output, session) {
  requireNamespace("Seurat")
  requireNamespace("ggplot2")
  requireNamespace("shinyWidgets")
  requireNamespace("shinydashboard")
  requireNamespace("SeuratObject")

  # 设置上传文件的大小限制
  options(shiny.maxRequestSize=5*1024^3)

  ## Dataset tab ----
  # reactiveValues: Create an object for storing reactive values,similar to a list,
  # but with special capabilities for reactive programming.
  data = reactiveValues(obj = NULL, loaded = FALSE, reduction_options = NULL, cluster_options = NULL, split_options = NULL, extra_qc_options = NULL)
  # reductions_options: 为可视化时的xy轴座标
  # cluster_options/split_options/extra_qc_options均为seurat object meta.data里的列名，后续绘图会作为可选参数被反复用到
  # 选择好数据后，读入数据
  observe({
    shiny::req(input$dataset_file) # req: Check for required values; dataset_file is a data.frame
    ext = tools::file_ext(input$dataset_file$datapath) # file_ext: returns the file (name) extensions
    validate(need(expr = ext == "rds", message = "Please upload a .rds file")) # validate + need：检查后缀是否为rds，否则抛出错误
    data$obj <- prepare_seurat_object(obj = Seurat::UpdateSeuratObject(readRDS(file = input$dataset_file$datapath)))
    data$reduction_options <- prepare_reduction_options(obj = data$obj, keywords = c("umap","tsne"))
    data$cluster_options <- prepare_cluster_options(df = data$obj@meta.data)
    data$split_options <- prepare_split_options(df = data$obj@meta.data, max.level = 6)
    data$extra_qc_options <- prepare_qc_options(df = data$obj@meta.data, types = c("double","integer","numeric"))
  })

  # 数据加载成功后，设置loaded为TRUE
  observe({
    req(data$obj)
    data$loaded = !is.null(data$obj)
  })

  ############################### Render metadata table
  # 可以下载全部，参考：https://stackoverflow.com/questions/50039186/add-download-buttons-in-dtrenderdatatable
  output$dataset_meta <- DT::renderDT(server=FALSE,{
    req(data$obj)
    # Show data
    DT::datatable(data$obj@meta.data, extensions = 'Buttons',
              options = list(scrollX=TRUE, lengthMenu = c(5,10,15),
                             paging = TRUE, searching = TRUE,
                             fixedColumns = TRUE, autoWidth = TRUE,
                             ordering = TRUE, dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel')))
  })

  # Conditional panel control based on loaded obj，条件panel,数据记载成功后，显示：dashboardSidebar -sidebarMenu - menuItem - Explorer和 dashboardBody - dataset - tabItem -  box - Cell Meta Info
  output$file_loaded = reactive({
    return(data$loaded)
  })

  # Disable suspend for output$file_loaded, 当被隐藏时，禁用暂停，conditionalpanel所需要要的参数
  outputOptions(output, 'file_loaded', suspendWhenHidden=FALSE)

  ############################# Dimension Reduction Plot
  # define reductions choices UI
  output$DimReductions.UI <- renderUI({
      selectInput("DimDimensionReduction", "Dimension Reduction:", choices = data$reduction_options) # set default reduction
  })

  # define Cluster Annotation choice
  output$DimClusterResolution.UI <- renderUI({
    selectInput("DimClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  })

  # define Split Choice UI
  output$DimSplit.UI <- renderUI({
    selectInput("DimSplit","Split by:", choices = c("None" = "None", data$split_options))
  })

  # Revise Split selection which will be appropriate for plot
  DimSplit.Revised <- reactive({
    req(input$DimSplit) # split值出现后，才会执行的代码
    # Revise the Split choice
    if(is.na(input$DimSplit) | input$DimSplit == "None") {
      return(NULL)
    }else{
      return(input$DimSplit)
    }
  })

  output$dimplot <- renderPlot({
    if (is.null(DimSplit.Revised())) { # not splited
      Seurat::DimPlot(data$obj, reduction = input$DimDimensionReduction, label = input$DimShowLabel, pt.size = input$DimPointSize, label.size = input$DimLabelSize,
                      group.by = input$DimClusterResolution)
    }else{ # splited
      plot_numbers <- length(levels(data$obj@meta.data[,DimSplit.Revised()]))
      Seurat::DimPlot(data$obj, reduction = input$DimDimensionReduction, label = input$DimShowLabel, pt.size = input$DimPointSize, label.size = input$DimLabelSize,
                      group.by = input$DimClusterResolution, split.by = DimSplit.Revised(), ncol = ceiling(sqrt(plot_numbers)))
    }
  }, height = function(){session$clientData$output_dimplot_width * input$DimPlotHWRatio}) # box plot: height = width default


  ################################ Feature Plot
  # define reductions choices UI
  output$FeatureReductions.UI <- renderUI({
    selectInput("FeatureDimensionReduction", "Dimension Reduction:", choices = data$reduction_options) # set default reduction
  })

  # # define Cluster Annotation choice
  # output$FeatureClusterResolution.UI <- renderUI({
  #   selectInput("FeatureClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  # })

  # define Split Choice UI
  output$FeatureSplit.UI <- renderUI({
    selectInput("FeatureSplit","Split by:", choices = c("None" = "None", data$split_options))
  })

  # 提示可用的qc选项作为Gene symbol
  output$Featurehints.UI <- renderUI({
    helpText(strong(paste("Multiple genes are separted by a comma, such as: CD4,CD8A; Also supports: ", paste(data$extra_qc_options, collapse = ", "), ".",sep = "")),style = "font-size:12px;")
  })


  # Revise Split selection which will be appropriate for DimPlot, FeaturePlot and Vlnplot functions.
  FeatureSplit.Revised <- reactive({
    req(input$FeatureSplit) # split值出现后，才会执行的代码
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
    ifelse(is.na(input$FeatureGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$FeatureGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  output$featureplot <- renderPlot({
    if (any(is.na(Featureplot.Gene.Revised()))) { # NA 值时
      ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else if(is.null(FeatureSplit.Revised())) { # not splited
      Seurat::FeaturePlot(data$obj, features = Featureplot.Gene.Revised(), pt.size = input$FeaturePointSize, reduction = input$FeatureDimensionReduction,
                          cols = c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor))
     }else{ # splited
      p <- Seurat::FeaturePlot(data$obj, features = Featureplot.Gene.Revised(), pt.size = input$FeaturePointSize, reduction = input$FeatureDimensionReduction,
                          cols =  c(input$FeaturePlotLowestExprColor,input$FeaturePlotHighestExprColor), split.by = FeatureSplit.Revised())
      if (length( Featureplot.Gene.Revised()) == 1) { # 仅仅一个基因时
        plot_numbers <- length(levels(data$obj@meta.data[,FeatureSplit.Revised()]))
        p + patchwork::plot_layout(ncol = ceiling(sqrt(plot_numbers)),nrow = ceiling(plot_numbers/ceiling(sqrt(plot_numbers))))
      }else{ # 多个基因时
        p
      }
    }
  }, height = function(){session$clientData$output_featureplot_width * input$FeaturePlotHWRatio}) # box plot: height = width default

  ################################ Violin Plot
  # Check the input gene
  Vlnplot.Gene.Revised <- reactive({
    req(input$VlnGeneSymbol)
    ifelse(is.na(input$VlnGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$VlnGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  # 提示可用的qc选项作为Gene symbol
  output$Vlnhints.UI <- renderUI({
    helpText(strong(paste("Multiple genes are separted by a comma, such as: CD4,CD8A; Also supports: ", paste(data$extra_qc_options, collapse = ", "), ".",sep = "")),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$VlnClusterResolution.UI <- renderUI({
    selectInput("VlnClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  })

  # define Split Choice UI
  output$VlnSplitBy.UI <- renderUI({
    selectInput("VlnSplitBy","Split by:", choices = c("None" = "None", data$split_options))
  })

  # Conditional panel: split.by被勾选，且level数目为2时，显示此panel
  output$Vlnplot_splitoption_twolevels = reactive({
    req(input$VlnSplitBy)
    if (input$VlnSplitBy == "None"){
      return(FALSE)
    }else if(length(levels(data$obj@meta.data[,input$VlnSplitBy])) == 2) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  # Disable suspend for output$file_loaded, 当被隐藏时，禁用暂停，conditionalpanel所需要要的参数
  outputOptions(output, 'Vlnplot_splitoption_twolevels', suspendWhenHidden = FALSE)

  # Conditional panel: symbol输入为多个基因时，显示此panel
  output$Vlnplot_multiple_genes = reactive({
    req(input$VlnGeneSymbol)
    if (length(Vlnplot.Gene.Revised()) > 1) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  # Disable suspend for output$file_loaded, 当被隐藏时，禁用暂停，conditionalpanel所需要要的参数
  outputOptions(output, 'Vlnplot_multiple_genes', suspendWhenHidden = FALSE)


  # Conditional panel: 输入为多个基因且stack设为TRUE时，显示此panel
  output$Vlnplot_StackPlot = reactive({
    req(input$VlnStackPlot)
    req(input$VlnGeneSymbol)
    if (length(Vlnplot.Gene.Revised()) > 1 & input$VlnStackPlot) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  # Disable suspend for output$file_loaded, 当被隐藏时，禁用暂停，conditionalpanel所需要要的参数
  outputOptions(output, 'Vlnplot_StackPlot', suspendWhenHidden = FALSE)

  # Revise Split selection which will be appropriate for DimPlot, FeaturePlot and Vlnplot functions.
  VlnSplit.Revised <- reactive({
    req(input$VlnSplitBy) # split值出现后，才会执行的代码
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
    updateCheckboxInput(session, "VlnSplitPlot", value = FALSE)
    updateCheckboxInput(session, "VlnStackPlot", value = FALSE)
    updateCheckboxInput(session, "VlnFlipPlot", value = FALSE)
    updateSelectInput(session, "VlnFillBy", selected = "feature")
  })

  # shiny related bug
  # 以后再debug吧！ 2024.05.15
  # how to make sure renderPlot run after the observe(input$VlnSplitBy)[Warning: Error in SingleExIPlot: Unknown plot type: splitViolin,
  # 因为此时的VlnSplitPlot还没有被更新]

  # seurat related bug
  # VlnPlot(cds,features = c("CD4","CD8A"),split.by = "orig.ident", stack = TRUE,group.by = "cca_clusters_res_0.2",flip = FALSE,split.plot = TRUE)
  # Error:
  # Error in `vln.geom()`:
  #   ! Problem while converting geom to grob.
  # ℹ Error occurred in the 1st layer.
  # Caused by error in `$<-.data.frame`:
  #   ! 替换数据里有0行，但数据有512
  # Run `rlang::last_trace()` to see where the error occurred
  # 经测试，与ggplot2,pathcwork,rlang等R包的版本无关！

  output$vlnplot <- renderPlot({
    req(Vlnplot.Gene.Revised())
    if (any(is.na(Vlnplot.Gene.Revised()))) { # NA 值时
      ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else if(length(Vlnplot.Gene.Revised()) == 1) { # only One Gene
      Seurat::VlnPlot(data$obj, features = Vlnplot.Gene.Revised(), group.by = input$VlnClusterResolution,
                      split.by = VlnSplit.Revised(), split.plot = input$VlnSplitPlot, pt.size = input$VlnPointSize, alpha = input$VlnPointAlpha) &
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                       axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize))
    }else{ # multiple genes
      Seurat::VlnPlot(data$obj, features = Vlnplot.Gene.Revised(), group.by = input$VlnClusterResolution,
                      split.by = VlnSplit.Revised(), split.plot = input$VlnSplitPlot, stack = input$VlnStackPlot,
                      flip = input$VlnFlipPlot, fill.by = input$VlnFillBy,
                      pt.size = input$VlnPointSize, alpha = input$VlnPointAlpha) &
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$VlnXlabelSize),
                       axis.text.y = ggplot2::element_text(size = input$VlnYlabelSize))
    }
  }, height = function(){session$clientData$output_vlnplot_width * input$VlnPlotHWRatio}) # box plot: height = width default

  ################################ Dot Plot

  # Check the input gene
  Dotplot.Gene.Revised <- reactive({
    req(input$DotGeneSymbol)
    ifelse(is.na(input$DotGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$DotGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  # 提示可用的qc选项作为Gene symbol
  output$Dothints.UI <- renderUI({
    helpText(strong(paste("Multiple genes are separted by a comma, such as: CD4,CD8A; Also supports: ", paste(data$extra_qc_options, collapse = ", "), ".",sep = "")),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$DotClusterResolution.UI <- renderUI({
    selectInput("DotClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  })

  # define the idents used
  output$DotIdentsSelected.UI <- renderUI({
    req(input$DotClusterResolution)
    # selectInput("DotIdentsSelected","Idents used:", choices = levels(data$obj@meta.data[,input$DotClusterResolution]))
    shinyWidgets::pickerInput(inputId = "DotIdentsSelected", label = "Idents used:",
      choices = levels(data$obj@meta.data[,input$DotClusterResolution]), selected = levels(data$obj@meta.data[,input$DotClusterResolution]),
      options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  # define Split Choice UI
  output$DotSplitBy.UI <- renderUI({
    selectInput("DotSplitBy","Split by:", choices = c("None" = "None", data$split_options))
  })

  # Revise Split selection which will be appropriate for DimPlot, FeaturePlot and Vlnplot functions.
  DotSplit.Revised <- reactive({
    req(input$DotSplitBy) # split值出现后，才会执行的代码
    # Revise the Split choice
    if(is.na(input$DotSplitBy) | input$DotSplitBy == "None") {
      return(NULL)
    }else{
      return(input$DotSplitBy)
    }
  })

  # Conditional panel: 当split为NULL时，可以自行设定最高和最低表达值对应的颜色。当split不为NULL时，需要软件自动使用ggplot2生成的颜色填充每组的点的颜色。
  output$DotPlot_Split_isNone <- reactive({
    req(input$DotSplitBy)
    if(is.na(input$DotSplitBy) | input$DotSplitBy == "None") {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  # Disable suspend for output$file_loaded, 当被隐藏时，禁用暂停，conditional panel所需要要的参数
  outputOptions(output, 'DotPlot_Split_isNone', suspendWhenHidden = FALSE)

  output$dotplot <- renderPlot({
    req(Dotplot.Gene.Revised())

    if (any(is.na(Dotplot.Gene.Revised()))) { # NA 值时
      ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj) # 不是一个优雅的做法，会使用额外的内存资源，另一个坏处是，可能对data$obj不在实时有反应？
      Idents(cds) <- input$DotClusterResolution

      if (is.null(DotSplit.Revised())) {
        p <- Seurat::DotPlot(cds, features = Dotplot.Gene.Revised(), group.by = input$DotClusterResolution,
                             idents = input$DotIdentsSelected, #不支持使用Group.by参数所定义的cluster，所以需要新见变量cds，修改Idents
                             split.by = DotSplit.Revised(), cluster.idents = input$DotClusterIdents, dot.scale = input$DotDotScale,
                             cols = c(input$DotPlotLowestExprColor, input$DotPlotHighestExprColor))
      }else{
        split.levels.length <- length(levels(cds@meta.data[,DotSplit.Revised()]))
        p <- Seurat::DotPlot(cds, features = Dotplot.Gene.Revised(), group.by = input$DotClusterResolution,
                             idents = input$DotIdentsSelected, #不支持使用Group.by参数所定义的cluster，所以需要新见变量cds，修改Idents
                             split.by = DotSplit.Revised(), cluster.idents = input$DotClusterIdents, dot.scale = input$DotDotScale,
                             cols = scales::hue_pal()(split.levels.length))
      }
      p <- p & ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$DotXlabelSize), axis.text.y = ggplot2::element_text(size = input$DotYlabelSize))
      if (input$DotRotateAxis) { p <- p + Seurat::RotatedAxis() }
      if (input$DotFlipCoordinate) { p <- p + ggplot2::coord_flip() }
      p
    }
  }, height = function(){session$clientData$output_dotplot_width * input$DotPlotHWRatio}) # box plot: height = width default


  ################################ Heatmap

  # Check the input gene
  Heatmap.Gene.Revised <- reactive({
    req(input$HeatmapGeneSymbol)
    ifelse(is.na(input$HeatmapGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$HeatmapGeneSymbol, GeneLibrary =  rownames(data$obj))))
  })

  # 提示可用的qc选项作为Gene symbol
  output$Heatmaphints.UI <- renderUI({
    helpText(strong("Multiple genes are separted by a comma, such as: CD4,CD8A."),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$HeatmapClusterResolution.UI <- renderUI({
    selectInput("HeatmapClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  })


  output$heatmap <- renderPlot({
    req(Heatmap.Gene.Revised())

    if (any(is.na(Heatmap.Gene.Revised()))) { # NA 值时
      ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      isolate(cds <- data$obj) # 不是一个优雅的做法，会使用额外的内存资源，另一个坏处是，可能对data$obj不在实时有反应？
      if (!all(Heatmap.Gene.Revised() %in% Seurat::VariableFeatures(cds))) { #问题： 每次修改任何绘图参数，都得要执行此代码！！！
        cds <- Seurat::ScaleData(object = cds, features = unique(c(Seurat::VariableFeatures(cds), Heatmap.Gene.Revised()))) # 不能只用一个基因去做scaledata，会报错。
      }

      Seurat::DoHeatmap(object = cds, features = Heatmap.Gene.Revised(), group.by = input$HeatmapClusterResolution, size = input$HeatmapTextSize,
                        hjust = input$HeatmapTextHjust, vjust = input$HeatmapTextVjust, angle = input$HeatmapTextRatateAngle,
                        group.bar.height = input$HeatmapGroupBarHeight, lines.width = input$HeatmapLineWidth) &
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = input$HeatmapFeatureTextSize))
    }
  }, height = function(){session$clientData$output_heatmap_width * input$HeatmapPlotHWRatio}) # box plot: height = width default

  ################################ Ridge Plot

  # Check the input gene
  Ridgeplot.Gene.Revised <- reactive({
    req(input$RidgeplotGeneSymbol)
    ifelse(is.na(input$RidgeplotGeneSymbol), yes = return(NA), no = return(CheckGene(InputGene = input$RidgeplotGeneSymbol, GeneLibrary =  c(rownames(data$obj), data$extra_qc_options))))
  })

  # 提示可用的qc选项作为Gene symbol
  output$Ridgeplothints.UI <- renderUI({
    helpText(strong(paste("Multiple genes are separted by a comma, such as: CD4,CD8A; Also supports: ", paste(data$extra_qc_options, collapse = ", "), ".",sep = "")),style = "font-size:12px;")
  })

  # define Cluster Annotation choice
  output$RidgeplotClusterResolution.UI <- renderUI({
    selectInput("RidgeplotClusterResolution","Cluster Resolution:", choices = data$cluster_options)
  })

  # Conditional panel: 输入为多个基因且stack设为TRUE时，显示此panel
  output$Ridgeplot_stack_show = reactive({
    req(input$RidgeplotGeneSymbol)
    if (length(Ridgeplot.Gene.Revised()) > 1) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })

  outputOptions(output, 'Ridgeplot_stack_show', suspendWhenHidden = FALSE)

  # Conditional panel: 输入为多个基因且stack设为TRUE时，显示此panel
  output$Ridgeplot_stack_NotSelected = reactive({
    if (input$RidgeplotStackPlot) {
      return(FALSE)
    }else{
      return(TRUE)
    }
  })

  outputOptions(output, 'Ridgeplot_stack_NotSelected', suspendWhenHidden = FALSE)

  # reset VlnSplitPlot value to FALSE when change the split options
  observe({
    req(input$RidgeplotGeneSymbol)
    updateCheckboxInput(session, "RidgeplotStackPlot", value = FALSE)
  })

  output$ridgeplot <- renderPlot({
    req(Ridgeplot.Gene.Revised())
    if (any(is.na(Ridgeplot.Gene.Revised()))) { # NA 值时
      ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::geom_blank() # when no symbol or wrong input, show a blank pic.
    }else{
      Seurat::RidgePlot(object = data$obj, features = Ridgeplot.Gene.Revised(), group.by = input$RidgeplotClusterResolution, ncol = input$RidgeplotNumberOfColumns,
                        stack = input$RidgeplotStackPlot, fill.by = input$RidgeplotFillBy) &
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = input$RidgeplotXlabelSize),
                       axis.text.y = ggplot2::element_text(size = input$RidgeplotYlabelSize))
    }
  }, height = function(){session$clientData$output_ridgeplot_width * input$RidgeplotHWRatio}) # box plot: height = width default

  ################################ DEGs analysis
  # Warning
  output$degs_warning = renderText({
    paste0('Differential Expression testing is computationally intensive and may take some time. 注意，请在点击"Analyze"之前，保存好先前的分析结果！')
  })

  DEGs <- reactiveValues(degs = NULL, degs_ready = FALSE)

  output$DEGs_ready <- reactive({
    return(DEGs$degs_ready)
  })

  # Disable suspend for output$file_loaded, 当被隐藏时，禁用暂停，conditionalpanel所需要要的参数
  outputOptions(output, 'DEGs_ready', suspendWhenHidden=FALSE)

  # Part-1: ClusterMarkers

  # define Cluster Annotation choice
  output$ClusterMarkersClusterResolution.UI <- renderUI({
    selectInput("ClusterMarkersClusterResolution","Choose A Cluster Resolution:", choices = data$cluster_options)
  })

  observeEvent(input$DEGsClusterMarkersAnalysis, {
    showModal(modalDialog(title = "Calculating Cluster Markers...", "Please wait for a few minutes!", footer= NULL, size = "l"))
    isolate(cds <- data$obj)
    Seurat::Idents(cds) <- input$ClusterMarkersClusterResolution
    check_dependency(test = input$testuse)
    cluster.markers <- Seurat::FindAllMarkers(cds, test.use = input$testuse, logfc.threshold = input$logfcthreshold,
                                      min.pct = input$minpct, min.diff.pct = ifelse(input$mindiffpct, input$mindiffpct, -Inf), only.pos = TRUE)
    removeModal()
    DEGs$degs <<- cluster.markers #修改全局变量，需不需要改为 <<-
    DEGs$degs_ready <<- TRUE
  })

  # 功能冗余
  # # Part-2: InterClusterDEGs
  #
  # # define Cluster Annotation choice
  # output$InterClusterDEGsClusterResolution.UI <- renderUI({
  #   selectInput("InterClusterDEGsClusterResolution","Choose A Cluster Resolution:", choices = data$cluster_options)
  # })
  #
  # # define the idents used
  # output$InterClusterDEGsGroupCase.UI <- renderUI({
  #   req(input$InterClusterDEGsClusterResolution)
  #   selectInput("InterClusterDEGsGroupCase","Choose Case Clusters:", choices = levels(data$obj@meta.data[,input$InterClusterDEGsClusterResolution]), multiple = TRUE)
  # })
  #
  # # define the idents used
  # output$InterClusterDEGsGroupControl.UI <- renderUI({
  #   req(input$InterClusterDEGsClusterResolution)
  #   req(input$InterClusterDEGsGroupCase)
  #   selectInput("InterClusterDEGsGroupControl","Choose control Clusters:", multiple = TRUE,
  #               choices = setdiff(levels(data$obj@meta.data[,input$InterClusterDEGsClusterResolution]),input$InterClusterDEGsGroupCase))
  # })
  #
  #
  # observeEvent(input$InterClusterDEGsAnalysis, {
  #   if (any(is.null(input$InterClusterDEGsGroupCase), is.null(input$InterClusterDEGsGroupControl))) {
  #     showModal(modalDialog(title = "Error:","Please specify the case and control clusters. Press ESC to close.",easyClose = TRUE,footer = NULL))
  #   }else{
  #     showModal(modalDialog(title = "Calculating Cluster Markers...", "Please wait for a few minutes!", footer= NULL, size = "l"))
  #     isolate(cds <- data$obj)
  #     Seurat::Idents(cds) <- input$InterClusterDEGsClusterResolution
  #     check_dependency(test = input$testuse)
  #     cluster.markers <- Seurat::FindMarkers(cds, ident.1 = input$InterClusterDEGsGroupCase, ident.2 = input$InterClusterDEGsGroupControl,
  #                                            test.use = input$testuse, logfc.threshold = input$logfcthreshold,
  #                                            min.pct = input$minpct, min.diff.pct = ifelse(input$mindiffpct, input$mindiffpct, -Inf))
  #     removeModal()
  #     DEGs$degs <<- cluster.markers #修改全局变量，需不需要改为 <<-
  #     DEGs$degs_ready <<- TRUE
  #   }
  # })


  # Part-3: IntraClusterDEGs
  # define Cluster Annotation choice
  output$IntraClusterDEGsCustomizedGroups.UI <- renderUI({
    selectInput("IntraClusterDEGsCustomizedGroups","Group samples by:", choices = data$cluster_options)
  })

  # define the idents used
  output$IntraClusterDEGsCustomizedGroupsCase.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    selectInput("IntraClusterDEGsCustomizedGroupsCase","Choose Case Samples:", choices = levels(data$obj@meta.data[,input$IntraClusterDEGsCustomizedGroups]), multiple = TRUE)
  })

  # define the idents used
  output$IntraClusterDEGsCustomizedGroupsControl.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    req(input$IntraClusterDEGsCustomizedGroupsCase)
    selectInput("IntraClusterDEGsCustomizedGroupsControl","Choose control Samples:", multiple = TRUE,
                choices = setdiff(levels(data$obj@meta.data[,input$IntraClusterDEGsCustomizedGroups]),input$IntraClusterDEGsCustomizedGroupsCase))
  })

  # define Cluster Annotation choice
  output$IntraClusterDEGsSubsetCells.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    selectInput("IntraClusterDEGsSubsetCells","Subset Cells By:", choices = setdiff(data$cluster_options, input$IntraClusterDEGsCustomizedGroups))
  })

  # define Cluster Annotation choice
  output$IntraClusterDEGsSubsetCellsSelectedClusters.UI <- renderUI({
    req(input$IntraClusterDEGsCustomizedGroups)
    req(input$IntraClusterDEGsSubsetCells)
    shinyWidgets::pickerInput(inputId = "IntraClusterDEGsSubsetCellsSelectedClusters", label = "Select Clusters:",
                              choices = levels(data$obj@meta.data[,input$IntraClusterDEGsSubsetCells]), selected = levels(data$obj@meta.data[,input$IntraClusterDEGsSubsetCells]),
                              options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10, selectedTextFormat = "count > 3"), multiple = TRUE)
  })

  # 计算自定义分组的差异基因，支持subset clusters
  observeEvent(input$IntraClusterDEGssAnalysis, {
    if (any(is.null(input$IntraClusterDEGsCustomizedGroupsCase), is.null(input$IntraClusterDEGsCustomizedGroupsControl), is.null(input$IntraClusterDEGsSubsetCellsSelectedClusters))) {
      showModal(modalDialog(title = "Error:","Please specify the case & control samples and clusters used. Press ESC to close.",easyClose = TRUE,footer = NULL))
    }else{
      showModal(modalDialog(title = "Calculating DEGs...", "Please wait for a few minutes!", footer= NULL, size = "l"))
      isolate(cds <- data$obj)
      Seurat::Idents(cds) <- input$IntraClusterDEGsSubsetCells
      cds <- SeuratObject:::subset.Seurat(cds, idents = input$IntraClusterDEGsSubsetCellsSelectedClusters)
      Seurat::Idents(cds) <- input$IntraClusterDEGsCustomizedGroups
      check_dependency(test = input$testuse)
      cluster.markers <- Seurat::FindMarkers(cds, ident.1 = input$IntraClusterDEGsCustomizedGroupsCase, ident.2 = input$IntraClusterDEGsCustomizedGroupsControl,
                                             test.use = input$testuse, logfc.threshold = input$logfcthreshold,
                                             min.pct = input$minpct, min.diff.pct = ifelse(input$mindiffpct, input$mindiffpct, -Inf))
      removeModal()
      DEGs$degs <<- cluster.markers #修改全局变量，需不需要改为 <<-
      DEGs$degs_ready <<- TRUE
    }
  })

 # part-4: 重置参数
  observeEvent(input$SetDefault, {
    updateSelectInput(session = session, inputId = "testuse", selected = "wilcox")
    updateSliderInput(session, "logfcthreshold", value = 0.1 )
    updateSliderInput(session, "minpct", value = 0.01 )
    updateSliderInput(session, "mindiffpct", value = 0 )
  })

  # part-5: 输出结果
  output$dataset_degs <-  DT::renderDT(server=FALSE,{
    req(DEGs$degs)
    # Show data
    DT::datatable(DEGs$degs, extensions = 'Buttons',
                  options = list(scrollX=TRUE, lengthMenu = c(5,10,15),
                                 paging = TRUE, searching = TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel')))
  })

}
