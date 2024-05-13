# server.R
## R shiny server side for SeuratExplorer

#' Server for SeuratExplorer shiny app
#' @import shiny
#' @param input Input from the UI
#' @param output Output to send back to UI
#' @param session from shiny server function
#' @export
server <- function(input, output, session) {
  requireNamespace("Seurat")

  # 设置上传文件的大小限制
  options(shiny.maxRequestSize=5*1024^3)

  ## Dataset tab ----
  # reactiveValues: Create an object for storing reactive values,similar to a list,
  # but with special capabilities for reactive programming.
  data = reactiveValues(obj = NULL, loaded=FALSE)

  # 选择好数据后，读入数据
  observe({
    shiny::req(input$dataset_file) # req: Check for required values; dataset_file is a data.frame
    ext = tools::file_ext(input$dataset_file$datapath) # file_ext: returns the file (name) extensions
    validate(need(expr = ext == "rds", message = "Please upload a .rds file")) # validate + need：检查后缀是否为rds，否则抛出错误
    data$obj <- prepare_seurat_object(obj = Seurat::UpdateSeuratObject(readRDS(file = input$dataset_file$datapath)))
  })

  # 数据加载成功后，设置loaded为TRUE
  observe({
    req(data$obj)
    data$loaded = !is.null(data$obj)
  })

  # Render metadata table
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

  # define reductions choices UI
  output$Reductions.UI <- renderUI({
    reductions.keywords <- c("umap","tsne")
    reduction.choice <- grep(paste0(paste0("(",reductions.keywords,")"),collapse = "|"), Seurat::Reductions(data$obj),value = TRUE)
    names(reduction.choice) <- toupper(reduction.choice)
    selectInput("DimensionReduction", "Dimension Reduction:", choices = reduction.choice) # set default reduction
  })

  # define Cluster Annotation choice
  output$ClusterResolution.UI <- renderUI({
    full.res <- colnames(data$obj@meta.data)[check_df_factor(data$obj@meta.data)]
    names(full.res) <- full.res
    selectInput("ClusterResolution","Cluster Resolution:", choices = full.res)
  })

  # define Split Choice UI
  output$Split.UI <- renderUI({
    SplitBy.levels.max <- 4 #最大的样本数目
    SplitBy.Choice <- df_factor_columns(data$obj@meta.data, max.level = SplitBy.levels.max)
    SplitBy.Choice <- SplitBy.Choice[!unname(apply(data$obj@meta.data[,SplitBy.Choice, drop = FALSE],2,function(x)any(is.na(x))))] #split choice列中不可以有NA值
    names(SplitBy.Choice) <- SplitBy.Choice
    selectInput("Split","Split by:", choices = c("None" = "None", SplitBy.Choice))
  })


  # # check if the split is set appropriate (not None or NA) 如何让此处的代码优雅一些！，然后如何优化这个R包！ 2024.05.11
  Split.Revised <- reactive({
    req(input$Split) # split值出现后，才会执行的代码
    # Revise the Split choice
    if(is.na(input$Split)) {
      return(NULL)
    }else if(input$Split == "None"){
      return(NULL)
    }else{ return(input$Split) }
  })

  # define all available features used for featureplot
  GeneLibrary <- reactive({
    c(rownames(data$obj),c("nCount_RNA", "nFeature_RNA","percent.mt","log10GenesPerUMI"))
  })

  # Check the input gene
  Gene.Revised <- reactive({
    CheckGene(InputGene = ifelse(is.na(input$GeneSymbol), NA, input$GeneSymbol), GeneLibrary = GeneLibrary())
  })


  # Scatter Plot
  PlotMain <- reactive({
    # Dimplot, Attention: wrong input gene will also output the dimplot!
    if(is.na(Gene.Revised())) {
      if (is.null(Split.Revised())) { # not splited
        Seurat::DimPlot(data$obj, reduction = input$DimensionReduction, label = input$ShowLabel, pt.size = input$PointSize, label.size = input$LabelSize, group.by = input$ClusterResolution)
      }else{ # splited
        plot_numbers <- length(levels(data$obj@meta.data[,Split.Revised()]))
        Seurat::DimPlot(data$obj, reduction = input$DimensionReduction, label = input$ShowLabel, pt.size = input$PointSize, label.size = input$LabelSize, group.by = input$ClusterResolution,
                split.by = Split.Revised(), ncol = ceiling(sqrt(plot_numbers)))
      }
    }else { # Feature plot
      if (is.null(Split.Revised())) { # not splited
        # Gene.Revised()[1]: could set only use the first gene (input DL, get Dl and dl)
        Seurat::FeaturePlot(data$obj, features = Gene.Revised(), pt.size = input$PointSize, label.size = input$LabelSize, reduction = input$DimensionReduction, cols = c("gray","red"))
      } else { # splited
        plot_numbers <- length(levels(data$obj@meta.data[,Split.Revised()]))
        Seurat::FeaturePlot(data$obj, features = Gene.Revised(), pt.size = input$PointSize, label.size = input$LabelSize, reduction = input$DimensionReduction, cols = c("gray","red"),
                    split.by = Split.Revised()) +
          patchwork::plot_layout(ncol = ceiling(sqrt(plot_numbers)),nrow = ceiling(plot_numbers/ceiling(sqrt(plot_numbers))))
      }
    }
  })

  # MainPlot
  output$MainPlot <- renderPlot({PlotMain()}, height = function(){session$clientData$output_MainPlot_width * input$MainPlotHWRatio}) # box plot: height = width default


}
