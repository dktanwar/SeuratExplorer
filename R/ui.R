# ui.R
# R shiny UI for SeuratExplorer

#' UI for shiny App interface
#' @import shiny
#' @export
ui <-  function(){
  # Header ----
  header = shinydashboard::dashboardHeader(title = "SeuratExplorer")

  # Sidebar ----
  sidebar = shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      menuItem("Dataset", tabName = "dataset", icon = icon("database")),
      conditionalPanel(
        condition = "output.file_loaded",
        shinydashboard::sidebarMenu(menuItem("Explorer", tabName = "explorer", icon = icon("dashboard")))
        )
     )
  )

  # BODY ----
  tab_list = list()

  tab_list[["dataset"]] = tabItem(tabName = "dataset",
                                  # 上传文件
                                  box(fileInput("dataset_file", "Choose Seurat .rds file:", accept = '.rds'), status = "primary", width = 12),
                                  conditionalPanel(
                                    condition = "output.file_loaded",
                                    box(title = "Cell Meta Info", collapsible = TRUE,
                                      shinycssloaders::withSpinner(DT::dataTableOutput('dataset_meta')), status = "primary", width = 12) # 计划：改成在上传完数据后，再显示这个box！
                                    )
                                  )

  tab_list[["explorer"]] = tabItem(tabName = "explorer",
                                  # 数据展示
                                  fluidRow(
                                    box(title = "Main Plot",
                                        shinycssloaders::withSpinner(plotOutput("MainPlot",height = "auto")), # Add a spinner that shows when an output is recalculating
                                        width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                    box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                        # proxy.height: spinner height
                                        shinycssloaders::withSpinner(uiOutput("Reductions.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("ClusterResolution.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("Split.UI"), proxy.height = "10px"),
                                        textInput("GeneSymbol", "Gene Symbol:", value = ""),
                                        # # helpText(strong(paste("Support: ", paste(Genes.qc, collapse = ", "), ".",sep = "")),style = "font-size:5px;"),
                                        # adjust the Ratio of width and height of plot.
                                        sliderInput("MainPlotHWRatio", label = "Adjust H/W Ratio of Main Plot", min = 0.1, max = 2, value = 0.9),
                                        checkboxInput("ShowLabel",label = "Show cluster label", TRUE),
                                        sliderInput("LabelSize", label = "Label Size:", min = 0, max = 10, value = 7),
                                        sliderInput("PointSize", label = "Point Size", min = 0.001, max = 2, value = 0.8),
                                        # downloadButton('DownloadPlot', label = 'Download Scater Plot - Not Recommmended'),
                                        # textInput("Group1","Group 1:", value = "1"),
                                        # textInput("Group2","Group 2:", value = "2"),
                                        # actionButton("CalculateDEG", "Calculate DEGs!"),
                                        # textInput("Cluster", "Cluster:", value = "0"),
                                        # actionButton("CalculateDEGWithin", "Calculate DEGs Within A Cluster"),
                                        # shinycssloaders::withSpinner(uiOutput("Vlnplot.Cluster.Selection.UI"), proxy.height = "10px"),
                                        # # adjust the Ratio of width and height of scater plot
                                        # sliderInput("ClusterAnnotationPlotHWRatio", label = "Cluster Annotation Plot Height/Width Ratio", min = 0.1, max = 2, value = 1)
                                      )
                                    )
                                  )
  body = shinydashboard::dashboardBody(
    div(class= "tab-content", tab_list)
  )

  # 整合到一起
  ui_out = shinydashboard::dashboardPage(header, sidebar, body)
  return(ui_out)
}



