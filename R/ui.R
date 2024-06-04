# ui.R
# R shiny UI for SeuratExplorer

#' sidebar - ui functions for Seurat explorer functions
#' @export
explorer_sidebar_ui <- function(){
  requireNamespace("shinydashboard")
  requireNamespace("shinyWidgets")
  conditionalPanel(
    condition = "output.file_loaded",
    sidebarMenu(menuItem(text = "Explorer", tabName = "explorer", icon = icon("dashboard"), startExpanded = TRUE,
                         menuSubItem(text = "DimPlot", tabName = "dimplot", icon = shiny::icon("angle-double-right")),
                         menuSubItem(text = "FeaturePlot", tabName = "featureplot", icon = shiny::icon("angle-double-right")),
                         menuSubItem(text = "VlnPlot", tabName = "vlnplot", icon = shiny::icon("angle-double-right")),
                         menuSubItem(text = "DotPlot", tabName = "dotplot", icon = shiny::icon("angle-double-right")),
                         menuSubItem(text = "HeatmapPlot", tabName = "heatmap", icon = shiny::icon("angle-double-right")),
                         menuSubItem(text = "RidgePlot", tabName = "ridgeplot", icon = shiny::icon("angle-double-right")),
                         menuSubItem(text = "DEG Analysis", tabName = "degs", icon = shiny::icon("angle-double-right"))))
  )
}

#' body - ui functions for Seurat explorer functions
#' @param tab_list settings for each seurat function
#' @export
explorer_body_ui <- function(tab_list){
  requireNamespace("shinydashboard")
  requireNamespace("shinyWidgets")
  tab_list[["dimplot"]] = tabItem(tabName = "dimplot",
                                  fluidRow(
                                    box(title = "Dimension Reduction Plot",
                                        shinycssloaders::withSpinner(plotOutput("dimplot",height = "auto")), # Add a spinner that shows when an output is recalculating
                                        width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                    box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                        shinycssloaders::withSpinner(uiOutput("DimReductions.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("DimClusterResolution.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("DimSplit.UI"), proxy.height = "10px"),
                                        sliderInput("DimPlotHWRatio", label = "Adjust H/W Ratio of DimPlot", min = 0.1, max = 2, value = 0.9),# adjust the Ratio of width and height of plot.
                                        checkboxInput("DimShowLabel",label = "Show cluster label", TRUE),
                                        sliderInput("DimLabelSize", label = "Label Size:", min = 0, max = 10, value = 7),
                                        sliderInput("DimPointSize", label = "Point Size", min = 0.001, max = 2, value = 0.8)
                                    )
                                  )
  )
  tab_list[["featureplot"]] = tabItem(tabName = "featureplot",
                                      fluidRow(
                                        box(title = "Visualize 'features' on a dimensional reduction plot",
                                            shinycssloaders::withSpinner(plotOutput("featureplot",height = "auto")), # Add a spinner that shows when an output is recalculating
                                            width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                        box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                            textInput("FeatureGeneSymbol", "Gene Symbol:", value = ""),
                                            shinycssloaders::withSpinner(uiOutput("Featurehints.UI"), proxy.height = "10px"),
                                            shinycssloaders::withSpinner(uiOutput("FeatureReductions.UI"), proxy.height = "10px"),
                                            # shinycssloaders::withSpinner(uiOutput("FeatureClusterResolution.UI"), proxy.height = "10px"),
                                            shinycssloaders::withSpinner(uiOutput("FeatureSplit.UI"), proxy.height = "10px"),
                                            # 拾色器参考： https://daattali.com/shiny/colourInput/
                                            colourpicker::colourInput("FeaturePlotLowestExprColor", "Pick Color for lowest expression:", "#E5E5E5", palette = "limited"),
                                            colourpicker::colourInput("FeaturePlotHighestExprColor", "Pick Color for highest expression:", "#FF0000",palette = "limited"),
                                            sliderInput("FeaturePlotHWRatio", label = "Adjust Height/Width Ratio:", min = 0.1, max = 2, value = 0.9), # adjust the Ratio of width and height of plot.
                                            # checkboxInput("FeatureShowLabel",label = "Show cluster label", TRUE),
                                            # sliderInput("FeatureLabelSize", label = "Label Size:", min = 0, max = 10, value = 7),
                                            sliderInput("FeaturePointSize", label = "Point Size:", min = 0.001, max = 2, value = 0.8),
                                        )
                                      )
  )
  tab_list[["vlnplot"]] = tabItem(tabName = "vlnplot",
                                  fluidRow(
                                    box(title = "Single cell violin plot",
                                        shinycssloaders::withSpinner(plotOutput("vlnplot",height = "auto")), # Add a spinner that shows when an output is recalculating
                                        width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                    box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                        textInput("VlnGeneSymbol", "Gene Symbols:", value = ""),
                                        shinycssloaders::withSpinner(uiOutput("Vlnhints.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("VlnClusterResolution.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("VlnSplitBy.UI"), proxy.height = "10px"),
                                        conditionalPanel(
                                          condition = "output.Vlnplot_splitoption_twolevels",
                                          checkboxInput("VlnSplitPlot",label = "Split Plot", FALSE)
                                        ),
                                        conditionalPanel(
                                          condition = "output.Vlnplot_multiple_genes",
                                          checkboxInput("VlnStackPlot",label = "Stack Plot", FALSE)
                                        ),
                                        conditionalPanel(
                                          condition = "output.Vlnplot_StackPlot",
                                          checkboxInput("VlnFlipPlot",label = "Flip Plot", FALSE)
                                        ),
                                        conditionalPanel(
                                          condition = "output.Vlnplot_StackPlot && input.VlnSplitBy == 'None'", # 仅对于split为NULL时有效。
                                          selectInput("VlnFillBy","Color By:", choices = c(Feature = "feature", Ident = "ident"))
                                        ),
                                        sliderInput("VlnPointSize", label = "Point Size:", min = 0, max = 4, value = 0),
                                        sliderInput("VlnPointAlpha", label = "Point Alpha:", min = 0, max = 1, value = 1),
                                        sliderInput("VlnXlabelSize", label = "x Axis Label Size:", min = 0, max = 20, value = 14),
                                        sliderInput("VlnYlabelSize", label = "Y Axis Label Size:", min = 0, max = 20, value = 10),
                                        sliderInput("VlnPlotHWRatio", label = "Adjust Height/Width Ratio:", min = 0.1, max = 2, value = 0.9)# adjust the Ratio of width and height of plot.
                                    )
                                  )
  )
  tab_list[["dotplot"]] = tabItem(tabName = "dotplot",
                                  fluidRow(
                                    box(title = "Dot plot visualization",
                                        shinycssloaders::withSpinner(plotOutput("dotplot",height = "auto")), # Add a spinner that shows when an output is recalculating
                                        width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                    box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                        textInput("DotGeneSymbol", "Gene Symbols:", value = ""),
                                        shinycssloaders::withSpinner(uiOutput("Dothints.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("DotClusterResolution.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("DotIdentsSelected.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("DotSplitBy.UI"), proxy.height = "10px"),
                                        checkboxInput("DotClusterIdents",label = "Cluster Idents", FALSE),
                                        checkboxInput("DotRotateAxis",label = "Rotate Axis", FALSE),
                                        checkboxInput("DotFlipCoordinate",label = "Flip Coordinate", FALSE),
                                        conditionalPanel(
                                          condition = "output.DotPlot_Split_isNone",
                                          colourpicker::colourInput("DotPlotLowestExprColor", "Pick Color for lowest expression:", "#E5E5E5", palette = "limited"),
                                          colourpicker::colourInput("DotPlotHighestExprColor", "Pick Color for highest expression:", "#0000FF",palette = "limited"),
                                        ),
                                        sliderInput("DotDotScale", label = "Dot Scale:", min = 1, max = 12, value = 6),
                                        sliderInput("DotXlabelSize", label = "x Axis Label Size:", min = 0, max = 20, value = 14),
                                        sliderInput("DotYlabelSize", label = "Y Axis Label Size:", min = 0, max = 20, value = 10),
                                        sliderInput("DotPlotHWRatio", label = "Adjust Height/Width Ratio:", min = 0.1, max = 2, value = 0.9) # adjust the Ratio of width and height of plot.
                                    )
                                  )
  )
  tab_list[["heatmap"]] = tabItem(tabName = "heatmap",
                                  fluidRow(
                                    box(title = "Feature expression heatmap",
                                        shinycssloaders::withSpinner(plotOutput("heatmap",height = "auto")), # Add a spinner that shows when an output is recalculating
                                        width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                    box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                        textInput("HeatmapGeneSymbol", "Gene Symbols:", value = ""),
                                        shinycssloaders::withSpinner(uiOutput("Heatmaphints.UI"), proxy.height = "10px"),
                                        shinycssloaders::withSpinner(uiOutput("HeatmapClusterResolution.UI"), proxy.height = "10px"),
                                        sliderInput("HeatmapTextSize", label = "Cluster Text Size:", min = 1, max = 12, value = 6),
                                        sliderInput("HeatmapTextHjust", label = "Cluster Text Hjust:", min = -10, max = 20, value = 0),
                                        sliderInput("HeatmapTextVjust", label = "Cluster Text Vjust:", min = -0.55, max = 0.55, value = 0),
                                        sliderInput("HeatmapTextRatateAngle", label = "Cluster Text Rotate Angle:", min = -90, max = 90, value = 0),
                                        sliderInput("HeatmapGroupBarHeight", label = "Cluster Group Bar Height:", min = 0, max = 0.1, value = 0.05),
                                        sliderInput("HeatmapLineWidth", label = "Line Width:", min = 1, max = 10, value = 1),
                                        sliderInput("HeatmapFeatureTextSize", label = "Feature Text Size:", min = 0, max = 20, value = 10),
                                        sliderInput("HeatmapPlotHWRatio", label = "Adjust Height/Width Ratio:", min = 0.1, max = 2, value = 0.9) # adjust the Ratio of width and height of plot.
                                    )
                                  )
  )
  tab_list[["ridgeplot"]] = tabItem(tabName = "ridgeplot",
                                    fluidRow(
                                      box(title = "Single cell ridge plot",
                                          shinycssloaders::withSpinner(plotOutput("ridgeplot",height = "auto")), # Add a spinner that shows when an output is recalculating
                                          width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                                      box(title = "Settings", solidHeader = TRUE, status = "primary", width = 3,
                                          textInput("RidgeplotGeneSymbol", "Gene Symbols:", value = ""),
                                          shinycssloaders::withSpinner(uiOutput("Ridgeplothints.UI"), proxy.height = "10px"),
                                          shinycssloaders::withSpinner(uiOutput("RidgeplotClusterResolution.UI"), proxy.height = "10px"),
                                          conditionalPanel(
                                            condition = "output.Ridgeplot_stack_NotSelected",
                                            sliderInput("RidgeplotNumberOfColumns", label = "Number of columns:", min = 1, max = 10, value = 1),
                                          ),
                                          conditionalPanel(
                                            condition = "output.Ridgeplot_stack_show",
                                            checkboxInput("RidgeplotStackPlot",label = "Stack Plot", FALSE),
                                          ),
                                          conditionalPanel(
                                            condition = "input.RidgeplotStackPlot",
                                            selectInput("RidgeplotFillBy","Color By:", choices = c(Feature = "feature", Ident = "ident"))
                                          ),
                                          sliderInput("RidgeplotXlabelSize", label = "x Axis Label Size:", min = 0, max = 20, value = 14),
                                          sliderInput("RidgeplotYlabelSize", label = "Y Axis Label Size:", min = 0, max = 20, value = 10),
                                          sliderInput("RidgeplotHWRatio", label = "Adjust Height/Width Ratio:", min = 0.1, max = 2, value = 0.9) # adjust the Ratio of width and height of plot.
                                      )
                                    )
  )
  tab_list[["degs"]] = tabItem(tabName = "degs",
                               fluidRow(
                                 box(textOutput("degs_warning"), title = "WARNING：", background = "orange", width = 12),
                                 tags$style(".nav-tabs {background: #f4f4f4;}
                                 .nav-tabs-custom .nav-tabs li.active:hover a, .nav-tabs-custom .nav-tabs li.active a {background-color: #fff;
                                 border-color: #fff;
                                 }
                                 .nav-tabs-custom .nav-tabs li.active {border-top-color:
                                 #314a6d;
                                 }"), # refer to: https://stackoverflow.com/questions/45247290/shiny-dashboard-tabbox-tabpanel-css
                                 tabBox(
                                   title = "Find Markers or DEGs",
                                   id = "tabset_degs", width = 12, # height = "250px",
                                   tabPanel("ClusterMarkers", strong(h3("Find Markers for All Clusters")),
                                            shinycssloaders::withSpinner(uiOutput("ClusterMarkersClusterResolution.UI"), proxy.height = "10px"),
                                            actionButton("DEGsClusterMarkersAnalysis", "Analyze", icon = icon("magnifying-glass-chart"), class = "btn-primary")),
                                   # 功能冗余
                                   # tabPanel("InterClusterDEGs", strong(h3("Find DEGs Between two Cluster Groups")),
                                   #          shinycssloaders::withSpinner(uiOutput("InterClusterDEGsClusterResolution.UI"), proxy.height = "10px"),
                                   #          shinycssloaders::withSpinner(uiOutput("InterClusterDEGsGroupCase.UI"), proxy.height = "10px"),
                                   #          shinycssloaders::withSpinner(uiOutput("InterClusterDEGsGroupControl.UI"), proxy.height = "10px"),
                                   #          actionButton("InterClusterDEGsAnalysis", "Analyze")),
                                   tabPanel("IntraClusterDEGs", strong(h3("Find DEGs by Customized Groups Within Selected Clusters")),
                                            shinycssloaders::withSpinner(uiOutput("IntraClusterDEGsCustomizedGroups.UI"), proxy.height = "10px"),
                                            shinycssloaders::withSpinner(uiOutput("IntraClusterDEGsCustomizedGroupsCase.UI"), proxy.height = "10px"),
                                            shinycssloaders::withSpinner(uiOutput("IntraClusterDEGsCustomizedGroupsControl.UI"), proxy.height = "10px"),
                                            tags$hr(style="border: none; border-top: 1px dashed #ccc;"),
                                            strong(h3("Modify the following parameters If you want to subset cells, otherwise ignore it.")),
                                            shinycssloaders::withSpinner(uiOutput("IntraClusterDEGsSubsetCells.UI"), proxy.height = "10px"),
                                            shinycssloaders::withSpinner(uiOutput("IntraClusterDEGsSubsetCellsSelectedClusters.UI"), proxy.height = "10px"),
                                            tags$hr(style="border: none; border-top: 1px dashed #ccc;"),
                                            actionButton("IntraClusterDEGssAnalysis", "Analyze", icon = icon("magnifying-glass-chart"), class = "btn-primary")),
                                   tabPanel("Parameters", strong(h3("Set Parameters")),
                                            sliderInput("logfcthreshold", label = "Logfc Threshold:", min = 0, max = 1, value = 0.1),
                                            selectInput("testuse","Test use:", choices = c(wilcox = "wilcox", wilcox_limma = "wilcox_limma",
                                                                                           T_test = "t", negbinom = "negbinom", poisson = "poisson",
                                                                                           LR = "LR", MAST = "MAST", DESeq2 = "DESeq2")),
                                            sliderInput("minpct", label = "Minimum Expression Percentage:", min = 0, max = 1, value = 0.01),
                                            sliderInput("mindiffpct", label = "Minimum Expression Percentage Difference:", min = 0, max = 1, value = 0),
                                            actionButton("SetDefault", "Set to Default", icon = icon("save"), class = "btn-primary"))

                                 ),
                                 conditionalPanel(
                                   condition = "output.DEGs_ready",
                                   box(title = "Analysis Results:", collapsible = TRUE, width = 12,
                                       shinycssloaders::withSpinner(DT::dataTableOutput('dataset_degs')))
                                 )
                               )
  )
  return(tab_list)
}



#' UI for shiny App interface
#' @import shiny
#' @import shinydashboard shinyWidgets
#' @export
ui <-  function(){
  requireNamespace("shinydashboard")
  requireNamespace("shinyWidgets")

  # notificationItem 默认函数无法在新页面打开链接; refer to: https://forum.posit.co/t/shinydashboard-notification-item-with-link-in-new-tab/37580/2
  notificationItemWithAttr <- function(text, icon = shiny::icon("warning"), status = "success", href = NULL, ...) {
    if (is.null(href))
      href <- "#"
    icon <- tagAppendAttributes(icon, class = paste0("text-",
                                                     status))
    tags$li(a(href = href, icon, text, ...))
  }

  # Header ----
  header = dashboardHeader(title = "SeuratExplorer",
                           dropdownMenu(type = "notifications", icon = icon("github"), headerText = "R packages on Github:",
                                        notificationItemWithAttr(icon = icon("github"), status = "info", text = "SeuratExplorer", href = "https://github.com/fentouxungui/SeuratExplorer", target = "_blank"),
                                        notificationItemWithAttr(icon = icon("github"), status = "info", text = "SeuratExplorerServer", href = "https://github.com/fentouxungui/SeuratExplorerServer", target = "_blank")))

  # Sidebar ----
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Dataset", tabName = "dataset", icon = icon("database")),
      explorer_sidebar_ui()
     )
  )

  # BODY ----
  tab_list = list()

  tab_list[["dataset"]] = tabItem(tabName = "dataset",
                                  fluidRow(
                                    # 上传文件
                                    box(status = "primary", title = "Upload Data", width = 12, collapsible = TRUE, solidHeader = TRUE,
                                        fileInput("dataset_file", "Choose A Seurat .rds file:", accept = '.rds')),
                                    conditionalPanel(
                                      condition = "output.file_loaded",
                                      box(title = "Metadata of Cells", collapsible = TRUE, width = 12,
                                          shinycssloaders::withSpinner(DT::dataTableOutput('dataset_meta')))
                                    ))
                                  )

  tab_list <- explorer_body_ui(tab_list = tab_list)

  body = dashboardBody(
    div(class= "tab-content", tab_list),
    tags$script(HTML(
      "document.querySelector('body > div.wrapper > header > nav > div > ul > li > a > span').style.visibility = 'hidden';"
    )) # 不显示dropdownMenu中notification的数目， refer to:https://stackoverflow.com/questions/65915414/alter-dropdown-menu-in-shiny
  )

  # 整合到一起
  ui_out = dashboardPage(header, sidebar, body)
  return(ui_out)
}



