# ui.R
# R shiny UI for SeuratExplorer

#' UI for shiny App interface
#' @import shiny
#' @import shinydashboard
#' @export
ui <-  function(){
  requireNamespace("shinydashboard")
  requireNamespace("shinyWidgets")

  # Header ----
  header = dashboardHeader(title = "SeuratExplorer")

  # Sidebar ----
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Dataset", tabName = "dataset", icon = icon("database")),
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
     )
  )

  # BODY ----
  tab_list = list()

  tab_list[["dataset"]] = tabItem(tabName = "dataset",
                                  # 上传文件
                                  box(fileInput("dataset_file", "Choose A Seurat .rds file:", accept = '.rds'), status = "primary", width = 12),
                                  conditionalPanel(
                                    condition = "output.file_loaded",
                                    box(title = "Cell Meta Info", collapsible = TRUE,
                                        shinycssloaders::withSpinner(DT::dataTableOutput('dataset_meta')), status = "primary", width = 12)
                                    )
                                  )

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

                                        # downloadButton('DownloadPlot', label = 'Download Scater Plot - Not Recommmended'),
                                        # textInput("Group1","Group 1:", value = "1"),
                                        # textInput("Group2","Group 2:", value = "2"),
                                        # actionButton("CalculateDEG", "Calculate DEGs!"),
                                        # textInput("Cluster", "Cluster:", value = "0"),
                                        # actionButton("CalculateDEGWithin", "Calculate DEGs Within A Cluster"),
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

                                        # downloadButton('DownloadPlot', label = 'Download Scater Plot - Not Recommmended'),
                                        # textInput("Group1","Group 1:", value = "1"),
                                        # textInput("Group2","Group 2:", value = "2"),
                                        # actionButton("CalculateDEG", "Calculate DEGs!"),
                                        # textInput("Cluster", "Cluster:", value = "0"),
                                        # actionButton("CalculateDEGWithin", "Calculate DEGs Within A Cluster"),
                                    )
                                  )
  )

  body = dashboardBody(
    div(class= "tab-content", tab_list)
  )

  # 整合到一起
  ui_out = dashboardPage(header, sidebar, body)
  return(ui_out)
}



