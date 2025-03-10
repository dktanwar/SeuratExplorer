# launch.R
## used to launch the shiny app in a web browser.


#' Launch shiny app
#'
#' @param verbose for debug use
#'
#' @import shiny
#' @import shinydashboard
#' @return In-browser Shiny Application launch
#' @examples
#' if(interactive()){launchSeuratExplorer()}
#' @export
launchSeuratExplorer <- function(verbose=FALSE){
  options(SeuratExplorerVerbose = verbose)
  app = shinyApp(ui, server)
  runApp(app, launch.browser = TRUE)
}
