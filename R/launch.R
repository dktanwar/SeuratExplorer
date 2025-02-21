# launch.R
# used to launch the shiny app in a web browser.


#' Launch shiny app
#'
#' @import shiny
#' @import shinydashboard
#' @return In-browser Shiny Application launch
#' @examples
#' # launchSeuratExplorer()
#' @export
launchSeuratExplorer <- function(){
  app = shinyApp(ui, server)
  runApp(app, launch.browser = TRUE)
}
