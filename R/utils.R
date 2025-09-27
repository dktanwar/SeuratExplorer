# Utility functions for SeuratExplorer

#' Render a plot and save it to a temporary file
#'
#' @param output The output object from the shiny server.
#' @param session The session object from the shiny server.
#' @param plot_type The type of the plot (e.g., "dimplot", "featureplot").
#' @param plot_width A reactive expression for the width of the plot.
#'
#' @param plot_height_ratio A reactive expression for the height/width ratio of the plot.
#' @param temp_dir The temporary directory to save the plot to.
#' @param plot_obj A reactive expression for the plot object.
#'
#' @return A rendered plot.
render_plot <- function(output, session, plot_type, plot_width, plot_height_ratio, temp_dir, plot_obj) {
  output[[plot_type]] <- renderPlot({
    p <- plot_obj()
    ggplot2::ggsave(paste0(temp_dir, "/", plot_type, ".pdf"),
                    p,
                    width = plot_width() * 0.03,
                    height = plot_width() * plot_height_ratio() * 0.03,
                    units = "cm",
                    limitsize = FALSE)
    return(p)
  }, height = function() {
    session$clientData[[paste0("output_", plot_type, "_width")]] * plot_height_ratio()
  })
}

#' Create a download handler for a plot
#'
#' @param output The output object from the shiny server.
#' @param plot_type The type of the plot (e.g., "dimplot", "featureplot").
#' @param temp_dir The temporary directory where the plot is saved.
#'
#' @return A download handler.
download_plot <- function(output, plot_type, temp_dir) {
  output[[paste0("download", plot_type)]] <- downloadHandler(
    filename = function() {
      paste0(plot_type, ".pdf")
    },
    content = function(file) {
      file.copy(paste0(temp_dir, "/", plot_type, ".pdf"), file, overwrite = TRUE)
    }
  )
}
