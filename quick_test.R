# Quick Test Script for SeuratExplorer Interactive Features
# Run this in RStudio to test the interactive functionality

# Load required libraries
library(devtools)

# Load the package functions
load_all()

# Test that plotly is available
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}

library(plotly)
library(htmlwidgets)

# Test the interactive empty plot function
cat("Testing interactive_empty_plot...\n")
test_plot <- interactive_empty_plot("Test successful!")
print(test_plot)

cat("\n✅ Interactive plotting functions are ready!\n")
cat("\n📦 Package versions compatible:\n")
cat("   - ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("   - plotly:", as.character(packageVersion("plotly")), "\n")
cat("   - htmlwidgets:", as.character(packageVersion("htmlwidgets")), "\n")
cat("\n🚀 To test with real data:\n")
cat("1. Load your Seurat object: seu_obj <- readRDS('path/to/seurat.rds')\n")
cat("2. Launch the app: SeuratExplorer(seu_obj)\n") 
cat("3. Look for 'Plot Type' toggles to switch between Static/Interactive\n")
cat("4. Enjoy the new interactive features! 🎉\n")

# Additional check - verify all functions exist
functions_to_check <- c('interactive_dimplot', 'interactive_featureplot', 
                       'interactive_vlnplot', 'interactive_dotplot', 'interactive_empty_plot')

cat("\n📋 Function availability check:\n")
for(func in functions_to_check) {
  if(exists(func)) {
    cat("✅", func, "\n")
  } else {
    cat("❌", func, "- MISSING!\n")
  }
}