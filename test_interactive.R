# Test script for SeuratExplorer interactive plots
# This script tests the interactive plotting functions

library(SeuratExplorer)
library(Seurat)
library(plotly)

# Test the interactive plotting functions with a mock ggplot
test_interactive_functions <- function() {
  cat("Testing SeuratExplorer interactive plotting functions...\n")
  
  # Test if the functions exist and are callable
  functions_to_test <- c(
    "interactive_dimplot",
    "interactive_featureplot", 
    "interactive_vlnplot",
    "interactive_dotplot",
    "interactive_empty_plot"
  )
  
  for (func_name in functions_to_test) {
    if (exists(func_name)) {
      cat("✓", func_name, "function exists\n")
    } else {
      cat("✗", func_name, "function NOT found\n")
    }
  }
  
  # Test empty plot function
  tryCatch({
    empty_plot <- interactive_empty_plot("Test message")
    if (inherits(empty_plot, "plotly")) {
      cat("✓ interactive_empty_plot works correctly\n")
    } else {
      cat("✗ interactive_empty_plot does not return plotly object\n")
    }
  }, error = function(e) {
    cat("✗ interactive_empty_plot failed:", e$message, "\n")
  })
  
  cat("\nInteractive plotting functions test completed!\n")
  cat("To fully test the interactive features:\n")
  cat("1. Load a Seurat object\n")
  cat("2. Launch SeuratExplorer()\n") 
  cat("3. Switch plot modes from 'Static' to 'Interactive'\n")
  cat("4. Verify plots are interactive with hover, zoom, and pan\n")
}

# Run the test
test_interactive_functions()

# Example of how to test with actual data (uncomment when you have a Seurat object)
# seu_obj <- readRDS("path/to/your/seurat_object.rds")
# SeuratExplorer(seu_obj)