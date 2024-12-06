# library(Seurat)
# cds <- readRDS("../SeuratExplorerServer/inst/extdata/demo/fly/Rds-file/G101_PC20res04.rds")
# features <- c("AstA","AstC")
# features %in% rownames(cds)
#
# SeuratObj <- cds
#
#
# requireNamespace("dplyr")
# requireNamespace("Seurat")
# # 如果是assay5类型
# if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
#   SeuratObj <- JoinLayers(SeuratObj)
#   normalized.expr <- as.matrix(SeuratObj@assays$RNA@layers$data)
# } else { # 如果是assay类型
#   normalized.expr <- as.matrix(SeuratObj@assays$RNA$data)
# }
#
# group.by = "RandomGroup"
# all.cell.types <- unique(SeuratObj@meta.data[,group.by])
#
# for (celltype in all.cell.types) {
#   cells.sub <- colnames(SeuratObj)[as.character(SeuratObj@meta.data[,group.by]) == celltype]
#   normalized.expr.sub <- normalized.expr[,cells.sub]
#   mean.expr <- apply(normalized.expr[features,,drop = FALSE], 1, mean)
#   median.expr <- apply(normalized.expr[features,,drop = FALSE], 1, median)
#   pct <- apply(normalized.expr[features,,drop = FALSE] > 0, 1, mean)
#   single.res <- data.frame(Gene = features, Expr.mean = mean.expr, Expr.median = median.expr, PCT = pct)
#   single.res$CellType <- celltype
#   single.res$TotalCells <- ncol(normalized.expr.sub)
#   single.res <- single.res[,c("CellType", "TotalCells","Gene", "PCT", "Expr.mean", "Expr.median")]
#   rownames(single.res) <- NULL
#   if (celltype == all.cell.types[1]) {
#     res <- single.res
#   }else{
#     res <- rbind(res, single.res)
#   }
# }
#
#
