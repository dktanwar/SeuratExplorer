# library(Seurat)
# cds <- readRDS("../SeuratExplorerServer/inst/extdata/demo/fly/Rds-file/G101_PC20res04.rds")
#
# SeuratObj <- cds
# if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
#   SeuratObj <- JoinLayers(SeuratObj)
#   normalized.expr <- as.matrix(SeuratObj@assays$RNA@layers$data)
# } else { # 如果是assay类型
#   normalized.expr <- as.matrix(SeuratObj@assays$RNA$data)
# }
#
# x <- normalized.expr[feature, ,drop = FALSE]
# y <- normalized.expr[rownames(normalized.expr)[rownames(normalized.expr) != feature], ]
# # 过滤细胞，表达值的和要大于细胞数目的1/10，设置比较随意！
# y <- y[apply(y, 1, sum) > 0.1 * ncol(SeuratObj),]
#
# cor.res = cor(x = t(as.matrix(x)), y =  t(as.matrix(y)), method = method)
# cor.res <- reshape2::melt(cor.res)
# colnames(cor.res) <- c("GeneA","GeneB","correlation")
# # cor.res <- cor.res[abs(cor.res$correlation) > 0.6,]
# cor.res <- cor.res[order(abs(cor.res$correlation),decreasing = TRUE),]
# rownames(cor.res) <- NULL
# return(cor.res)




