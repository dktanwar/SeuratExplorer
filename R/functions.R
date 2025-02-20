prepare_seurat_object <- function(obj){
  requireNamespace("Seurat")
  # 将meta.data中的部分非factor类型的列，转为factor类型
  # 如果unique数目少于细胞总数的1/20，并且不大于50种，并且数据类型为chr或num类型，会强制转为factor类型。
  # 可能的问题： unique_max_percent = 0.05可能不适合只有100个细胞但由大于5群的的数据
  obj@meta.data <- modify_columns_types(df = obj@meta.data, types_to_check = c("numeric", "character"), unique_max_counts = 50, unique_max_percent = 0.05)
  # for splited object, join layers
  if (sum(grepl("^counts",Layers(object = obj))) > 1 | sum(grepl("^data",Layers(object = obj))) > 1) {
    obj <- SeuratObject::JoinLayers(object = obj)
  }
  message("SeuratExplorer: prepare_seurat_object runs successfully!")
  return(obj)
}

# 把符合条件的非因子列，转为因子类型, 并且把可能是数字的字符串转为数字
modify_columns_types <- function(df, types_to_check = c("numeric", "character"), unique_max_counts = 50, unique_max_percent = 0.05){
  candidates.types.logic <- sapply(df, class) %in% types_to_check
  cutoff <- min(unique_max_counts, round(nrow(df) * unique_max_percent))
  candidates.unique.logic <- sapply(df, FUN = function(x)length(unique(x))) <= cutoff
  candidates.logic <- candidates.types.logic & candidates.unique.logic & !sapply(df, is.factor)
  # before trans to factor, trans numeric character vector to numeric vector: c('1','2','1') to c(1,2,1)
  char2numeric_columns.logic <- suppressWarnings(unlist(lapply(df[candidates.logic], function(x)any(!is.na(as.numeric(unique(x)))))))
  char2numeric_columns <- names(char2numeric_columns.logic)[char2numeric_columns.logic]
  df[char2numeric_columns] <- lapply(df[char2numeric_columns], as.numeric)
  # finally trans all char and numeric to factor
  df[candidates.logic] <- lapply(df[candidates.logic], as.factor)
  message("SeuratExplorer: modify_columns_types runs successfully!")
  return(df)
}

# 通过关键字换取reduction options
prepare_reduction_options <- function(obj, keywords = c("umap","tsne")){
  requireNamespace("Seurat")
  reduction.choice <- grep(paste0(paste0("(", keywords,")"),collapse = "|"), Seurat::Reductions(obj), value = TRUE, ignore.case = TRUE)
  names(reduction.choice) <- toupper(reduction.choice)
  message("SeuratExplorer: prepare_reduction_options runs successfully!")
  return(reduction.choice)
}


# 将所有meta.data中所有类型为因子的列名作为cluster options
prepare_cluster_options <- function(df){
  cluster.options <- colnames(df)[sapply(df, is.factor)]
  names(cluster.options) <- cluster.options
  message("SeuratExplorer: prepare_cluster_options runs successfully!")
  return(cluster.options)
}

# 将所有meta.data中level数目少于max_level的因子列作为split options
prepare_split_options <- function(df, max.level = 4){
  cluster.options <- colnames(df)[sapply(df, is.factor)]
  leve.counts <- unname(sapply(df[cluster.options],FUN = function(x)length(levels(x))))
  split.options <- cluster.options[leve.counts <= max.level]
  names(split.options) <- split.options
  message("SeuratExplorer: prepare_split_options runs successfully!")
  return(split.options)
}

# 添加额外的来自meta data的列名为qc options
prepare_qc_options <- function(df, types = c("double","integer","numeric")){
  qc_options <- colnames(df)[sapply(df, class) %in% types]
  message("SeuratExplorer: prepare_qc_options runs successfully!")
  return(qc_options)
}


# Check the input gene, return the revised gene, which can be used for FeaturePlot, Vlnplot ect.
CheckGene <- function(InputGene, GeneLibrary){
  InputGenes <- trimws(unlist(strsplit(InputGene,split = "\n")))
  InputGenes <- InputGenes[InputGenes != ""]
  revised.genes <- sapply(InputGenes, FUN = function(x)ReviseGene(x, GeneLibrary = GeneLibrary))
  revised.genes <- unique(unname(revised.genes[!is.na(revised.genes)]))
  message("SeuratExplorer: CheckGene runs successfully!")
  ifelse(length(revised.genes) == 0, yes = return(NA), no = return(revised.genes))
}

ReviseGene <- function(Agene, GeneLibrary){
  if (Agene %in% GeneLibrary) { # when input gene is absolutely right
    return(Agene)
  }else if(tolower(Agene) %in% tolower(GeneLibrary)){ # when not case sensitive
    return(GeneLibrary[tolower(GeneLibrary) %in% tolower(Agene)][1]) # gene list length can > 1
  }else{ # when not match
    return(NA)
  }
}


# 检查差异分析时所需要的R包依赖
check_dependency <- function(test){
  if (test == "wilcox") { # 检查所需要的R包,暂时这样吧，其它几个test先不测试了
    if(!require(devtools, quietly = TRUE)){
      utils::install.packages("devtools")
    }
    if(!require(presto)){
      devtools::install_github('immunogenomics/presto', upgrade = "never")
    }
  }else if(test == "DESeq2"){
    if (!require("BiocManager", quietly = TRUE))
      utils::install.packages("BiocManager")
    if (!require("DESeq2", quietly = TRUE)) {
      BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
    }
  } else if(test == "MAST"){
    if (!require("BiocManager", quietly = TRUE))
      utils::install.packages("BiocManager")
    if (!require("MAST", quietly = TRUE)) {
      BiocManager::install("MAST", update = FALSE, ask = FALSE)
    }
  }
}

color_list <- list(stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                                "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                                "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                                "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D"),

                   stallion2 = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                                 "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                                 "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                                 "#6E4B9E","#0C727C", "#7E1416","#D8A767"),

                   calm = c("#7DD06F", "#844081","#688EC1", "#C17E73", "#484125",
                            "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
                            "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D",
                            "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736"),

                   kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020",
                             "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
                             "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800",
                             "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16"),

                   bear = c("#faa818", "#41a30d","#fbdf72", "#367d7d", "#d33502",
                            "#6ebcbc", "#37526d","#916848", "#f5b390", "#342739",
                            "#bed678","#a6d9ee", "#0d74b6",
                            "#60824f","#725ca5", "#e0598b"),

                   #15-colors
                   ironMan = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',
                               "#F88D51","#FAD510","#FFFF5F",'#88CFA4','#238B45',
                               "#02401B", "#0AD7D3","#046C9A", "#A2A475", 'grey35'),

                   circus = c("#D52126","#88CCEE", "#FEE52C", "#117733", "#CC61B0",
                              "#99C945", "#2F8AC4", "#332288","#E68316", "#661101",
                              "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"),

                   #12-colors
                   paired = c("#A6CDE2","#1E78B4","#74C476","#34A047","#F59899","#E11E26",
                              "#FCBF6E","#F47E1F","#CAB2D6","#6A3E98","#FAF39B","#B15928"),

                   #11-colors
                   grove = c("#1a1334","#01545a","#017351","#03c383","#aad962",
                             "#fbbf45","#ef6a32","#ed0345","#a12a5e","#710162","#3B9AB2"),

                   #7-colors
                   summerNight = c("#2a7185","#a64027","#fbdf72","#60824f","#9cdff0",
                                   "#022336","#725ca5"),

                   #5-colors
                   zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"),
                   darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"),
                   rushmore = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" , "#F2300F"),
                   captain = c("grey","#A1CDE1","#12477C","#EC9274","#67001E"),

                   #---------------------------------------------------------------
                   # Primarily Continuous Palettes
                   #---------------------------------------------------------------
                   #10-colors
                   horizon = c('#000075','#2E00FF', '#9408F7', '#C729D6', '#FA4AB5',
                               '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'),

                   #9-colors
                   horizonExtra =c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5",
                                   "#FD619D","#FF9965","#FFD32B","#FFFC5A"),

                   blueYellow = c("#352A86","#343DAE","#0262E0","#1389D2","#2DB7A3",
                                  "#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"),

                   sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                                  '#A2E700','#FFFF00','#FFD200','#FFA500'),

                   solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC',
                                  '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'),

                   whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                                   '#8c6bb1','#88419d','#810f7c','#4d004b'),

                   whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                                 '#3690c0','#0570b0','#045a8d','#023858'),


                   comet = c("#E6E7E8","#3A97FF","#8816A7","black"),

                   #7-colors
                   greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                                 '#0868ac','#084081'),

                   #6-colors
                   beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30",
                             "#F7962E","#FCEE2B"),

                   #5-colors
                   coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"),
                   fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),
                   greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"),
                   fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A"),
                   purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"))



color_choice_vector <- names(color_list)
names(color_choice_vector) <- paste(names(color_list), unlist(lapply(color_list,length)),sep = "-")
color_choice_vector <- color_choice_vector[names(color_choice_vector)[order(unlist(lapply(color_list,length)),decreasing = TRUE)]]

#' @title getColors
#'
#' @param color.platte predefined color list
#' @param choice color name
#' @param n The color numbers to use
#'
#' @return a color list
#' @export
#'
#' @examples
#' # null
getColors <- function(color.platte = NULL,
                      choice = NULL,
                      n = NULL){
  return(color.platte[[choice]][1:n])
}

globalVariables(c("num"))

#' @title 绘制细胞比例柱状图
#' @description
#' 绘制细胞比例柱状图，支持分面，源代码参考了: https://github.com/junjunlab/scRNAtoolVis/blob/master/R/cellRatioPlot.R
#' 并做了修改。
#'
#' @param object an seurat object
#' @param sample.name x axis
#' @param sample.order order for x axis
#' @param celltype.name column fill by
#' @param celltype.order order for fill by
#' @param facet.name 分面所有的列
#' @param facet.order 分面的order
#' @param col.width column width, from 0-1
#' @param flow.alpha 中间连接线的透明度
#' @param flow.curve 中间连接线的弯曲度
#' @param color.choice fill by所用的颜色系列
#' @import dplyr
#' @importFrom dplyr %>%
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' #NULL
cellRatioPlot <- function(object = NULL,
                          sample.name = NULL,
                          sample.order = NULL,
                          celltype.name = NULL,
                          celltype.order = NULL,
                          facet.name = NULL,
                          facet.order = NULL,
                          col.width = 0.7,
                          flow.alpha = 0.25,
                          flow.curve = 0,
                          color.choice = NULL) {
  requireNamespace("dplyr")
  # get meta info
  meta <- object@meta.data

  # order
  if(!is.null(sample.order)){
    meta[,sample.name] <- factor(meta[,sample.name], levels = sample.order)
  }

  if(!is.null(celltype.order)){
    meta[,celltype.name] <- factor(meta[,celltype.name], levels = celltype.order)
  }
  if(!is.null(facet.order)){
    meta[,facet.name] <- factor(meta[,facet.name], levels = facet.order)
  }
  # calculate percent ratio
  if (is.null(facet.name)) {
    ratio.info <- meta %>%
      dplyr::group_by(.data[[sample.name]], .data[[celltype.name]]) %>%
      dplyr::summarise(num = dplyr::n()) %>%
      dplyr::mutate(rel_num = num / sum(num))
  }else{
    ratio.info <- meta %>%
      dplyr::group_by(.data[[sample.name]], .data[[facet.name]],.data[[celltype.name]]) %>%
      dplyr::summarise(num = dplyr::n()) %>%
      dplyr::mutate(rel_num = num / sum(num))
  }

  # color
  fill.col <- getColors(color.platte = color_list, choice = color.choice, n = length(unique(meta[, celltype.name])))


  # plot
  p <-
    ggplot2::ggplot(
      ratio.info,
      ggplot2::aes_string(x = sample.name, y = "rel_num")
    ) +
    ggplot2::geom_col(
      ggplot2::aes_string(fill = celltype.name),
      width = col.width
    ) +
    ggalluvial::geom_flow(
      ggplot2::aes_string(
        stratum = celltype.name,
        alluvium = celltype.name,
        fill = celltype.name
      ),
      width = col.width,
      alpha = flow.alpha,
      knot.pos = flow.curve
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_fill_manual(
      values = fill.col,
      name = "Cell Type"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.2), color = "black"),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.5), color = "black"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), color = "black"),
      legend.title = ggplot2::element_text(size = ggplot2::rel(1.5), color = "black")
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("Cell percent ratio")
  if (!is.null(facet.name)) {
    p <- p +
      ggplot2::facet_grid(stats::as.formula(paste0(facet.name, "~ ."))) +
      ggplot2::theme(panel.spacing = grid::unit(0.8, "cm", data = NULL))
  }
  return(p)
}




#' @title 各个细胞内的高表达基因
#' @description
#'  对各群里的各个单细胞，分别统计top表达的基因，即基因的reads比例大于所设定的阈值expr.cut。
#'
#'
#' @param SeuratObj Seurat object
#' @param expr.cut 在单个细胞中，如果某个基因的UMI比例高于此值，则会被认为是被高表达的基因。
#' @param group.by 细胞分群参数
#'
#' @return data frame
#' @export
#'
#' @examples
#' # top_genes(cds, expr.cut = 0.01, group.by = "orig.ident")
#' # top_genes(cds, expr.cut = 0.1, group.by = "res.0.4")
#' # top_genes(cds, expr.cut = 0.1, group.by = "RandomGroup")
top_genes <- function(SeuratObj, expr.cut = 0.01, group.by) {
  requireNamespace("dplyr")
  requireNamespace("Seurat")
  # 如果是assay5类型
  if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
    SeuratObj <- JoinLayers(SeuratObj)
    counts.expr <- as.matrix(SeuratObj@assays$RNA@layers$counts)
  } else { # 如果是assay类型
    counts.expr <- as.matrix(SeuratObj@assays$RNA$counts)
  }
  colnames(counts.expr) <- colnames(SeuratObj)
  rownames(counts.expr) <- rownames(SeuratObj)
  # 分celltype计算统计值
  all.cell.types <- unique(SeuratObj@meta.data[,group.by])
  for (celltype in all.cell.types) {
    cells.sub <- colnames(SeuratObj)[as.character(SeuratObj@meta.data[,group.by]) == celltype]
    if (length(cells.sub) < 3) {
      next
    }
    counts.expr.sub <- counts.expr[,cells.sub]
    for (i in 1:ncol(counts.expr.sub)) {
      values <- sort(counts.expr.sub[, i], decreasing = TRUE)
      rates <- values/sum(values)
      top <- rates[rates > expr.cut]
      if (i == 1) {
        res <- data.frame(Gene = names(top), Expr = unname(top))
      }else{
        res <- rbind(res, data.frame(Gene = names(top), Expr = unname(top)))
      }
    }
    genes.statics <- dplyr::group_by(res, Gene) %>%
      dplyr::summarise(cut.pct.mean = round(mean(Expr),digits = 4), cut.pct.median = round(median(Expr),digits = 4), cut.Cells = length(Expr))
    genes.statics$total.pos.cells <- apply(counts.expr.sub[genes.statics$Gene,,drop = FALSE] > 0, 1, sum)
    genes.statics$total.UMI.pct <- round(apply(counts.expr.sub[genes.statics$Gene,,drop = FALSE], 1, sum)/sum(counts.expr),digits = 4)
    genes.statics$total.cells <- ncol(counts.expr.sub)
    genes.statics$celltype <- celltype
    genes.statics <- genes.statics[,c("celltype", "total.cells", "Gene", "total.pos.cells", "total.UMI.pct", "cut.Cells", "cut.pct.mean", "cut.pct.median")]
    genes.statics <- genes.statics[order(genes.statics$total.UMI.pct, decreasing = TRUE),]
    if (celltype == all.cell.types[1]) {
      results.statics <- genes.statics
    }else{
      results.statics <- rbind(results.statics, genes.statics)
    }
  }
  rownames(results.statics) <- NULL
  return(results.statics)
}

top_accumulated_genes <- function(SeuratObj, top, group.by){
  requireNamespace("dplyr")
  requireNamespace("Seurat")
  # 如果是assay5类型
  if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
    SeuratObj <- JoinLayers(SeuratObj)
    counts.expr <- as.matrix(SeuratObj@assays$RNA@layers$counts)
  } else { # 如果是assay类型
    counts.expr <- as.matrix(SeuratObj@assays$RNA$counts)
  }
  colnames(counts.expr) <- colnames(SeuratObj)
  rownames(counts.expr) <- rownames(SeuratObj)
  # 分celltype计算统计值
  all.cell.types <- unique(SeuratObj@meta.data[,group.by])
  all.cell.types <- all.cell.types[!is.na(all.cell.types)] # 特殊情况，有些celltype中有NA值
  for (celltype in all.cell.types) {
    cells.sub <- Cells(SeuratObj)[(as.character(SeuratObj@meta.data[, group.by]) == as.character(celltype)) & !is.na(SeuratObj@meta.data[, group.by])] # 特殊情况，有些celltype中有NA值
    if (length(cells.sub) < 3) {
      next
    }
    counts.expr.sub <- counts.expr[,cells.sub,drop = FALSE]
    sss <- sort(apply(counts.expr.sub, 1, sum),decreasing = TRUE)
    if (length(sss) > top) {
      sss <- sss[1:top]
    }
    res <- data.frame(Gene = names(sss), NormalizedAccumulatedUMICounts = unname(sss)/ncol(counts.expr.sub), PCT = round(unname(sss)/sum(counts.expr.sub),digits = 4))
    res$total.pos.cells <- unname(apply(counts.expr.sub[res$Gene,,drop = FALSE] > 0, 1, sum))
    res$total.cells <- ncol(counts.expr.sub)
    res$celltype <- celltype
    res <- res[,c("celltype", "total.cells", "Gene", "total.pos.cells", "NormalizedAccumulatedUMICounts", "PCT")]
    res <- res[order(res$AccumulatedUMICounts, decreasing = TRUE),]
    if (celltype == all.cell.types[1]) {
      results.statics <- res
    }else{
      results.statics <- rbind(results.statics, res)
    }
  }
  rownames(results.statics) <- NULL
  return(results.statics)
}


summary_features <- function(SeuratObj, features, group.by){
  requireNamespace("dplyr")
  requireNamespace("Seurat")
  # 如果是assay5类型
  if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
    SeuratObj <- JoinLayers(SeuratObj)
    normalized.expr <- as.matrix(SeuratObj@assays$RNA@layers$data)
  } else { # 如果是assay类型
    normalized.expr <- as.matrix(SeuratObj@assays$RNA$data)
  }
  all.cell.types <- unique(SeuratObj@meta.data[,group.by])
  for (celltype in all.cell.types) {
    cells.sub <- colnames(SeuratObj)[as.character(SeuratObj@meta.data[,group.by]) == celltype]
    if (length(cells.sub) < 3) {
      next
    }
    normalized.expr.sub <- normalized.expr[,cells.sub]
    mean.expr <- apply(normalized.expr.sub[features,,drop = FALSE], 1, mean)
    median.expr <- apply(normalized.expr.sub[features,,drop = FALSE], 1, median)
    pct <- apply(normalized.expr.sub[features,,drop = FALSE] > 0, 1, mean)
    single.res <- data.frame(Gene = features, Expr.mean = round(mean.expr,digits = 4), Expr.median = round(median.expr, digits = 4), PCT = round(pct,digits = 4))
    single.res$CellType <- celltype
    single.res$TotalCells <- ncol(normalized.expr.sub)
    single.res <- single.res[,c("CellType", "TotalCells","Gene", "PCT", "Expr.mean", "Expr.median")]
    rownames(single.res) <- NULL
    if (celltype == all.cell.types[1]) {
      res <- single.res
    }else{
      res <- rbind(res, single.res)
    }
  }
  rownames(res) <- NULL
  return(res)
}

calculate_top_correlations <- function(SeuratObj, method, top = 1000){
  if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
    SeuratObj <- JoinLayers(SeuratObj)
    normalized.expr <- as.matrix(SeuratObj@assays$RNA@layers$data)
  } else { # 如果是assay类型
    normalized.expr <- as.matrix(SeuratObj@assays$RNA$data)
  }
  # 过滤细胞，表达值的和要大于细胞数目的1/10，设置比较随意！
  normalized.expr <- normalized.expr[apply(normalized.expr, 1, sum) > 0.1 * ncol(SeuratObj),,drop = FALSE]
  cor.res = cor(t(as.matrix(normalized.expr)), method = method)
  cor.res[lower.tri(cor.res, diag = TRUE)] <- 0
  cor.res <- reshape2::melt(cor.res)
  colnames(cor.res) <- c("GeneA","GeneB","correlation")
  cor.res <- cor.res[order(abs(cor.res$correlation),decreasing = TRUE),]
  cor.res <- cor.res[cor.res$correlation != 0,,drop = FALSE]
  cor.res$correlation <- round(cor.res$correlation, digits = 4)
  if (nrow(cor.res) > top) {
    cor.res <- cor.res[1:top, ]
  }
  rownames(cor.res) <- NULL
  return(cor.res)
}

calculate_most_correlated <- function(SeuratObj, feature, method){
  if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
    SeuratObj <- JoinLayers(SeuratObj)
    normalized.expr <- as.matrix(SeuratObj@assays$RNA@layers$data)
  } else { # 如果是assay类型
    normalized.expr <- as.matrix(SeuratObj@assays$RNA$data)
  }
  x <- normalized.expr[feature, ,drop = FALSE]
  y <- normalized.expr[rownames(normalized.expr)[rownames(normalized.expr) != feature],,drop = FALSE]
  # 过滤细胞，表达值的和要大于细胞数目的1/10，设置比较随意！
  y <- y[apply(y, 1, sum) > 0.1 * ncol(SeuratObj),,drop = FALSE]
  cor.res = cor(x = t(as.matrix(x)), y =  t(as.matrix(y)), method = method)
  cor.res <- reshape2::melt(cor.res)
  colnames(cor.res) <- c("GeneA","GeneB","correlation")
  cor.res <- cor.res[order(abs(cor.res$correlation),decreasing = TRUE),]
  cor.res$correlation <- round(cor.res$correlation, digits = 4)
  rownames(cor.res) <- NULL
  return(cor.res)
}

calculate_correlation <- function(SeuratObj, features, method){
  if (class(SeuratObj@assays$RNA)[1] == "Assay5") {
    SeuratObj <- JoinLayers(SeuratObj)
    normalized.expr <- as.matrix(SeuratObj@assays$RNA@layers$data)
  } else { # 如果是assay类型
    normalized.expr <- as.matrix(SeuratObj@assays$RNA$data)
  }
  normalized.expr <- normalized.expr[rownames(normalized.expr) %in% features,,drop = FALSE]
  # 过滤细胞，表达值的和要大于细胞数目的1/10，设置比较随意！
  normalized.expr <- normalized.expr[apply(normalized.expr, 1, sum) > 0.1 * ncol(SeuratObj),,drop = FALSE]
  cor.res = cor(t(as.matrix(normalized.expr)), method = method)
  cor.res[lower.tri(cor.res, diag = TRUE)] <- 0
  cor.res <- reshape2::melt(cor.res)
  colnames(cor.res) <- c("GeneA","GeneB","correlation")
  cor.res <- cor.res[order(abs(cor.res$correlation),decreasing = TRUE),]
  cor.res <- cor.res[cor.res$correlation != 0, ,drop = FALSE]
  cor.res$correlation <- round(cor.res$correlation, digits = 4)
  rownames(cor.res) <- NULL
  return(cor.res)
}


# refer to: https://github.com/junjunlab/scRNAtoolVis/blob/master/R/averageHeatmap.R
AverageHeatmap <- function(
    object = NULL,
    markerGene = NULL,
    group.by = "ident",
    assays = "RNA",
    slot = "data",
    htCol = c("#0099CC", "white", "#CC0033"),
    colseed = 666,
    htRange = c(-2, 0, 2),
    annoCol = FALSE,
    myanCol = NULL,
    annoColType = "light",
    annoColTypeAlpha = 0,
    clusterAnnoName = TRUE,
    showRowNames = TRUE,
    row_names_side = "left",
    markGenes = NULL,
    border = FALSE,
    cluster.fontsize = 12,
    feature.fontsize = 10,
    column_names_rot = 45,
    width = NULL,
    height = NULL,
    cluster.order = NULL,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    gene.order = NULL,
    ...) {
  # BiocManager::install("ComplexHeatmap")
  requireNamespace("ComplexHeatmap")

  # get cells mean gene expression
  # check Seurat version first
  vr <- utils::compareVersion(as.character(utils::packageVersion("Seurat")),"5")
  if(vr == 1){
    mean_gene_exp <- as.matrix(
      data.frame(
        Seurat::AverageExpression(object,
                                  features = markerGene,
                                  group.by = group.by,
                                  assays = assays,
                                  layer = slot
        )
      )
    )
  }else{
    mean_gene_exp <- as.matrix(
      data.frame(
        Seurat::AverageExpression(object,
                                  features = markerGene,
                                  group.by = group.by,
                                  assays = assays,
                                  slot = slot
        )
      )
    )
  }

  # add colnames
  # name1 <- gsub(
  #   pattern = paste0(assays, ".", sep = ""),
  #   replacement = "",
  #   colnames(mean_gene_exp)
  # )
  #
  # colnames(mean_gene_exp) <- gsub(
  #   pattern = "\\.",
  #   replacement = " ", name1
  # )

  colnames(mean_gene_exp) <- levels(Seurat::Idents(object))

  # Z-score
  htdf <- t(scale(t(mean_gene_exp), scale = TRUE, center = TRUE))

  # cluster order
  if (!is.null(cluster.order)) {
    htdf <- htdf[, cluster.order]
  }

  # gene order
  if (!is.null(gene.order)) {
    htdf <- htdf[gene.order, ]
  }

  # color
  col_fun <- circlize::colorRamp2(htRange, htCol)

  # anno color
  if (annoCol == FALSE) {
    set.seed(colseed)
    anno_col <- circlize::rand_color(
      ncol(htdf),
      luminosity = annoColType,
      transparency = annoColTypeAlpha
    )
    print(c("Your cluster annotation color is:", anno_col))
  } else if (annoCol == TRUE) {
    # give your own color vectors
    anno_col <- myanCol
  } else {
    print("Give TRUE or FALSE paramters!")
  }
  names(anno_col) <- colnames(htdf)

  # top annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    cluster = colnames(htdf),
    show_legend = FALSE,
    show_annotation_name = clusterAnnoName,
    col = list(cluster = anno_col)
  )

  # whether mark your genes on plot
  if (!is.null(markGenes)) {
    # all genes
    rowGene <- rownames(htdf)

    # tartget gene
    annoGene <- markGenes

    # get target gene index
    index <- match(annoGene, rowGene)

    # some genes annotation
    geneMark <- ComplexHeatmap::rowAnnotation(
      gene = ComplexHeatmap::anno_mark(
        at = index,
        labels = annoGene,
        labels_gp = grid::gpar(
          fontface = "italic",
          fontsize = feature.fontsize
        )
      ),
      ...
    )

    right_annotation <- geneMark
  } else {
    right_annotation <- NULL
  }

  # control heatmap width and height
  if (is.null(width) || is.null(height)) {
    # plot
    ComplexHeatmap::Heatmap(
      htdf,
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      name = "Z-score",
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows,
      # column_title = "Clusters",
      right_annotation = right_annotation,
      show_row_names = showRowNames,
      column_names_gp = grid::gpar(fontsize = cluster.fontsize),
      row_names_gp = grid::gpar(
        fontface = "italic",
        fontsize = feature.fontsize
      ),
      row_names_side = row_names_side,
      border = border,
      column_names_side = "top",
      column_names_rot = column_names_rot,
      top_annotation = column_ha,
      col = col_fun,
      ...
    )
  } else {
    # plot
    ComplexHeatmap::Heatmap(
      htdf,
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      name = "Z-score",
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      # column_title = "Clusters",
      right_annotation = right_annotation,
      show_row_names = showRowNames,
      column_names_gp = grid::gpar(fontsize = cluster.fontsize),
      row_names_gp = grid::gpar(
        fontface = "italic",
        fontsize = feature.fontsize
      ),
      row_names_side = row_names_side,
      border = border,
      column_names_side = "top",
      column_names_rot = column_names_rot,
      top_annotation = column_ha,
      col = col_fun,
      width = ggplot2::unit(width, "cm"),
      height = ggplot2::unit(height, "cm"),
      ...
    )
  }
}

readSeurat <- function(path){
  message("Reading in data...")
  # read data
  if (tools::file_ext(path) == 'qs2') {
    seu_obj <- qs2::qs_read(path)
  }else(
    seu_obj <- readRDS(path)
  )
  # check version
  if (SeuratObject::Version(seu_obj) < packageVersion('SeuratObject')) {
    seu_obj <- SeuratObject::UpdateSeuratObject(seu_obj)
  }else{
    message('Update Seurat Object escaped for it has been the latest version!')
  }
  return(seu_obj)
}

