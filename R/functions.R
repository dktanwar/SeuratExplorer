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
  InputGenes <- unlist(strsplit(InputGene,split = " "))
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


#' @title useMyCol
#' @name useMyCol
#' @author Junjun Lao
#' @description This function is used to produce avaliable colors for plot.
#' @param platte The platte name. Default("stallion").
#' @param n The color numbers to use. Default(NULL).
#' @param showAll Whether to show all plattes. Default(FALSE).
#'
#' @return Return the color names you have choosed.
#' @export
#'
#' @examples
#' useMyCol(platte = 'stallion2',n = 5)
#' useMyCol(showAll = TRUE)

# define functions
useMyCol <- function(platte = NULL,
                     n = NULL,
                     showAll = FALSE){
  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------
  # ============================================================================
  #20-colors
  stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
               "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
               "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
               "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")

  stallion2 = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                "#6E4B9E","#0C727C", "#7E1416","#D8A767")

  calm = c("#7DD06F", "#844081","#688EC1", "#C17E73", "#484125",
           "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
           "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D",
           "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736")

  kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020",
            "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
            "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800",
            "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16")

  #16-colors
  bear = c("#faa818", "#41a30d","#fbdf72", "#367d7d", "#d33502",
           "#6ebcbc", "#37526d","#916848", "#f5b390", "#342739",
           "#bed678","#a6d9ee", "#0d74b6",
           "#60824f","#725ca5", "#e0598b")

  #15-colors
  ironMan = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',
              "#F88D51","#FAD510","#FFFF5F",'#88CFA4','#238B45',
              "#02401B", "#0AD7D3","#046C9A", "#A2A475", 'grey35')

  circus = c("#D52126","#88CCEE", "#FEE52C", "#117733", "#CC61B0",
             "#99C945", "#2F8AC4", "#332288","#E68316", "#661101",
             "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74")

  #12-colors
  paired = c("#A6CDE2","#1E78B4","#74C476","#34A047","#F59899","#E11E26",
             "#FCBF6E","#F47E1F","#CAB2D6","#6A3E98","#FAF39B","#B15928")

  #11-colors
  grove = c("#1a1334","#01545a","#017351","#03c383","#aad962",
            "#fbbf45","#ef6a32","#ed0345","#a12a5e","#710162","#3B9AB2")

  #7-colors
  summerNight = c("#2a7185","#a64027","#fbdf72","#60824f","#9cdff0",
                  "#022336","#725ca5")

  #5-colors
  zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
  darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
  rushmore = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" , "#F2300F")
  captain = c("grey","#A1CDE1","#12477C","#EC9274","#67001E")

  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------
  #10-colors
  horizon = c('#000075','#2E00FF', '#9408F7', '#C729D6', '#FA4AB5',
              '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60')

  #9-colors
  horizonExtra =c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5",
                  "#FD619D","#FF9965","#FFD32B","#FFFC5A")

  blueYellow = c("#352A86","#343DAE","#0262E0","#1389D2","#2DB7A3",
                 "#A5BE6A","#F8BA43","#F6DA23","#F8FA0D")

  sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                 '#A2E700','#FFFF00','#FFD200','#FFA500')

  solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC',
                 '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D')

  whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b')

  whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                '#3690c0','#0570b0','#045a8d','#023858')


  comet = c("#E6E7E8","#3A97FF","#8816A7","black")

  #7-colors
  greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                '#0868ac','#084081')

  #6-colors
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30",
            "#F7962E","#FCEE2B")

  #5-colors
  coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29")
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A")
  greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF")
  fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A")
  purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F")

  # ============================================================================
  if(showAll == FALSE){
    if(platte == "stallion"){
      col <- stallion[1:n]
      message(paste0('This palatte have ',length(stallion),' colors!'))
    }else if(platte == "stallion2"){
      col <- stallion2[1:n]
      message(paste0('This palatte have ',length(stallion2),' colors!'))
    }else if(platte == "calm"){
      col <- calm[1:n]
      message(paste0('This palatte have ',length(calm),' colors!'))
    }else if(platte == "kelly"){
      col <- kelly[1:n]
      message(paste0('This palatte have ',length(kelly),' colors!'))
    }else if(platte == "bear"){
      col <- bear[1:n]
      message(paste0('This palatte have ',length(bear),' colors!'))
    }else if(platte == "ironMan"){
      col <- ironMan[1:n]
      message(paste0('This palatte have ',length(ironMan),' colors!'))
    }else if(platte == "circus"){
      col <- circus[1:n]
      message(paste0('This palatte have ',length(circus),' colors!'))
    }else if(platte == "paired"){
      col <- paired[1:n]
      message(paste0('This palatte have ',length(paired),' colors!'))
    }else if(platte == "grove"){
      col <- grove[1:n]
      message(paste0('This palatte have ',length(grove),' colors!'))
    }else if(platte == "summerNight"){
      col <- summerNight[1:n]
      message(paste0('This palatte have ',length(summerNight),' colors!'))
    }else if(platte == "zissou"){
      col <- zissou[1:n]
      message(paste0('This palatte have ',length(zissou),' colors!'))
    }else if(platte == "darjeeling"){
      col <- darjeeling[1:n]
      message(paste0('This palatte have ',length(darjeeling),' colors!'))
    }else if(platte == "rushmore"){
      col <- rushmore[1:n]
      message(paste0('This palatte have ',length(rushmore),' colors!'))
    }else if(platte == "captain"){
      col <- captain[1:n]
      message(paste0('This palatte have ',length(captain),' colors!'))
    }else if(platte == "horizon"){
      col <- horizon[1:n] # continues colors
      message(paste0('This palatte have ',length(horizon),' colors!'))
    }else if(platte == "horizonExtra"){
      col <- horizonExtra[1:n]
      message(paste0('This palatte have ',length(horizonExtra),' colors!'))
    }else if(platte == "blueYellow"){
      col <- blueYellow[1:n]
      message(paste0('This palatte have ',length(blueYellow),' colors!'))
    }else if(platte == "sambaNight"){
      col <- sambaNight[1:n]
      message(paste0('This palatte have ',length(sambaNight),' colors!'))
    }else if(platte == "solarExtra"){
      col <- solarExtra[1:n]
      message(paste0('This palatte have ',length(solarExtra),' colors!'))
    }else if(platte == "whitePurple"){
      col <- whitePurple[1:n]
      message(paste0('This palatte have ',length(whitePurple),' colors!'))
    }else if(platte == "whiteBlue"){
      col <- whiteBlue[1:n]
      message(paste0('This palatte have ',length(whiteBlue),' colors!'))
    }else if(platte == "comet"){
      col <- comet[1:n]
      message(paste0('This palatte have ',length(comet),' colors!'))
    }else if(platte == "greenBlue"){
      col <- greenBlue[1:n]
      message(paste0('This palatte have ',length(greenBlue),' colors!'))
    }else if(platte == "beach"){
      col <- beach[1:n]
      message(paste0('This palatte have ',length(beach),' colors!'))
    }else if(platte == "coolwarm"){
      col <- coolwarm[1:n]
      message(paste0('This palatte have ',length(coolwarm),' colors!'))
    }else if(platte == "fireworks"){
      col <- fireworks[1:n]
      message(paste0('This palatte have ',length(fireworks),' colors!'))
    }else if(platte == "greyMagma"){
      col <- greyMagma[1:n]
      message(paste0('This palatte have ',length(greyMagma),' colors!'))
    }else if(platte == "fireworks2"){
      col <- fireworks2[1:n]
      message(paste0('This palatte have ',length(fireworks2),' colors!'))
    }else if(platte == "purpleOrange"){
      col <- purpleOrange[1:n]
      message(paste0('This palatte have ',length(purpleOrange),' colors!'))
    }else{
      message('Please give the correct name!')
    }
    return(col)
  }else{
    dicrete_col <- c('stallion', 'stallion2', 'calm', 'kelly', 'bear' ,
                     'ironMan', 'circus', 'paired', 'grove', 'summerNight' ,
                     'zissou', 'darjeeling', 'rushmore', 'captain')

    continues_col <- c('horizon', 'horizonExtra', 'blueYellow', 'sambaNight', 'solarExtra',
                       'whitePurple', 'whiteBlue', 'comet', 'greenBlue', 'beach' ,
                       'coolwarm', 'fireworks', 'greyMagma', 'fireworks2', 'purpleOrange')

    message(cat("discrete colors: ",dicrete_col))
    message(cat("continues colors: ",continues_col))
  }
}

#' @title 绘制细胞比例柱状图
#' @description
#' 绘制细胞比例柱状图，支持分面，源代码参考了: https://github.com/junjunlab/scRNAtoolVis/blob/master/R/cellRatioPlot.R，并做了修改。
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

#' @param fill.col fill by所用的颜色
#'
#' @importFrom dplyr %>%
#' @return
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
                          fill.col = NULL) {
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
      dplyr::summarise(num = n()) %>%
      dplyr::mutate(rel_num = num / sum(num))
  }else{
    ratio.info <- meta %>%
      dplyr::group_by(.data[[sample.name]], .data[[facet.name]],.data[[celltype.name]]) %>%
      dplyr::summarise(num = n()) %>%
      dplyr::mutate(rel_num = num / sum(num))
  }

  # color
  if (is.null(fill.col)) {
    fill.col <- useMyCol("paired", n = length(unique(meta[, celltype.name])))
  } else {
    fill.col <- fill.col
  }

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
    p <- p + ggplot2::facet_grid(stats::as.formula(paste0(facet.name, "~ .")))
  }
  return(p)
}

