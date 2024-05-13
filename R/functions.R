prepare_seurat_object <- function(obj){
  requireNamespace("Seurat")
  # 1. meta.data中，分析非factor类型的列，如果unique数目少于细胞总数的1/50，并且不大于50种，也强制转为factor类型。
  # change some columns in meta.data to factors
  columns.to.factor <- colnames(obj@meta.data)[check_df_level(obj@meta.data) & !check_df_factor(obj@meta.data)]
  obj@meta.data[columns.to.factor] <- lapply(obj@meta.data[columns.to.factor],as.factor)
  return(obj)
}

# check data.frame all columns levels, return a logic vector
check_df_level <- function(df, max.counts = 50, max.pct = 1/50){
  cutoff <- min(max.counts, round(nrow(df) * max.pct))
  unname(sapply(df,function(x){length(unique(as.character(x)))})) <= cutoff
}
# check if data.frame all columns are factor, return a logic vector
check_df_factor <- function(df){
  unname(sapply(df,is.factor))
}

# data.frame中，类型为factor，且level数目不超过指定值的列有哪些，返回列名
df_factor_columns <- function(df, max.level = 4){
  df <- df[,check_df_factor(df)]
  colnames(df)[unname(sapply(df,FUN = function(x)length(levels(x)))) <= max.level]
}

# Chcek the input gene, return the revised gene
CheckGene <- function(InputGene, GeneLibrary){
  InputGene <- trimws(InputGene)
  if (InputGene %in% GeneLibrary) { # when input gene is absolutely right
    return(InputGene)
  }else if(tolower(InputGene) %in% tolower(GeneLibrary)){ # when not case sensitive
    return(GeneLibrary[tolower(GeneLibrary) %in% tolower(InputGene)]) # gene list length can > 1
  }else{ # when not match
    return(NA)
  }
}
