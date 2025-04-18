#' Genes information extracted from genome gtf file
#'
#' A data set containing the gene information from species of human, mouse and fly.
#'
#' @format A list:
#' \describe{
#'   \item{fly}{a data frame, with columns: "Ensembl","Symbol", "EntrezID", "UniProt", "Name", "Alias"}
#'   \item{mouse}{a data frame, with columns: "Ensembl","Symbol", "MGI", "EntrezID", "UniProt", "Name", "Alias"}
#'   \item{human}{a data frame, with columns: "Ensembl","Symbol", "HGNC", "EntrezID", "UniProt", "Name", "Alias"}
#' }
"GenesDB"
