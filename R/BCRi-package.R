#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @import scoper
#' @import igraph
#' @import      ggplot2
#' @import      graphics
#' @import      methods
#' @import      utils
#' @importFrom  alakazam    pairwiseDist checkColumns getDNAMatrix getAAMatrix
#'                          padSeqEnds progressBar groupGenes baseTheme translateDNA
#'                          getLocus
#' @importFrom  data.table  as.data.table .I rbindlist
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  dplyr       n %>% do
#'                          filter select arrange bind_rows
#'                          group_by ungroup group_indices
#'                          mutate summarize slice  distinct left_join
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  rlang       sym syms
#' @importFrom  Rcpp        evalCpp
#' @importFrom  scales      pretty_breaks
#' @importFrom  shazam      consensusSequence
#' @importFrom  stats       density kmeans sd cor
#'                          as.dist hclust cutree setNames
#' @importFrom  stringi     stri_split_fixed stri_length stri_count stri_join
#' @importFrom  tidyr       separate
#' @useDynLib BCRi, .registration = TRUE
## usethis namespace: end
NULL
