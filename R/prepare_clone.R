#' Prepare Clone
#'
#' @param db db
#' @param junction junciotn seq
#' @param v_call v_call
#' @param j_call j_call
#' @param first first allele or union
#' @param cdr3 cdr3 only
#' @param fields additional fields eg subject
#' @param cell_id cell id
#' @param locus locus column
#' @param only_heavy default true
#' @param mod3 mod3 seq only
#' @param max_n max
#' @returns list
#' @export
prepare_clone <- function(db,
                          junction = "junction", v_call = "v_call", j_call = "j_call",
                          first = FALSE, cdr3 = FALSE, fields = NULL,
                          cell_id = NULL, locus = "locus", only_heavy = TRUE,
                          mod3 = FALSE, max_n = 0) {
  # add junction length column
  db$junction_l <- stringi::stri_length(db[[junction]])
  junction_l <- "junction_l"

  ### check for mod3
  # filter mod 3 junction lengths
  if (mod3) {
    n_rmv_mod3 <- sum(db[[junction_l]]%%3 != 0)
    db <- db %>%
      dplyr::filter(!!rlang::sym(junction_l)%%3 == 0)
  } else {
    n_rmv_mod3 <- 0
  }

  ### check for cdr3
  # filter junctions with length > 6
  if (cdr3) {
    n_rmv_cdr3 <- sum(db[[junction_l]] <= 6)
    db <- db %>%
      dplyr::filter(!!rlang::sym(junction_l) > 6)
    # add cdr3 column
    db$cdr3_col <- substr(db[[junction]], 4, db[[junction_l]]-3)
    cdr3_col <- "cdr3_col"
    #cdr3_col = db$cdr3_col
  } else {
    n_rmv_cdr3 <- 0
    cdr3_col <- NA
  }

  ### check for degenerate characters (non-ATCG's)
  # Count the number of non-ATCG's in junction
  if (!is.null(max_n)) {
    n_rmv_N <- sum(stringi::stri_count(db[[junction]], regex = "[^ATCG]") > max_n)
    db <- db %>%
      dplyr::filter(stringi::stri_count(!!rlang::sym(junction), regex = "[^ATCG]") <= max_n)
  } else {
    n_rmv_N <- 0
  }

  ### Parse V and J columns to get gene
  if (!is.null(fields)) {
    . <- NULL
    db <- db %>%
      dplyr::group_by(!!!rlang::syms(fields)) %>%
      do(alakazam::groupGenes(.,
                              v_call = v_call,
                              j_call = j_call,
                              junc_len = NULL,
                              cell_id = cell_id,
                              locus = locus,
                              only_heavy = only_heavy,
                              first = first))
  } else {
    db <- alakazam::groupGenes(db,
                               v_call = v_call,
                               j_call = j_call,
                               junc_len = NULL,
                               cell_id = cell_id,
                               locus = locus,
                               only_heavy = only_heavy,
                               first = first)
  }

  ### groups to use
  groupBy <- c("vj_group", junction_l, fields)

  ### assign group ids to db
  db$vjl_group <- db %>%
    dplyr::group_by(!!!rlang::syms(groupBy)) %>%
    dplyr::group_indices()

  ### return results
  return_list <- list("db" = db,
                      "n_rmv_mod3" = n_rmv_mod3,
                      "n_rmv_cdr3" = n_rmv_cdr3,
                      "n_rmv_N" = n_rmv_N,
                      "junction_l" = junction_l,
                      "cdr3_col" =  cdr3_col)
  return(return_list)
}





bcrCounts <- function(data, inds) {
  dbj = data %>% dplyr::group_by(.data[[inds]]) %>%
    dplyr::summarise(n = n())
  dbj$p = dbj$n / sum(dbj$n)
  return(dbj)
}


bcrCounts_pheno <- function(data, inds, pheno) {
  dbj = data %>% dplyr::group_by(.data[[inds]], .data[[pheno]]) %>%
    dplyr::summarise(n = n())
  dbj$p = dbj$n / sum(dbj$n)
  return(dbj)
}


#
# phenocounts <- function(data, pheno, inds) {
#   dbp = data %>% dplyr::group_by(.data[[inds]], .data[[pheno]]) %>%
#     dplyr::summarise(n = n()) %>%
#     dplyr::mutate(p = n / sum(n)) %>%
#     tidyr::spread(key = {{pheno}}, value = p, fill=0)
#   return(dbp)
# }






phenocounts <- function(data, pheno, inds) {
  dbp = data %>% dplyr::group_by(.data[[inds]], .data[[pheno]]) %>%
    dplyr::summarise(n = n()) %>%
    #dplyr::mutate(p = n / sum(n)) %>%
    tidyr::pivot_wider(names_from = {{pheno}}, values_from = n, fill=0)
  return(dbp)
}



get_jd <- function(db, pheno, indVar="ind", w="w") {
  dbp = db %>% dplyr::group_by(.data[[indVar]], .data[[w]],  .data[[pheno]]) %>%
    dplyr::summarise(n = n()) %>%
    tidyr::spread(key = {{pheno}}, value = n, fill=0)
  return(dbp)
}




get_jd_p <- function(db, pheno, indVar="ind", w="w") {
  dbp = db %>% dplyr::group_by(.data[[indVar]], .data[[w]], .data[[pheno]]) %>%
    dplyr::summarise(n = n() )
  dbp$p = dbp$n / sum(dbp$n)
  dbp = dbp %>%  tidyr::spread(key = {{pheno}}, value = p, fill=0)
  return(dbp)
}



