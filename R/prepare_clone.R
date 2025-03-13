#' Prepare Clone
#'
#' @param db
#' @param junction
#' @param v_call
#' @param j_call
#' @param first
#' @param cdr3
#' @param fields
#' @param cell_id
#' @param locus
#' @param only_heavy
#' @param mod3
#' @param max_n
#'
#' @returns
#' @export
#'
#' @examples
prepare_clone <- function(db,
                          junction = "junction", v_call = "v_call_genotyped", j_call = "j_call",
                          first = FALSE, cdr3 = FALSE, fields = NULL,
                          cell_id = NULL, locus = NULL, only_heavy = TRUE,
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

