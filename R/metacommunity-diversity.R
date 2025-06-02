#' Raw metacommunity alpha diversity
#'
#' Calculates similarity-sensitive raw metacommunity alpha diversity (the
#' naive-community metacommunity diversity). This measure may be calculated
#' for a series of orders, represented as a vector of \code{qs}.
#'
#' @param meta object of class \code{metacommunity}
#' @param qs \code{vector} of mode \code{numeric} containing \emph{q} values
#'
#' @return \code{raw_meta_alpha} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate raw metacommunity alpha diversity
#' raw_meta_alpha(meta, 0:2)
#'
raw_meta_alpha <- function(meta, qs)
  metadiv(raw_alpha(meta), qs)

#' Normalised metacommunity alpha diversity
#'
#' Calculates similarity-sensitive normalised metacommunity alpha diversity
#' (the average similarity-sensitive diversity of subcommunities). This
#' measure may be calculated for a series of orders, represented as a vector
#' of \code{qs}.
#'
#' @inheritParams raw_meta_alpha
#'
#' @return \code{norm_meta_alpha} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate normalised metacommunity alpha diversity
#' norm_meta_alpha(meta, 0:2)
#'
norm_meta_alpha <- function(meta, qs)
  metadiv(norm_alpha(meta), qs)

#' Raw metacommunity beta diversity
#'
#' Calculates similarity-sensitive raw metacommunity beta diversity (the
#' average distinctiveness of subcommunities). This  measure may be
#' calculated for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_meta_alpha
#'
#' @return \code{raw_meta_beta} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate raw metacommunity beta diversity
#' raw_meta_beta(meta, 0:2)
#'
raw_meta_beta <- function(meta, qs)
  metadiv(raw_beta(meta), qs)

#' Normalised metacommunity beta diversity
#'
#' Calculates similarity-sensitive normalised metacommunity beta diversity
#' (the effective number of distinct subcommunities. This measure may be
#' calculated for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_meta_alpha
#'
#' @return \code{norm_meta_beta} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate normalised metacommunity beta diversity
#' norm_meta_beta(meta, 0:2)
#'
norm_meta_beta <- function(meta, qs)
  metadiv(norm_beta(meta), qs)

#' Raw metacommunity rho diversity
#'
#' Calculates similarity-sensitive raw metacommunity rho diversity (the
#' average redundancy of subcommunities. This measure may be calculated
#' for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_meta_alpha
#'
#' @return \code{raw_meta_rho} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate metacommunity rho diversity
#' raw_meta_rho(meta, 0:2)
#'
raw_meta_rho <- function(meta, qs)
  metadiv(raw_rho(meta), qs)

#' Normalised metacommunity rho diversity
#'
#' Calculates similarity-sensitive normalised metacommunity rho diversity (the
#' average representativeness of subcommunities. This measure may be
#' calculated for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_meta_alpha
#'
#' @return \code{norm_meta_rho} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate normalised metacommunity rho diversity
#' norm_meta_rho(meta, 0:2)
#'
norm_meta_rho <- function(meta, qs)
  metadiv(norm_rho(meta), qs)

#' Metacommunity gamma diversity
#'
#' Calculates similarity-sensitive metacommunity gamma diversity (the
#' metacommunity similarity-sensitive diversity). This measure may be
#' calculated for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_meta_alpha
#'
#' @return \code{meta_gamma} returns a standard output of class \code{rdiv}
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate metacommunity gamma diversity
#' meta_gamma(meta, 0:2)
#'
meta_gamma <- function(meta, qs)
  metadiv(raw_gamma(meta), qs)