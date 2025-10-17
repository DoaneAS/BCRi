#' Check partition matrix
#'
#' \code{check_partition()} is used to validate partition matrices.
#'
#' @param partition two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#'
#' @return Returns a two-dimensions \code{matrix} of mode \code{numeric}. If
#' the partition matrix was valid, this should be identical to that which was
#' input as an argument.
#'
#' @noRd
#'
check_partition <- function(partition) {
  if (is.vector(partition)) partition <- as.matrix(partition)
  if (is.data.frame(partition)) partition <- as.matrix(partition)

  # normalise partition if it does not sum to 1
  if (!isTRUE(all.equal(sum(partition), 1))) {
    partition <- partition / sum(partition)
    ##message("Metacommunity matrix was normalised to sum to 1.")
  }

  partition
}

#' Check similarity matrix
#'
#' \code{check_similarity()} is used to validate similarity matrices.
#'
#' @param similarity two-dimensional \code{matrix} of mode \code{numeric};
#' contains pair-wise similarity between types.
#' @param partition two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#'
#' @return Returns a two-dimensions \code{matrix} of mode \code{numeric}. If
#' the similarity matrix was valid, this should be identical to that which was
#' input as an argument.
#'
#' @noRd
#'
check_similarity <- function(similarity, partition) {
  if (is.data.frame(similarity)) similarity <- as.matrix(similarity)

  if (any(similarity[!is.na(similarity)] < 0))
    stop("similarity matrix elements must take positive values.")

  if (ncol(similarity) != nrow(similarity))
    stop("similarity matrix must be square.")

  if (!missing(partition)) {
    if (nrow(similarity) != nrow(partition))
      stop("similarity and partition matrices must have equal types.")

    if (is.null(row.names(similarity))){
      row.names(similarity) <- row.names(partition)
      colnames(similarity) <- row.names(partition)
    }
  }

  similarity
}

#' Power mean of vector elements
#'
#' \code{power_mean()} calculates the power mean of a set of values.
#'
#' Calculates the order-th power mean of a single
#' set of non-negative values, weighted by weights; by default, weights are
#' equal and order is 1, so this is just the arithmetic mean. Equal weights
#' and a order of 0 gives the geometric mean, and an order of -1 gives the
#' harmonic mean.
#'
#' @param values Values for which to calculate mean.
#' @param order Order of power mean.
#' @param weights Weights of elements, normalised to 1 inside function.
#'
#' @return Weighted power mean
#' @export
#'
#' @examples
#' values <- sample(1:50, 5)
#' power_mean(values)
#'
power_mean <-
  function(values, order = 1, weights = rep(1, length(values))) {
    # Number of values must equal the number of weights
    stopifnot(length(values) == length(weights))

    # Values and weights must be greater than zero
    if (any(values[!is.nan(values)] < 0)) stop("values must be greater than 0")
    if (any(weights[!is.nan(weights)] < 0)) stop("weights must be greater than 0")

    # Normalise weights to sum to 1 (as per Rényi)
    proportions <- weights / sum(weights)

    # Check whether all proportions are NaN - happens in normalisation when all
    # weights are zero in group. In that case we want to propagate the NaN
    if (all(is.nan(proportions)))
      return(NaN)

    # Otherwise NaNs should only occur (in values) when weight is 0 and so will be
    # stripped out here as we have to eliminate non-zero weights
    non.zero <- weights > 0
    values <- values[non.zero]
    proportions <- proportions[non.zero]

    # Avoid rounding errors for order 0
    if (abs(order) < .Machine$double.eps ^ 0.5) {
      prod(values ^ proportions)
    } else if (is.infinite(order)) {
      if (order > 0)
        max(values)
      else
        min(values)
    } else {
      sum(proportions * values ^ order) ^ (1 / order)
    }
  }

#' Summary function
#'
#' This function converts columns of an array (each representing community
#' counts) into proportions, so that each column sums to 1.
#'
#' @param populations An S x N array whose columns are counts of individuals.
#' @param normalise Normalise probability distribution to sum to 1 for each
#' column rather than just along each set.
#'
#' @return Returns an array whose columns are proportions.
#'
#' @noRd
#'
summarise <-
  function(populations, normalise = TRUE) {
    totals <- array(rowSums(populations), dim = c(dim(populations)[1], 1))

    if (normalise) {
      total <- sum(totals)
      totals <- totals / total
      proportions <- populations / total
      weights <- colSums(proportions)
      proportions <- proportions %*% diag(1 / (weights))
    } else {
      proportions <- populations
      weights <- colSums(proportions)
    }
    num <- length(weights)

    list(proportions = proportions, totals = totals,
         weights = array(weights, dim = c(1, num)), num = num)
  }
