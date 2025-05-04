#' Calculate the weighted Gini-Simpson diversity index
#'
#' @description
#' The weighted Gini-Simpson index extends the standard Gini-Simpson index by
#' incorporating pairwise weights. It is defined as:
#' \deqn{WGS = \sum_{i=1}^{n} \sum_{j=1}^{n} w_{ij} p_i p_j}
#' where p_i is the probability of category i, and w_{ij} is the weight
#' between categories i and j.
#'
#' For the standard Gini-Simpson index, w_{ij} = 1 - \delta_{ij}, where
#' \delta_{ij} is the Kronecker delta (1 if i=j, 0 otherwise).
#'
#' @param probabilities A numeric vector of probabilities that should sum to 1
#' @param weights A symmetric matrix of weights between categories.
#'                If NULL, uses standard Gini-Simpson weights (1-δ_ij)
#' @param check_symmetry Logical. Whether to verify that the weight matrix is symmetric
#'
#' @return The weighted Gini-Simpson index value
#'
#' @examples
#' # Basic usage with standard weights
#' probabilities <- c(0.2, 0.3, 0.5)
#' weighted_gini_simpson_index(probabilities)
#'
#' # With custom weights
#' weights <- matrix(c(
#'   0.0, 0.5, 0.7,
#'   0.5, 0.0, 0.6,
#'   0.7, 0.6, 0.0
#' ), nrow = 3, byrow = TRUE)
#' weighted_gini_simpson_index(probabilities, weights)
#'
#' @export
weighted_gini_simpson_index <- function(probabilities, weights = NULL, check_symmetry = TRUE) {
  # Validate probabilities
  if (any(probabilities < 0)) {
    stop("Probabilities cannot be negative")
  }

  if (!isTRUE(all.equal(sum(probabilities), 1, tolerance = 1e-10))) {
    stop("Probabilities must sum to 1")
  }

  n <- length(probabilities)

  # If no weights provided, use standard Gini-Simpson weights
  if (is.null(weights)) {
    # Create weights matrix using 1-δ_ij (1 minus Kronecker delta)
    weights <- matrix(1, nrow = n, ncol = n)
    diag(weights) <- 0
  } else {
    # Validate provided weights
    if (!is.matrix(weights)) {
      stop("Weights must be a matrix")
    }

    if (nrow(weights) != n || ncol(weights) != n) {
      stop("Weights matrix dimensions must match length of probabilities vector")
    }

    # Check symmetry if requested
    if (check_symmetry && !isSymmetric(weights, tol = 1e-10)) {
      stop("Weights matrix must be symmetric")
    }
  }

  # Calculate weighted Gini-Simpson index
  # We compute the double summation: sum_i sum_j w_ij * p_i * p_j

  # Create a column vector of probabilities
  p_col <- matrix(probabilities, ncol = 1)

  # Calculate p_i * p_j for all i,j pairs
  p_outer <- p_col %*% t(p_col)

  # Element-wise multiply with weights and sum all elements
  result <- sum(weights * p_outer)

  return(result)

  # Alternative Method Using explicit loops
  # result <- 0
  # for (i in 1:n) {
  #   for (j in 1:n) {
  #     result <- result + weights[i, j] * probabilities[i] * probabilities[j]
  #   }
  # }
  # return(result)
}



#' Calculate the weighted Rich-Gini-Simpson diversity index
#'
#' @description
#' The weighted Rich-Gini-Simpson index combines the concepts of richness
#' (number of categories) and evenness (distribution of probabilities), while
#' also incorporating weights between categories. It is defined as:
#' \deqn{WRGS = n \times \sum_{i=1}^{n} \sum_{j=1}^{n} w_{ij} p_i p_j}
#' where n is the number of categories, p_i is the probability of category i,
#' and w_{ij} is the weight between categories i and j.
#'
#' For the standard Rich-Gini-Simpson index, w_{ij} = 1 - \delta_{ij}, where
#' \delta_{ij} is the Kronecker delta (1 if i=j, 0 otherwise).
#'
#' @param probabilities A numeric vector of probabilities that should sum to 1
#' @param weights A symmetric matrix of weights between categories.
#'                If NULL, uses standard Gini-Simpson weights (1-δ_ij)
#' @param check_symmetry Logical. Whether to verify that the weight matrix is symmetric
#' @param zero_handling How to handle categories with zero probability:
#'                     "keep" (default) - include zeros in richness count
#'                     "remove" - remove categories with zero probability
#'                     "ignore" - don't multiply by n (use regular weighted Gini-Simpson)
#'
#' @return The weighted Rich-Gini-Simpson index value
#'
#' @examples
#' # Basic usage
#' probabilities <- c(0.2, 0.3, 0.5)
#' weighted_rich_gini_simpson_index(probabilities)
#'
#' # With zero probabilities
#' probs_with_zeros <- c(0.3, 0.0, 0.2, 0.5)
#' weighted_rich_gini_simpson_index(probs_with_zeros)
#' weighted_rich_gini_simpson_index(probs_with_zeros, zero_handling = "remove")
#'
#' # With custom weights
#' weights <- matrix(c(
#'   0.0, 0.5, 0.7,
#'   0.5, 0.0, 0.6,
#'   0.7, 0.6, 0.0
#' ), nrow = 3, byrow = TRUE)
#' weighted_rich_gini_simpson_index(probabilities, weights)
#'
#' @export
weighted_rich_gini_simpson_index <- function(probabilities, weights = NULL,
                                             check_symmetry = TRUE,
                                             zero_handling = c("keep", "remove", "ignore"),
                                             normalize = TRUE) {
  # Process zero_handling parameter
  zero_handling <- match.arg(zero_handling)

  # Validate probabilities
  if (any(probabilities < 0)) {
    stop("Probabilities cannot be negative")
  }

  if (!isTRUE(all.equal(sum(probabilities), 1, tolerance = 1e-10))) {
    stop("Probabilities must sum to 1")
  }

  # Handle zeros if requested
  if (zero_handling == "remove") {
    non_zero <- probabilities > 0
    probabilities <- probabilities[non_zero]
    # Need to renormalize if we've removed values
    probabilities <- probabilities / sum(probabilities)
  }

  n <- length(probabilities)

  # If no weights provided, use standard Gini-Simpson weights
  if (is.null(weights)) {
    # Create weights matrix using 1-δ_ij (1 minus Kronecker delta)
    weights <- matrix(1, nrow = n, ncol = n)
    diag(weights) <- 0
  } else {
    # Validate provided weights
    if (!is.matrix(weights)) {
      stop("Weights must be a matrix")
    }

    if (zero_handling == "remove") {
      # Adjust weights matrix if we removed zero probabilities
      non_zero <- probabilities > 0
      weights <- weights[non_zero, non_zero]
    }

    if (nrow(weights) != n || ncol(weights) != n) {
      stop("Weights matrix dimensions must match length of probabilities vector")
    }

    # Check symmetry if requested
    if (check_symmetry && !isSymmetric(weights, tol = 1e-10)) {
      stop("Weights matrix must be symmetric")
    }
  }

  # Calculate weighted Gini-Simpson component
  # We compute the double summation: sum_i sum_j w_ij * p_i * p_j
  p_col <- matrix(probabilities, ncol = 1)
  p_outer <- p_col %*% t(p_col)

  ## weights are given weights * richness
  weights = weights * n


  gs_component <- sum(weights * p_col)

  if (normalize){
    maxweight <- max(weights)
    beta1 = maxweight * (1 - (1/n))
    ngs_component = gs_component / beta1
  }
  else {
    ngs_component = NULL
  }

  returnList = list(
    gs_component = gs_component,
    beta1 = beta1,
    ngs_component = ngs_component
  )

  return(returnList)
}
