#' Calculate Similarity-Sensitive Shannon Entropy and Diversity
#'
#' @description
#' This function computes a generalized version of Shannon entropy and the corresponding
#' diversity measure that accounts for similarities between categories. It also calculates
#' the theoretical maximum values of both measures.
#'
#' @details
#' The similarity-sensitive Shannon entropy is defined as:
#' \deqn{H_S(p) = -\sum_{i=1}^{n} \sum_{j=1}^{n} s_{ij} p_i \log(p_j)}
#' The corresponding diversity measure is:
#' \deqn{D_S(p) = e^{H_S(p)}}
#' where p_i are the probabilities, s_{ij} is the similarity between categories i and j,
#' and n is the number of categories.
#'
#' @param probabilities A numeric vector of probabilities that should sum to 1
#' @param similarity_matrix A square matrix of similarities between categories.
#'        Diagonal elements are typically 1 (self-similarity).
#'        Off-diagonal elements typically range from 0 (no similarity) to 1 (identical).
#' @param base The logarithm base to use (default: natural logarithm, base e)
#' @param normalize_similarity Logical. Whether to normalize the similarity matrix
#'        so that rows sum to 1. Default is FALSE.
#' @param optimization_method Method to find maximum entropy:
#'        "analytical" - For standard Shannon entropy or simple similarity matrices
#'        "numerical" - Uses optimization for complex similarity matrices (default)
#' @param check_inputs Logical. Whether to validate inputs. Default is TRUE.
#' @param return_optimal_probs Logical. Whether to return the probability distribution
#'        that achieves maximum entropy/diversity. Default is TRUE.
#'
#' @return A list containing:
#'         - entropy: The similarity-sensitive Shannon entropy for the provided probabilities
#'         - diversity: The corresponding diversity measure (exp(entropy))
#'         - max_entropy: The maximum possible similarity-sensitive Shannon entropy
#'         - max_diversity: The maximum possible diversity measure
#'         - normalized_entropy: The ratio of entropy to max_entropy (ranges from 0 to 1)
#'         - normalized_diversity: The ratio of diversity to max_diversity (ranges from 0 to 1)
#'         - optimal_probs: (if return_optimal_probs=TRUE) The probability distribution
#'                         that achieves maximum entropy/diversity
#'
#' @examples
#' # Create a custom similarity matrix
#' custom_sim <- matrix(c(
#'   1.0, 0.3, 0.1,
#'   0.3, 1.0, 0.6,
#'   0.1, 0.6, 1.0
#' ), nrow = 3, byrow = TRUE)
#'
#' # Calculate entropy, diversity and their maximum values
#' probs <- c(0.5, 0.3, 0.2)
#' result <- similarity_sensitive_entropy(probs, custom_sim)
#' print(result)
#'
#' @export
similarity_sensitive_entropy <- function(probabilities, similarity_matrix,
                                         base = exp(1),
                                         normalize_similarity = FALSE,
                                         optimization_method = c("numerical", "analytical"),
                                         check_inputs = TRUE,
                                         return_optimal_probs = FALSE) {

  optimization_method <- match.arg(optimization_method)

  if (check_inputs) {
    # Validate probability vector
    if (!is.numeric(probabilities) || any(is.na(probabilities))) {
      stop("Probabilities must be numeric with no NA values")
    }

    if (any(probabilities < 0)) {
      stop("Probabilities cannot be negative")
    }

    if (!isTRUE(all.equal(sum(probabilities), 1, tolerance = 1e-8))) {
      warning("Probabilities do not sum to 1. Normalizing automatically.")
      probabilities <- probabilities / sum(probabilities)
    }

    # Validate similarity matrix
    if (!is.matrix(similarity_matrix)) {
      stop("Similarity matrix must be a matrix")
    }

    if (nrow(similarity_matrix) != ncol(similarity_matrix)) {
      stop("Similarity matrix must be square")
    }

    if (length(probabilities) != nrow(similarity_matrix)) {
      stop("Number of probabilities must match dimensions of similarity matrix")
    }

    if (any(is.na(similarity_matrix))) {
      stop("Similarity matrix contains NA values")
    }

    # Check if similarity matrix is valid (typically between 0 and 1)
    if (any(similarity_matrix < 0)) {
      warning("Similarity matrix contains negative values")
    }

    # Check diagonal elements (typically 1 for self-similarity)
    if (!all(diag(similarity_matrix) == 1)) {
      warning("Diagonal elements of similarity matrix are not all 1. This is unusual for self-similarity.")
    }
  }

  # Normalize similarity matrix if requested
  if (normalize_similarity) {
    row_sums <- rowSums(similarity_matrix)
    if (any(row_sums == 0)) {
      stop("Cannot normalize similarity matrix: some rows sum to zero")
    }
    similarity_matrix <- sweep(similarity_matrix, 1, row_sums, "/")
  }

  # Get number of categories
  n <- length(probabilities)

  #--------------------------------------------------------------------------
  # Calculate the similarity-sensitive Shannon entropy for given probabilities
  #--------------------------------------------------------------------------
  calculate_entropy <- function(probs, sim_matrix) {
    # Handle zero probabilities in log calculation to avoid -Inf
    # Using the convention that 0 * log(0) = 0
    log_probs <- rep(0, length(probs))
    non_zero <- probs > 0
    log_probs[non_zero] <- log(probs[non_zero], base = base)

    # Calculate using matrix operations
    p_outer <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      if (probs[i] > 0) {  # Skip when p_i = 0
        p_outer[i, ] <- sim_matrix[i, ] * probs[i]
      }
    }

    # Matrix multiply with log_probs to get the final sum
    -sum(p_outer * rep(log_probs, each = n))
  }

  # Calculate entropy for the provided probabilities
  entropy <- calculate_entropy(probabilities, similarity_matrix)

  # Calculate diversity as exp(entropy)
  diversity <- exp(entropy)

  #--------------------------------------------------------------------------
  # Find the maximum possible similarity-sensitive Shannon entropy
  #--------------------------------------------------------------------------

  # For standard Shannon entropy (identity similarity matrix),
  # the maximum is achieved with uniform distribution
  is_identity <- all(diag(similarity_matrix) == 1) &&
    all(similarity_matrix[lower.tri(similarity_matrix)] == 0) &&
    all(similarity_matrix[upper.tri(similarity_matrix)] == 0)

  # Check if similarity matrix is symmetric (important for optimization)
  is_symmetric <- isSymmetric(similarity_matrix, tol = 1e-8)
  if (!is_symmetric) {
    warning("Similarity matrix is not symmetric. This may affect the maximum entropy calculation.")
  }

  if (is_identity && optimization_method == "analytical") {
    # For identity matrix, maximum entropy is with uniform distribution
    optimal_probs <- rep(1/n, n)
    max_entropy <- calculate_entropy(optimal_probs, similarity_matrix)
  } else {
    # For complex similarity matrices, use numerical optimization

    # Objective function for optimization (negative entropy to minimize)
    objective <- function(probs) {
      # Ensure probabilities sum to 1
      probs <- probs / sum(probs)
      -calculate_entropy(probs, similarity_matrix)
    }

    # Initial guess: uniform distribution
    initial_guess <- rep(1/n, n)

    # Use constrOptim for constrained optimization
    require(stats)

    # Inequality constraints: all probabilities ≥ 0
    ui <- diag(n)  # Identity matrix for p_i ≥ 0
    ci <- rep(0, n)  # Lower bound of 0 for all p_i

    result <- try({
      constrOptim(
        theta = initial_guess,
        f = objective,
        grad = NULL,  # Use numerical gradients
        ui = ui,
        ci = ci,
        method = "BFGS"
      )
    }, silent = TRUE)

    if (inherits(result, "try-error")) {
      # Fall back to optim with projection
      projection <- function(probs) {
        probs[probs < 0] <- 0  # Ensure non-negativity
        probs / sum(probs)     # Ensure sum to 1
      }

      result <- optim(
        par = initial_guess,
        fn = function(p) {
          objective(projection(p))
        },
        method = "BFGS"
      )
      optimal_probs <- projection(result$par)
    } else {
      optimal_probs <- result$par / sum(result$par)  # Normalize to ensure sum=1
    }

    max_entropy <- calculate_entropy(optimal_probs, similarity_matrix)
  }

  # Calculate maximum diversity
  max_diversity <- exp(max_entropy)

  # Calculate normalized entropy and diversity
  if (max_entropy > 0) {
    normalized_entropy <- entropy / max_entropy
    normalized_diversity <- diversity / max_diversity
  } else {
    normalized_entropy <- NA
    normalized_diversity <- NA
    warning("Maximum entropy is zero or negative, unable to normalize")
  }

  # Prepare result
  result <- list(
    entropy = entropy,
    diversity = diversity,
    max_entropy = max_entropy,
    max_diversity = max_diversity,
    normalized_entropy = normalized_entropy,
    normalized_diversity = normalized_diversity
  )

  if (return_optimal_probs) {
    result$optimal_probs <- optimal_probs
  }

  return(result)
}

#' Calculate Similarity-Sensitive Cross-Entropy Between Probability Distributions
#'
#' @description
#' This function computes the similarity-sensitive cross-entropy between two
#' probability distributions, which generalizes the standard cross-entropy by
#' incorporating similarities between categories.
#'
#' @details
#' The similarity-sensitive cross-entropy between distributions P and Q is defined as:
#' \deqn{H_S(P,Q) = -\sum_{i=1}^{n} \sum_{j=1}^{n} s_{ij} p_i \log(q_j)}
#' where p_i and q_j are the probabilities in distributions P and Q,
#' s_{ij} is the similarity between categories i and j, and n is the number of categories.
#'
#' @param p A numeric vector representing the first probability distribution (should sum to 1)
#' @param q A numeric vector representing the second probability distribution (should sum to 1)
#' @param similarity_matrix A square matrix of similarities between categories.
#' @param base The logarithm base to use (default: natural logarithm, base e)
#' @param normalize_similarity Logical. Whether to normalize the similarity matrix
#'        so that rows sum to 1. Default is FALSE.
#'
#' @return The similarity-sensitive cross-entropy value
#'
#' @examples
#' # Create two probability distributions
#' p <- c(0.5, 0.3, 0.2)
#' q <- c(0.2, 0.3, 0.5)
#'
#' # Create a similarity matrix
#' sim_matrix <- matrix(c(
#'   1.0, 0.3, 0.1,
#'   0.3, 1.0, 0.6,
#'   0.1, 0.6, 1.0
#' ), nrow = 3, byrow = TRUE)
#'
#' # Calculate cross-entropy
#' cross_ent <- similarity_sensitive_cross_entropy(p, q, sim_matrix)
#' print(cross_ent)
#'
#' @export
similarity_sensitive_cross_entropy <- function(p, q, similarity_matrix,
                                               base = exp(1),
                                               normalize_similarity = FALSE) {

  # Validate inputs
  if (!is.numeric(p) || !is.numeric(q)) {
    stop("Probability distributions must be numeric vectors")
  }

  if (length(p) != length(q)) {
    stop("Probability distributions must have the same length")
  }

  if (!is.matrix(similarity_matrix)) {
    stop("Similarity matrix must be a matrix")
  }

  if (nrow(similarity_matrix) != length(p) || ncol(similarity_matrix) != length(p)) {
    stop("Dimensions of similarity matrix must match length of probability vectors")
  }

  # Normalize probability vectors if needed
  if (!isTRUE(all.equal(sum(p), 1, tolerance = 1e-8))) {
    warning("First probability vector does not sum to 1. Normalizing automatically.")
    p <- p / sum(p)
  }

  if (!isTRUE(all.equal(sum(q), 1, tolerance = 1e-8))) {
    warning("Second probability vector does not sum to 1. Normalizing automatically.")
    q <- q / sum(q)
  }

  # Normalize similarity matrix if requested
  if (normalize_similarity) {
    row_sums <- rowSums(similarity_matrix)
    if (any(row_sums == 0)) {
      stop("Cannot normalize similarity matrix: some rows sum to zero")
    }
    similarity_matrix <- sweep(similarity_matrix, 1, row_sums, "/")
  }

  # Handle zero probabilities in log calculation to avoid -Inf
  log_q <- rep(0, length(q))
  non_zero_q <- q > 0
  log_q[non_zero_q] <- log(q[non_zero_q], base = base)

  # Calculate cross-entropy using matrix operations
  n <- length(p)
  p_outer <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    if (p[i] > 0) {  # Skip when p_i = 0
      p_outer[i, ] <- similarity_matrix[i, ] * p[i]
    }
  }

  # If any q_j = 0 and the sum of p_i * s_ij > 0 for that j, cross-entropy is infinity
  for (j in 1:n) {
    if (q[j] == 0 && sum(p_outer[, j]) > 0) {
      return(Inf)
    }
  }

  # Matrix multiply with log_q to get the final sum
  cross_entropy <- -sum(p_outer * rep(log_q, each = n))

  return(cross_entropy)
}

#' Calculate Similarity-Sensitive Relative Entropy and Cross-Entropy for Multiple Distributions
#'
#' @description
#' This function computes similarity-sensitive relative entropy (Kullback-Leibler divergence)
#' and cross-entropy between multiple probability distributions organized as column vectors.
#'
#' @details
#' The similarity-sensitive relative entropy between distributions P and Q is defined as:
#' \deqn{D_{KL,S}(P||Q) = \sum_{i=1}^{n} \sum_{j=1}^{n} s_{ij} p_i \log(p_j/q_j)}
#' The similarity-sensitive cross-entropy is:
#' \deqn{H_S(P,Q) = -\sum_{i=1}^{n} \sum_{j=1}^{n} s_{ij} p_i \log(q_j)}
#'
#' @param probability_matrix A matrix where each column represents a probability distribution.
#'        Each column should sum to 1.
#' @param similarity_matrix A square matrix of similarities between categories.
#'        Dimensions should match the number of rows in the probability_matrix.
#' @param reference_distribution The column index of the distribution to use as the reference (Q).
#'        If NULL, all pairwise entropy measures are calculated.
#' @param base The logarithm base to use (default: natural logarithm, base e)
#' @param normalize_similarity Logical. Whether to normalize the similarity matrix
#'        so that rows sum to 1. Default is FALSE.
#' @param symmetric Logical. If TRUE, calculate a symmetric version of the divergence.
#'        Default is FALSE (standard asymmetric KL divergence).
#' @param return_cross_entropy Logical. Whether to return cross-entropy values along with
#'        relative entropy. Default is TRUE.
#' @param distribution_names Optional character vector of names for the distributions.
#'
#' @return A list containing:
#'         - relative_entropy: Matrix or vector of relative entropy values
#'         - cross_entropy: (if return_cross_entropy=TRUE) Matrix or vector of cross-entropy values
#'
#' @examples
#' # Create probability distributions as column vectors
#' probs_matrix <- cbind(
#'   c(0.5, 0.3, 0.2),  # Distribution 1
#'   c(0.2, 0.3, 0.5),  # Distribution 2
#'   c(1/3, 1/3, 1/3)   # Distribution 3 (uniform)
#' )
#'
#' # Create similarity matrix
#' sim_matrix <- matrix(c(
#'   1.0, 0.3, 0.1,
#'   0.3, 1.0, 0.6,
#'   0.1, 0.6, 1.0
#' ), nrow = 3)
#'
#' # Calculate relative entropy and cross-entropy
#' result <- similarity_sensitive_relative_entropy(probs_matrix, sim_matrix)
#' print(result)
#'
#' @export
similarity_sensitive_relative_entropy <- function(probability_matrix,
                                                  similarity_matrix,
                                                  reference_distribution = NULL,
                                                  base = exp(1),
                                                  normalize_similarity = FALSE,
                                                  symmetric = FALSE,
                                                  return_cross_entropy = TRUE,
                                                  distribution_names = NULL) {

  # Validate inputs
  if (!is.matrix(probability_matrix)) {
    stop("probability_matrix must be a matrix where each column is a probability distribution")
  }

  if (!is.matrix(similarity_matrix)) {
    stop("similarity_matrix must be a matrix")
  }

  n_categories <- nrow(probability_matrix)
  n_distributions <- ncol(probability_matrix)

  if (nrow(similarity_matrix) != n_categories || ncol(similarity_matrix) != n_categories) {
    stop("similarity_matrix dimensions must match the number of categories (rows) in probability_matrix")
  }

  # Check if columns sum to 1
  col_sums <- colSums(probability_matrix)
  if (!all(abs(col_sums - 1) < 1e-8)) {
    warning("Some columns in probability_matrix do not sum to 1. Normalizing automatically.")
    probability_matrix <- sweep(probability_matrix, 2, col_sums, "/")
  }

  # Normalize similarity matrix if requested
  if (normalize_similarity) {
    row_sums <- rowSums(similarity_matrix)
    if (any(row_sums == 0)) {
      stop("Cannot normalize similarity matrix: some rows sum to zero")
    }
    similarity_matrix <- sweep(similarity_matrix, 1, row_sums, "/")
  }

  # Initialize results matrices
  relative_entropy_matrix <- matrix(0, nrow = n_distributions, ncol = n_distributions)

  if (return_cross_entropy) {
    cross_entropy_matrix <- matrix(0, nrow = n_distributions, ncol = n_distributions)
  }

  # Function to calculate entropy and cross-entropy for distribution pair
  calculate_entropies <- function(p, q) {
    # Calculate self-entropy of p
    self_entropy <- similarity_sensitive_entropy(p, similarity_matrix,
                                                 base = base,
                                                 normalize_similarity = FALSE,
                                                 check_inputs = FALSE,
                                                 return_optimal_probs = FALSE)$entropy

    # Calculate cross-entropy between p and q
    cross_entropy <- similarity_sensitive_cross_entropy(p, q, similarity_matrix,
                                                        base = base,
                                                        normalize_similarity = FALSE)

    # Calculate relative entropy (KL divergence)
    relative_entropy <- cross_entropy - self_entropy

    return(list(
      relative_entropy = relative_entropy,
      cross_entropy = cross_entropy
    ))
  }

  # If reference_distribution is provided, calculate entropies relative to it
  if (!is.null(reference_distribution)) {
    if (reference_distribution < 1 || reference_distribution > n_distributions) {
      stop("reference_distribution index is out of bounds")
    }

    reference_dist <- probability_matrix[, reference_distribution]
    relative_entropy_vector <- numeric(n_distributions)

    if (return_cross_entropy) {
      cross_entropy_vector <- numeric(n_distributions)
    }

    # The entropy of reference compared to itself is 0
    relative_entropy_vector[reference_distribution] <- 0

    if (return_cross_entropy) {
      # Self-entropy for cross-entropy calculation
      self_entropy <- similarity_sensitive_entropy(reference_dist, similarity_matrix,
                                                   base = base,
                                                   normalize_similarity = FALSE,
                                                   check_inputs = FALSE,
                                                   return_optimal_probs = FALSE)$entropy
      cross_entropy_vector[reference_distribution] <- self_entropy
    }

    # Calculate for all other distributions
    for (i in setdiff(1:n_distributions, reference_distribution)) {
      p <- probability_matrix[, i]

      if (symmetric) {
        # For symmetric version, calculate both directions and average
        forward <- calculate_entropies(p, reference_dist)
        reverse <- calculate_entropies(reference_dist, p)

        relative_entropy_vector[i] <- (forward$relative_entropy + reverse$relative_entropy) / 2

        if (return_cross_entropy) {
          cross_entropy_vector[i] <- (forward$cross_entropy + reverse$cross_entropy) / 2
        }
      } else {
        # Standard asymmetric calculation
        result <- calculate_entropies(p, reference_dist)
        relative_entropy_vector[i] <- result$relative_entropy

        if (return_cross_entropy) {
          cross_entropy_vector[i] <- result$cross_entropy
        }
      }
    }

    # Add names if provided
    if (!is.null(distribution_names)) {
      if (length(distribution_names) == n_distributions) {
        names(relative_entropy_vector) <- distribution_names
        if (return_cross_entropy) {
          names(cross_entropy_vector) <- distribution_names
        }
      } else {
        warning("Length of distribution_names doesn't match number of distributions. Names not applied.")
      }
    }

    # Return results
    if (return_cross_entropy) {
      return(list(
        relative_entropy = relative_entropy_vector,
        cross_entropy = cross_entropy_vector
      ))
    } else {
      return(list(
        relative_entropy = relative_entropy_vector
      ))
    }
  } else {
    # Calculate all pairwise entropies
    for (i in 1:n_distributions) {
      for (j in 1:n_distributions) {
        if (i == j) {
          # The entropy of a distribution relative to itself is 0
          relative_entropy_matrix[i, j] <- 0

          if (return_cross_entropy) {
            # For cross-entropy, when i=j it's just the self-entropy
            p_i <- probability_matrix[, i]
            self_entropy <- similarity_sensitive_entropy(p_i, similarity_matrix,
                                                         base = base,
                                                         normalize_similarity = FALSE,
                                                         check_inputs = FALSE,
                                                         return_optimal_probs = FALSE)$entropy
            cross_entropy_matrix[i, j] <- self_entropy
          }
        } else {
          p_i <- probability_matrix[, i]
          p_j <- probability_matrix[, j]

          if (symmetric) {
            # For symmetric version, calculate both directions and average
            forward <- calculate_entropies(p_i, p_j)
            reverse <- calculate_entropies(p_j, p_i)

            relative_entropy_matrix[i, j] <- (forward$relative_entropy + reverse$relative_entropy) / 2

            if (return_cross_entropy) {
              cross_entropy_matrix[i, j] <- (forward$cross_entropy + reverse$cross_entropy) / 2
            }
          } else {
            # Standard asymmetric calculation
            result <- calculate_entropies(p_i, p_j)
            relative_entropy_matrix[i, j] <- result$relative_entropy

            if (return_cross_entropy) {
              cross_entropy_matrix[i, j] <- result$cross_entropy
            }
          }
        }
      }
    }

    # Add names if provided
    if (!is.null(distribution_names)) {
      if (length(distribution_names) == n_distributions) {
        dimnames(relative_entropy_matrix) <- list(distribution_names, distribution_names)
        if (return_cross_entropy) {
          dimnames(cross_entropy_matrix) <- list(distribution_names, distribution_names)
        }
      } else {
        warning("Length of distribution_names doesn't match number of distributions. Names not applied.")
      }
    }

    # Return results
    if (return_cross_entropy) {
      return(list(
        relative_entropy = relative_entropy_matrix,
        cross_entropy = cross_entropy_matrix
      ))
    } else {
      return(list(
        relative_entropy = relative_entropy_matrix
      ))
    }
  }
}



#' @description helper function for log calculation
 logf = function(x, base=log_base){
    log(x) / log(base)
  }
#'
#'
#'
