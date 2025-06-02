#' #' Calculate Similarity-Sensitive Shannon Entropy
#' #'
#' #' @description
#' #' This function computes a generalized version of Shannon entropy that accounts
#' #' for similarities between categories. The classical Shannon entropy assumes all
#' #' categories are equally distinct, while this version adjusts the entropy based
#' #' on a provided similarity matrix.
#' #'
#' #' @details
#' #' The similarity-sensitive Shannon entropy is defined as:
#' #' \deqn{H_S(p) = -\sum_{i=1}^{n} \sum_{j=1}^{n} s_{ij} p_i \log(p_j)}
#' #' where p_i are the probabilities, s_{ij} is the similarity between categories i and j,
#' #' and n is the number of categories.
#' #'
#' #' When the similarity matrix is the identity matrix (1 on diagonal, 0 elsewhere),
#' #' this reduces to the classical Shannon entropy: \deqn{H(p) = -\sum_{i=1}^{n} p_i \log(p_i)}
#' #'
#' #' @param probabilities A numeric vector of probabilities that sum to 1
#' #' @param similarity_matrix A square matrix of similarities.
#' #' @param base The logarithm base to use (default: natural logarithm, base e)
#' #' @param normalize_similarity Logical. Whether to normalize the similarity matrix
#' #'        so that rows sum to 1. Default is FALSE.
#' #' @param check_inputs Logical. Whether to validate inputs. Default is TRUE.
#' #'
#' #' @return The similarity-sensitive Shannon entropy value
#' #'
#' #' @examples
#' #' # Calculate regular Shannon entropy (using identity similarity matrix)
#' #' probs <- c(0.5, 0.25, 0.25)
#' #' sim_matrix <- diag(3)
#' #' similarity_sensitive_shannon(probs, sim_matrix)
#' #'
#' #' # With a custom similarity matrix
#' #' custom_sim <- matrix(c(
#' #'   1.0, 0.3, 0.1,
#' #'   0.3, 1.0, 0.6,
#' #'   0.1, 0.6, 1.0
#' #' ), nrow = 3, byrow = TRUE)
#' #' similarity_sensitive_shannon(probs, custom_sim)
#' #'
#' #' @export
#' similarity_sensitive_shannon_ent <- function(probabilities, similarity_matrix,
#'                                          base = exp(1),
#'                                          normalize_similarity = FALSE,
#'                                          check_inputs = TRUE) {
#'
#'   if (check_inputs) {
#'     # Validate probability vector
#'     if (!is.numeric(probabilities) || any(is.na(probabilities))) {
#'       stop("Probabilities must be numeric with no NA values")
#'     }
#'
#'     if (any(probabilities < 0)) {
#'       stop("Probabilities cannot be negative")
#'     }
#'
#'     if (!isTRUE(all.equal(sum(probabilities), 1, tolerance = 1e-8))) {
#'       warning("Probabilities do not sum to 1. Normalizing automatically.")
#'       probabilities <- probabilities / sum(probabilities)
#'     }
#'
#'     # Validate similarity matrix
#'     if (!is.matrix(similarity_matrix)) {
#'       stop("Similarity matrix must be a matrix")
#'     }
#'
#'     if (nrow(similarity_matrix) != ncol(similarity_matrix)) {
#'       stop("Similarity matrix must be square")
#'     }
#'
#'     if (length(probabilities) != nrow(similarity_matrix)) {
#'       stop("Number of probabilities must match dimensions of similarity matrix")
#'     }
#'
#'     if (any(is.na(similarity_matrix))) {
#'       stop("Similarity matrix contains NA values")
#'     }
#'
#'     # Check if similarity matrix is valid (typically between 0 and 1)
#'     if (any(similarity_matrix < 0)) {
#'       warning("Similarity matrix contains negative values")
#'     }
#'
#'     # Check diagonal elements (typically 1 for self-similarity)
#'     if (!all(diag(similarity_matrix) == 1)) {
#'       warning("Diagonal elements of similarity matrix are not all 1. This is unusual for self-similarity.")
#'     }
#'   }
#'
#'   # Normalize similarity matrix if requested
#'   if (normalize_similarity) {
#'     row_sums <- rowSums(similarity_matrix)
#'     if (any(row_sums == 0)) {
#'       stop("Cannot normalize similarity matrix: some rows sum to zero")
#'     }
#'     similarity_matrix <- sweep(similarity_matrix, 1, row_sums, "/")
#'   }
#'
#'   # Get number of categories
#'   n <- length(probabilities)
#'
#'   # Handle zero probabilities in log calculation to avoid -Inf
#'   # use the convention that 0 * log(0) = 0
#'   log_probs <- rep(0, n)
#'   non_zero <- probabilities > 0
#'   log_probs[non_zero] <- log(probabilities[non_zero], base = base)
#'
#'   # Calculate the similarity-sensitive Shannon entropy
#'   entropy <- 0
#'
#'   # Using nested loops (more explicit but slower)
#'   # for (i in 1:n) {
#'   #   for (j in 1:n) {
#'   #     if (probabilities[i] > 0) {  # Handle 0 * log(0) = 0 convention
#'   #       entropy <- entropy - (similarity_matrix[i, j] * probabilities[i] * log_probs[j])
#'   #     }
#'   #   }
#'   # }
#'
#'   # Using matrix operations (more efficient)
#'   # Calculate p_i * s_{ij} for all i,j pairs
#'   p_outer <- matrix(0, nrow = n, ncol = n)
#'   for (i in 1:n) {
#'     if (probabilities[i] > 0) {  # Skip when p_i = 0
#'       p_outer[i, ] <- similarity_matrix[i, ] * probabilities[i]
#'     }
#'   }
#'
#'   # Matrix multiply with log_probs to get the final sum
#'   entropy <- -sum(p_outer * rep(log_probs, each = n))
#'
#'   return(entropy)
#' }
#'
#' #' @description helper function for log calculation
#' logf = function(x, base=2){
#'   log(x) / log(base)
#' }
#'
#'
#' #' Calculate Similarity-Sensitive Shannon Entropy with Maximum Value
#' #'
#' #' @description
#' #' This function computes a generalized version of Shannon entropy that accounts
#' #' for similarities between categories. It also calculates the theoretical maximum
#' #' entropy possible for the given similarity matrix.
#' #'
#' #' @details
#' #' The similarity-sensitive Shannon entropy is defined as:
#' #' \deqn{H_S(p) = -\sum_{i=1}^{n} \sum_{j=1}^{n} s_{ij} p_i \log(p_j)}
#' #' where p_i are the probabilities, s_{ij} is the similarity between categories i and j,
#' #' and n is the number of categories.
#' #'
#' #' Finding the maximum entropy requires optimization, as the probability distribution
#' #' that maximizes similarity-sensitive entropy depends on the structure of the
#' #' similarity matrix.
#' #'
#' #' @param probabilities A numeric vector of probabilities that should sum to 1
#' #' @param similarity_matrix A square matrix of similarities between categories.
#' #'        Diagonal elements are typically 1 (self-similarity).
#' #'        Off-diagonal elements typically range from 0 (no similarity) to 1 (identical).
#' #' @param base The logarithm base to use (default: natural logarithm, base e)
#' #' @param normalize_similarity Logical. Whether to normalize the similarity matrix
#' #'        so that rows sum to 1. Default is FALSE.
#' #' @param optimization_method Method to find maximum entropy:
#' #'        "analytical" - For standard Shannon entropy or simple similarity matrices
#' #'        "numerical" - Uses optimization for complex similarity matrices (default)
#' #' @param check_inputs Logical. Whether to validate inputs. Default is TRUE.
#' #' @param return_optimal_probs Logical. Whether to return the probability distribution
#' #'        that achieves maximum entropy. Default is TRUE.
#' #'
#' #' @return A list containing:
#' #'         - entropy: The similarity-sensitive Shannon entropy for the provided probabilities
#' #'         - max_entropy: The maximum possible similarity-sensitive Shannon entropy
#' #'         - normalized_entropy: The ratio of entropy to max_entropy (ranges from 0 to 1)
#' #'         - optimal_probs: (if return_optimal_probs=TRUE) The probability distribution
#' #'                         that achieves maximum entropy
#' #'
#' #' @examples
#' #' # Create a custom similarity matrix
#' #' custom_sim <- matrix(c(
#' #'   1.0, 0.3, 0.1,
#' #'   0.3, 1.0, 0.6,
#' #'   0.1, 0.6, 1.0
#' #' ), nrow = 3, byrow = TRUE)
#' #'
#' #' # Calculate entropy and maximum entropy
#' #' probs <- c(0.5, 0.3, 0.2)
#' #' result <- similarity_sensitive_shannon(probs, custom_sim)
#' #' print(result)
#' #'
#' #' @export
#' similarity_sensitive_shannon_opt <- function(probabilities, similarity_matrix,
#'                                          base = exp(1),
#'                                          normalize_similarity = FALSE,
#'                                          optimization_method = c("numerical", "analytical"),
#'                                          check_inputs = TRUE,
#'                                          return_optimal_probs = FALSE) {
#'
#'   optimization_method <- match.arg(optimization_method)
#'
#'   if (check_inputs) {
#'     # Validate probability vector
#'     if (!is.numeric(probabilities) || any(is.na(probabilities))) {
#'       stop("Probabilities must be numeric with no NA values")
#'     }
#'
#'     if (any(probabilities < 0)) {
#'       stop("Probabilities cannot be negative")
#'     }
#'
#'     if (!isTRUE(all.equal(sum(probabilities), 1, tolerance = 1e-8))) {
#'       warning("Probabilities do not sum to 1. Normalizing automatically.")
#'       probabilities <- probabilities / sum(probabilities)
#'     }
#'
#'     # Validate similarity matrix
#'     if (!is.matrix(similarity_matrix)) {
#'       stop("Similarity matrix must be a matrix")
#'     }
#'
#'     if (nrow(similarity_matrix) != ncol(similarity_matrix)) {
#'       stop("Similarity matrix must be square")
#'     }
#'
#'     if (length(probabilities) != nrow(similarity_matrix)) {
#'       stop("Number of probabilities must match dimensions of similarity matrix")
#'     }
#'
#'     if (any(is.na(similarity_matrix))) {
#'       stop("Similarity matrix contains NA values")
#'     }
#'
#'     # Check if similarity matrix is valid (typically between 0 and 1)
#'     if (any(similarity_matrix < 0)) {
#'       warning("Similarity matrix contains negative values")
#'     }
#'
#'     # Check diagonal elements (typically 1 for self-similarity)
#'     if (!all(diag(similarity_matrix) == 1)) {
#'       warning("Diagonal elements of similarity matrix are not all 1. This is unusual for self-similarity.")
#'     }
#'   }
#'
#'   # Normalize similarity matrix if requested
#'   if (normalize_similarity) {
#'     row_sums <- rowSums(similarity_matrix)
#'     if (any(row_sums == 0)) {
#'       stop("Cannot normalize similarity matrix: some rows sum to zero")
#'     }
#'     similarity_matrix <- sweep(similarity_matrix, 1, row_sums, "/")
#'   }
#'
#'   # Get number of categories
#'   n <- length(probabilities)
#'
#'   #--------------------------------------------------------------------------
#'   # Calculate the similarity-sensitive Shannon entropy for given probabilities
#'   #--------------------------------------------------------------------------
#'   calculate_entropy <- function(probs, sim_matrix) {
#'     # Handle zero probabilities in log calculation to avoid -Inf
#'     # Using the convention that 0 * log(0) = 0
#'     log_probs <- rep(0, length(probs))
#'     non_zero <- probs > 0
#'     log_probs[non_zero] <- log(probs[non_zero], base = base)
#'
#'     # Calculate using matrix operations
#'     p_outer <- matrix(0, nrow = n, ncol = n)
#'     for (i in 1:n) {
#'       if (probs[i] > 0) {  # Skip when p_i = 0
#'         p_outer[i, ] <- sim_matrix[i, ] * probs[i]
#'       }
#'     }
#'
#'     # Matrix multiply with log_probs to get the final sum
#'     -sum(p_outer * rep(log_probs, each = n))
#'   }
#'
#'   # Calculate entropy for the provided probabilities
#'   entropy <- calculate_entropy(probabilities, similarity_matrix)
#'
#'   #--------------------------------------------------------------------------
#'   # Find the maximum possible similarity-sensitive Shannon entropy
#'   #--------------------------------------------------------------------------
#'
#'   # For standard Shannon entropy (identity similarity matrix),
#'   # the maximum is achieved with uniform distribution
#'   is_identity <- all(diag(similarity_matrix) == 1) &&
#'     all(similarity_matrix[lower.tri(similarity_matrix)] == 0) &&
#'     all(similarity_matrix[upper.tri(similarity_matrix)] == 0)
#'
#'   # Check if similarity matrix is symmetric (important for optimization)
#'   is_symmetric <- isSymmetric(similarity_matrix, tol = 1e-8)
#'   if (!is_symmetric) {
#'     warning("Similarity matrix is not symmetric. This may affect the maximum entropy calculation.")
#'   }
#'
#'   if (is_identity && optimization_method == "analytical") {
#'     # For identity matrix, maximum entropy is with uniform distribution
#'     optimal_probs <- rep(1/n, n)
#'     max_entropy <- calculate_entropy(optimal_probs, similarity_matrix)
#'   } else {
#'     # For complex similarity matrices, use numerical optimization
#'
#'     # Objective function for optimization (negative entropy to minimize)
#'     objective <- function(probs) {
#'       # Ensure probabilities sum to 1
#'       probs <- probs / sum(probs)
#'       -calculate_entropy(probs, similarity_matrix)
#'     }
#'
#'     # Constraint: probabilities sum to 1
#'     sum_constraint <- function(probs) {
#'       sum(probs) - 1
#'     }
#'
#'     # Initial guess: uniform distribution
#'     initial_guess <- rep(1/n, n)
#'
#'     # Use constrOptim for constrained optimization
#'     require(stats)
#'
#'     # Inequality constraints: all probabilities ≥ 0
#'     ui <- diag(n)  # Identity matrix for p_i ≥ 0
#'     ci <- rep(0, n)  # Lower bound of 0 for all p_i
#'
#'     result <- try({
#'       constrOptim(
#'         theta = initial_guess,
#'         f = objective,
#'         grad = NULL,  # Use numerical gradients
#'         ui = ui,
#'         ci = ci,
#'         method = "BFGS"
#'       )
#'     }, silent = TRUE)
#'
#'     if (inherits(result, "try-error")) {
#'       # Fall back to optim with projection
#'       projection <- function(probs) {
#'         probs[probs < 0] <- 0  # Ensure non-negativity
#'         probs / sum(probs)     # Ensure sum to 1
#'       }
#'
#'       result <- optim(
#'         par = initial_guess,
#'         fn = function(p) {
#'           objective(projection(p))
#'         },
#'         method = "BFGS"
#'       )
#'       optimal_probs <- projection(result$par)
#'     } else {
#'       optimal_probs <- result$par / sum(result$par)  # Normalize to ensure sum=1
#'     }
#'
#'     max_entropy <- calculate_entropy(optimal_probs, similarity_matrix)
#'   }
#'
#'   # Calculate normalized entropy (ratio of entropy to maximum entropy)
#'   if (max_entropy > 0) {
#'     normalized_entropy <- entropy / max_entropy
#'   } else {
#'     normalized_entropy <- NA
#'     warning("Maximum entropy is zero or negative, unable to normalize entropy")
#'   }
#'
#'   # Prepare result
#'   result <- list(
#'     entropy = entropy,
#'     max_entropy = max_entropy,
#'     normalized_entropy = normalized_entropy
#'   )
#'
#'   if (return_optimal_probs) {
#'     result$optimal_probs <- optimal_probs
#'   }
#'
#'   return(result)
#' }
