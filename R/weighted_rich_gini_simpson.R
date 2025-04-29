#' Calculate the weighted Rich-Gini-Simpson diversity index with affinity weights
#'
#' @description
#' This function computes the weighted Rich-Gini-Simpson index for seqs
#' abundance data collected across multiple phenos, while incorporating seqs
#' affinity importance. It also calculates the theoretical maximum value
#' of the index for comparison.
#'
#' @param abundance_matrix A matrix where rows represent seqs and columns
#'                        represent phenos. Values are seqs abundances or proportions.
#' @param sequence_weights A symmetric matrix of weights between seqs.
#'                                   If NULL, uses standard weights (1-δ_ij).
#' @param affinity_weights A vector of affinity importance weights for each seqs.
#'                            If NULL, all seqs are weighted equally.
#' @param normalize_phenos Logical. Whether to normalize abundances at each pheno
#'                       to sum to 1. Default is TRUE.
#' @param pheno_weights Optional vector of weights for each pheno. If NULL, all
#'                     phenos are weighted equally.
#' @param zero_handling How to handle seqs with zero abundance:
#'                     "keep" (default) - include zeros in richness count
#'                     "remove" - remove seqs with zero abundance at each pheno
#'                     "ignore" - don't multiply by n (use regular weighted Gini-Simpson)
#' @param return_pheno_values Logical. If TRUE, returns individual pheno values
#'                          along with the overall index. Default is FALSE.
#'
#' @return A list containing:
#'         - index: The overall weighted Rich-Gini-Simpson index
#'         - max_index: The theoretical maximum value of the index
#'         - normalized_index: The index divided by its maximum value (0-1 scale)
#'         - affinity_weighted_index: The index weighted by seqs affinity importance
#'         - pheno_indices: (if return_pheno_values=TRUE) Individual pheno indices
#'         - pheno_max_indices: (if return_pheno_values=TRUE) Maximum possible indices per pheno
#'         - sequence_weighted_pheno_indices: (if return_pheno_values=TRUE) affinity-weighted indices per pheno
#'
#' @examples
#' # Create a sample abundance matrix (5 seqs across 3 phenos)
#' abundance <- matrix(c(
#'   10, 5,  0,   # seqs 1 abundances across 3 phenos
#'   2,  15, 8,   # seqs 2
#'   0,  3,  12,  # seqs 3
#'   20, 10, 5,   # seqs 4
#'   5,  0,  2    # seqs 5
#' ), nrow = 5, byrow = TRUE)
#'
#' # Create affinity weights (higher values = more affinity importance)
#' affinity_weights <- c(2.0, 1.0, 3.0, 0.5, 1.5)
#'
#' # Get index with affinity weighting
#' result <- weighted_rich_gini_simpson_affinity(
#'   abundance,
#'   affinity_weights = affinity_weights
#' )
#' print(result)
#'
#' @export
weighted_rich_gini_simpson_affinity_pairs <- function(abundance_matrix,
                                                    sequence_weights = NULL,
                                                    affinity_weights = NULL,
                                                    #useRichness= TRUE,
                                                    normalize_phenos = TRUE,
                                                    pheno_weights = NULL,
                                                    zero_handling = c("keep", "remove", "ignore"),
                                                    return_pheno_values = FALSE) {
  # Process zero_handling parameter
  zero_handling <- match.arg(zero_handling)

  # Validate input matrix
  if (!is.matrix(abundance_matrix)) {
    stop("Abundance data must be provided as a matrix")
  }

  n_seqs <- nrow(abundance_matrix)
  n_phenos <- ncol(abundance_matrix)

  if (n_seqs == 0 || n_phenos == 0) {
    stop("Abundance matrix cannot have zero rows or columns")
  }

  # Check for negative values
  if (any(abundance_matrix < 0)) {
    stop("Abundance values cannot be negative")
  }

  # Initialize affinity weights if not provided
  if (is.null(affinity_weights)) {
    affinity_weights <- rep(1, n_seqs)
  } else {
    if (length(affinity_weights) != n_seqs) {
      stop("Length of affinity_weights must match the number of seqs (rows) in the abundance matrix")
    }
    if (any(affinity_weights < 0)) {
      stop("affinity weights cannot be negative")
    }
  }

  # Initialize pheno weights if not provided
  if (is.null(pheno_weights)) {
    pheno_weights <- rep(1/n_phenos, n_phenos)
  } else {
    if (length(pheno_weights) != n_phenos) {
      stop("Length of pheno_weights must match the number of phenos (columns) in the abundance matrix")
    }
    if (any(pheno_weights < 0)) {
      stop("pheno weights cannot be negative")
    }
    # Normalize pheno weights to sum to n_phenos
    #pheno_weights <- pheno_weights * n_phenos / sum(pheno_weights)
    # # Normalize pheno weights to sum to 1
    #pheno_weights <- pheno_weights  / sum(pheno_weights)

  }

  # Initialize seqs relationship weights if not provided
  if (is.null(sequence_weights)) {
    sequence_weights <- matrix(1, nrow = n_seqs, ncol = n_seqs)
    diag(sequence_weights) <- 0
  } else {
    if (!is.matrix(sequence_weights) ||
        nrow(sequence_weights) != n_seqs ||
        ncol(sequence_weights) != n_seqs) {
      stop("seqs relationship weights matrix dimensions must match the number of seqs")
    }
    if (!isSymmetric(sequence_weights, tol = 1e-10)) {
      stop("seqs relationship weights matrix must be symmetric")
    }
  }
  nullsequence_weights <- matrix(1, nrow = n_seqs, ncol = n_seqs)
  diag(nullsequence_weights) <- 0
  # Initialize storage for pheno-specific indices and their maximum values
  gsw_indices <- numeric(n_phenos)
  rgsw_indices <- numeric(n_phenos)
  rgs_indices <- numeric(n_phenos)

  vgsw_indices <- numeric(n_phenos)
  vrgsw_indices <- numeric(n_phenos)
  vrgs_indices <- numeric(n_phenos)

  gsw_weight_indices <- numeric(n_phenos)
  rgsw_weight_indices <- numeric(n_phenos)
  rgs_weight_indices <- numeric(n_phenos)
  ncells_indices <- numeric(n_phenos)

  pheno_max_indices <- numeric(n_phenos)
  sequence_weighted_pheno_indices <- numeric(n_phenos)
  pheno_probs <- list()
  weighted_pheno_probs <- list()


  pheno_names = colnames(abundance_matrix)
  ## convert abundance to probabilities
  probability_mat =  normalize_columns(abundance_matrix, handling_zero_cols = "zero", na.rm = FALSE)
  #probability_mat = abundance_matrix / sum(abundance_matrix)
  ##
  # Calculate index for each pheno
  for (pheno in 1:n_phenos) {
    # Extract abundances for this pheno
    ncells_indices[pheno] <- sum(abundance_matrix[, pheno])
    probs <- probability_mat[, pheno]

    # Handle zeros if necessary
    if (zero_handling == "remove") {
      non_zero <- probs > 0
      if (sum(non_zero) == 0) {
        pheno_indices[pheno] <- 0
        pheno_max_indices[pheno] <- 0
        sequence_weighted_pheno_indices[pheno] <- 0
        next
      }
      pheno_abundances <- probs[non_zero]
      pheno_rel_weights <- sequence_weights[non_zero, non_zero]
      pheno_cons_weights <- affinity_weights[non_zero]
    } else {
      pheno_abundances <- probs
      pheno_rel_weights <- sequence_weights
      pheno_cons_weights <- affinity_weights
    }

    # Normalize abundances to probabilities if requested
    # if (normalize_phenos) {
    #   total_abundance <- sum(pheno_abundances)
    #   if (total_abundance > 0) {
    #     probabilities <- pheno_abundances / total_abundance
    #   } else {
    #     pheno_indices[pheno] <- 0
    #     pheno_max_indices[pheno] <- 0
    #     sequence_weighted_pheno_indices[pheno] <- 0
    #     next
    #   }
    # } else {
    probabilities <- probability_mat[, pheno]
    #probabilities <- pheno_abundances
    #  }

    # Count seqs richness for this pheno
    if (zero_handling == "remove") {
      richness <- sum(probabilities > 0)
    } else {
      richness <- length(probabilities)
    }

    richness <- sum(probabilities > 0)
    ## calculate weights
    ##
    # Calculate Gini-Simpson component
    p_col <- matrix(probabilities, ncol = 1)
    p_outer <- p_col %*% t(p_col)
    pheno_probs[[pheno]] <- p_outer
    weighted_pheno_probs[[pheno]] <- p_outer * pheno_weights[pheno]

    gs_component <- sum(p_outer * sequence_weights * (1 - p_outer))
    pheno_indices[pheno] <- gs_component
    vweights = outer(as.vector(pheno_cons_weights), as.vector(pheno_cons_weights), FUN = "+")
    nweights = (richness * (richness - 1) / 2)

    w_gs_component <- upper_triangle_sum( sequence_weights * vweights * p_outer * (1-p_outer)) * nweights

    # Calculate final affinity-weighted index for this pheno
    #if (zero_handling == "ignore") {
    sequence_weighted_pheno_indices[pheno] <- w_gs_component
    #} else

    # # Count seqs richness for this pheno
    # if (zero_handling == "remove") {
    #   richness <- sum(probabilities > 0)
    # } else {
    #   richness <- length(probabilities)
    # }

    prichness <- sum(probabilities > 0)


    ## calculate weights
    ##
    # Calculate Gini-Simpson component
    p_col <- matrix(probabilities, ncol = 1)
    p_outer <- p_col %*% t(p_col)
    pheno_probs[[pheno]] <- p_outer
    weighted_pheno_probs[[pheno]] <- p_outer * pheno_weights[pheno]

    gsw_component <- sum(p_outer * sequence_weights * (1 - p_outer))
    gsw_indices[pheno] <- gsw_component


    # #} else {
    #  # pheno_indices[pheno] <- richness * gs_component
    # #}
    #
    # #------------------------------------------------------------------
    # # Calculate the maximum possible value of the index for this pheno
    # #------------------------------------------------------------------
    #
    # # For the maximum value, we assume equal probabilities
    # if (richness <= 1) {
    #   # With only 1 seqs, the index is always 0
    #   max_val <- 0
    # } else {
    #   # Calculate maximum value with equal probabilities
    #   equal_probs <- rep(1/richness, richness)
    #   p_col <- matrix(equal_probs, ncol = 1)
    #   p_outer <- p_col %*% t(p_col)
    #
    #   # The maximum GS component with equal probabilities
    #   max_gs_component <- sum(pheno_rel_weights * p_outer)
    #
    #   if (zero_handling == "ignore") {
    #     max_val <- max_gs_component
    #   } else {
    #     max_val <- richness * max_gs_component
    #   }
    # }
    #
    # pheno_max_indices[pheno] <- max_val

    #------------------------------------------------------------------
    # Calculate affinity-weighted index for this pheno
    #------------------------------------------------------------------

    vweights = outer(as.vector(pheno_cons_weights), as.vector(pheno_cons_weights), FUN = "+")
    nweights = (richness * (richness - 1) / 2)

    w_gs_component <- upper_triangle_sum( sequence_weights * vweights * p_outer * (1-p_outer)) * nweights

    vrgsw_component <- upper_triangle_sum( sequence_weights * vweights * p_outer * (1-p_outer)) * nweights

    ## no affinity weights
    rgsw_component <- sum( sequence_weights  * p_outer * (1-p_outer)) * nweights

    # Calculate final affinity-weighted index for this pheno
    #if (zero_handling == "ignore") {
    rgsw_indices[pheno] <- rgsw_component
    vrgsw_indices[pheno] <- vrgsw_component
    #} else {
    # sequence_weighted_pheno_indices[pheno] <- richness * cons_weighted_gs_component
    #}
    vrgs_component <- upper_triangle_sum(nullsequence_weights * vweights * p_outer * (1-p_outer)) * nweights
    vrgs_indices[pheno] <- vrgs_component

    rgs_component <- sum(  p_outer * (1-p_outer)) * nweights
    rgs_indices[pheno] <- rgs_component


    # Calculate final affinity-weighted index for this pheno
    #if (zero_handling == "ignore") {
    #rgs_indices[pheno] <- rgs_component

    #rgs_weight_indices[pheno]


  }

  richness <- nrow(probability_mat)
  ## calcuate beta
  nweights = (richness * (richness - 1) / 2)
  vweights = outer(as.vector(pheno_cons_weights), as.vector(pheno_cons_weights), FUN = "+")

  wij = sequence_weights * vweights #* nweights

  w_ilj = upper_triangle_sum(wij, diag_include = FALSE) #* nweights
  b = ((richness * (richness - 1) / 2) -1)^2

  #b = (nweights - 1)^2

  bc = 1 / (upper_triangle_sum( (1/wij), diag_include = FALSE))
  pbeta =  (w_ilj - (b * bc))

  vwij = sequence_weights * vweights #* nweights

  vw_ilj = upper_triangle_sum(vwij, diag_include = FALSE) #* nweights
  b = ((richness * (richness - 1) / 2) -1)^2

  #b = (nweights - 1)^2

  #### revised beta
  b = ((richness^2) -2)^2

  bc = 1 / (sum( (1/vwij)))

  beta =  (sum(vwij) - (b * bc)) * (1 / 4)


  ## beta with affinity weights without sequence weights
  ##
  nwij = nullsequence_weights * vweights #* nwe
  nw_ilj = upper_triangle_sum(nwij, diag_include = FALSE) #* nweights

  nbc = 1 / (sum( (1/nwij)))
  nbeta =  (sum(nwij) - (b * nbc)) * (1 / 4)


  nullaffinity_weights <- rep(1, n_seqs)
  nullvweights = outer(as.vector(nullaffinity_weights), as.vector(nullaffinity_weights), FUN = "+")
  wij = sequence_weights * nullvweights #* nwe

  bc = 1 / (sum( (1/wij)))
  nvbeta =  (sum(wij) - (b * bc)) * (1 / 4)




  #b = (nweights - 1)^2



  #calculate overall
  #
  #######
  #########
  #########
  #########
  #########
  #te the maximum possible value of the index for this pheno

  tprob <- Reduce('+', weighted_pheno_probs)
  p_outer <- tprob %*% t(tprob)



  gs_component <- sum(tprob * sequence_weights * (1 - tprob))

  vweights = outer(as.vector(pheno_cons_weights), as.vector(pheno_cons_weights), FUN = "+")
  nweights = (richness * (richness - 1) / 2)

  gamme_div <- upper_triangle_sum( sequence_weights * vweights * tprob * (1-tprob)) * nweights

  # pheno_indices[pheno] <- gs_component


  alpha_div = sum(pheno_weights * sequence_weighted_pheno_indices)

  vrgs_norm = vrgs_indices / nbeta
  vrgsw_norm = vrgsw_indices / beta
  rgsw_norm = rgsw_indices / nvbeta

  resdf = data.frame(
    pheno = pheno_names[ 1:n_phenos],
    rgs = rgs_indices,
    vrgs_norm = vrgs_norm,
    rgsw_norm = rgsw_norm,
    vrgsw_norm = vrgsw_norm,
    weighted_gini_simpson = sequence_weighted_pheno_indices ,
    ncells = ncells_indices,
    total_cells = sum(abundance_matrix))

  resdf$weighted_gini_simpsonn = resdf$weighted_gini_simpson / pbeta
  resdf$alpha_div = alpha_div
  resdf$gamme_div = gamme_div
  resdf$beta = pbeta


  # Calculate overall weighted indices across phenos
  #overall_index <- sum(pheno_indices * pheno_weights) / sum(pheno_weights)
  #overall_max_index <- sum(pheno_max_indices * pheno_weights) / sum(pheno_weights)
  #affinity_weighted_index <- sum(sequence_weighted_pheno_indices * pheno_weights) / sum(pheno_weights)

  # Calculate overall normalized index
  # if (overall_max_index > 0) {
  #   overall_normalized_index <- overall_index / overall_max_index
  # } else {
  #   overall_normalized_index <- NA
  # }
  #
  # Return results
  # result <- list(
  #   index = overall_index,
  #   max_index = overall_max_index,
  #   normalized_index = overall_normalized_index,
  #   affinity_weighted_index = affinity_weighted_index
  # )
  #
  # if (return_pheno_values) {
  #   result$pheno_indices <- pheno_indices
  #   result$pheno_max_indices <- pheno_max_indices
  #   result$sequence_weighted_pheno_indices <- sequence_weighted_pheno_indices
  #   result$pheno_weights <- pheno_weights
  # }

  return(resdf)
}


#' Calculate the weighted Rich-Gini-Simpson diversity index with affinity weights
#'
#' @description
#' This function computes the weighted Rich-Gini-Simpson index for seqs
#' abundance data collected across multiple phenos, while incorporating seqs
#' affinity importance. It also calculates the theoretical maximum value
#' of the index for comparison.
#'
#' @param abundance_matrix A matrix where rows represent seqs and columns
#'                        represent phenos. Values are seqs abundances or proportions.
#' @param sequence_weights A symmetric matrix of weights between seqs.
#'                                   If NULL, uses standard weights (1-δ_ij).
#' @param affinity_weights A vector of affinity importance weights for each seqs.
#'                            If NULL, all seqs are weighted equally.
#' @param normalize_phenos Logical. Whether to normalize abundances at each pheno
#'                       to sum to 1. Default is TRUE.
#' @param pheno_weights Optional vector of weights for each pheno. If NULL, all
#'                     phenos are weighted equally.
#' @param zero_handling How to handle seqs with zero abundance:
#'                     "keep" (default) - include zeros in richness count
#'                     "remove" - remove seqs with zero abundance at each pheno
#'                     "ignore" - don't multiply by n (use regular weighted Gini-Simpson)
#' @param return_pheno_values Logical. If TRUE, returns individual pheno values
#'                          along with the overall index. Default is FALSE.
#'
#' @return A list containing:
#'         - index: The overall weighted Rich-Gini-Simpson index
#'         - max_index: The theoretical maximum value of the index
#'         - normalized_index: The index divided by its maximum value (0-1 scale)
#'         - affinity_weighted_index: The index weighted by seqs affinity importance
#'         - pheno_indices: (if return_pheno_values=TRUE) Individual pheno indices
#'         - pheno_max_indices: (if return_pheno_values=TRUE) Maximum possible indices per pheno
#'         - sequence_weighted_pheno_indices: (if return_pheno_values=TRUE) affinity-weighted indices per pheno
#'
#' @examples
#' # Create a sample abundance matrix (5 seqs across 3 phenos)
#' abundance <- matrix(c(
#'   10, 5,  0,   # seqs 1 abundances across 3 phenos
#'   2,  15, 8,   # seqs 2
#'   0,  3,  12,  # seqs 3
#'   20, 10, 5,   # seqs 4
#'   5,  0,  2    # seqs 5
#' ), nrow = 5, byrow = TRUE)
#'
#' # Create affinity weights (higher values = more affinity importance)
#' affinity_weights <- c(2.0, 1.0, 3.0, 0.5, 1.5)
#'
#' # Get index with affinity weighting
#' result <- weighted_rich_gini_simpson_affinity(
#'   abundance,
#'   affinity_weights = affinity_weights
#' )
#' print(result)
#'
#' @export
weighted_rich_gini_simpson_affinity <- function(abundance_matrix,
                                                          affinity_weights = NULL,
                                                          #useRichness= TRUE,
                                                          normalize_phenos = TRUE,
                                                          zero_handling = c("keep", "remove", "ignore"),
                                                          return_pheno_values = FALSE) {
  # Process zero_handling parameter
  zero_handling <- match.arg(zero_handling)

  # Validate input matrix
  if (!is.matrix(abundance_matrix)) {
    stop("Abundance data must be provided as a matrix")
  }

  n_seqs <- nrow(abundance_matrix)
  n_phenos <- ncol(abundance_matrix)

  if (n_seqs == 0 || n_phenos == 0) {
    stop("Abundance matrix cannot have zero rows or columns")
  }

  # Check for negative values
  if (any(abundance_matrix < 0)) {
    stop("Abundance values cannot be negative")
  }

  # Initialize affinity weights if not provided
  if (is.null(affinity_weights)) {
    affinity_weights <- rep(1, n_seqs)
  } else {
    if (length(affinity_weights) != n_seqs) {
      stop("Length of affinity_weights must match the number of seqs (rows) in the abundance matrix")
    }
    if (any(affinity_weights < 0)) {
      stop("affinity weights cannot be negative")
    }
  }


    # Normalize pheno weights to sum to n_phenos
    #pheno_weights <- pheno_weights * n_phenos / sum(pheno_weights)
    # # Normalize pheno weights to sum to 1
    #pheno_weights <- pheno_weights  / sum(pheno_weights)




  # Initialize storage for pheno-specific indices and their maximum values
  pheno_indices <- numeric(n_phenos)
  weight_indices = numeric(n_phenos)
  pheno_max_indices <- numeric(n_phenos)
  gs_pheno_indices <- numeric(n_phenos)
  rgs_pheno_indices <- numeric(n_phenos)
  wrgs_pheno_indices <- numeric(n_phenos)
  sequence_weighted_pheno_indices <- numeric(n_phenos)
  pheno_probs <- list()
  weighted_pheno_probs <- list()


  pheno_names = colnames(abundance_matrix)
  ## convert abundance to probabilities
  probability_mat =  normalize_columns(abundance_matrix, handling_zero_cols = "zero", na.rm = FALSE)
  #probability_mat = abundance_matrix / sum(abundance_matrix)
  ##
  # Calculate index for each pheno
  for (pheno in 1:n_phenos) {
    # Extract abundances for this pheno
    probs <- probability_mat[, pheno]

    # Handle zeros if necessary
    if (zero_handling == "remove") {
      non_zero <- probs > 0
      if (sum(non_zero) == 0) {
        pheno_indices[pheno] <- 0
        pheno_max_indices[pheno] <- 0
        sequence_weighted_pheno_indices[pheno] <- 0
        next
      }
      pheno_abundances <- probs[non_zero]
      pheno_rel_weights <- sequence_weights[non_zero, non_zero]
      pheno_cons_weights <- affinity_weights[non_zero]
    } else {
      pheno_abundances <- probs
      pheno_cons_weights <- affinity_weights
    }

    # Normalize abundances to probabilities if requested
    # if (normalize_phenos) {
    #   total_abundance <- sum(pheno_abundances)
    #   if (total_abundance > 0) {
    #     probabilities <- pheno_abundances / total_abundance
    #   } else {
    #     pheno_indices[pheno] <- 0
    #     pheno_max_indices[pheno] <- 0
    #     sequence_weighted_pheno_indices[pheno] <- 0
    #     next
    #   }
    # } else {
    probabilities <- probability_mat[, pheno]
    #probabilities <- pheno_abundances
    #  }

    # Count seqs richness for this pheno
    if (zero_handling == "remove") {
      richness <- sum(probabilities > 0)
    } else {
      richness <- length(probabilities)
    }


    ## calculate weights
    ##
    # Calculate Gini-Simpson component
    p_col <- matrix(probabilities, ncol = 1)
    pheno_probs[[pheno]] <- p_col
    weighted_pheno_probs[[pheno]] <- p_col

    ## GS
    gs_component <- sum(p_col * (1 - p_col))
    # RGS component
    rgs_component <- sum(p_col * richness * (1 - p_col))

    weights = affinity_weights * richness
    #weight_indices[pheno] <- weights
      # weighted RGS component
    wrgs_component <- sum(p_col * weights * (1 - p_col))

    gs_pheno_indices[pheno] <- gs_component
    rgs_pheno_indices[pheno] <- rgs_component
    wrgs_pheno_indices[pheno] <- wrgs_component
    #}
  }

  ## calcuate beta
  ##
  maxw =  max(weights)

  beta = maxw * ( 1 - (1 / richness))



  #calculate overall
  #
  #######
  #########
  #########
  #########
  #########
  #te the maximum possible value of the index for this pheno

  tprob <- Reduce('+', weighted_pheno_probs)


  resdf = data.frame(
    pheno = pheno_names[ 1:n_phenos],
    gini_simpson = gs_pheno_indices,
    rich_gini_simpson = rgs_pheno_indices,
    weighted_rich_gini_simpson = wrgs_pheno_indices,
    weighted_rich_gini_simpsonn = wrgs_pheno_indices / beta)


  resdf$beta = beta
  #affinity_weighted_index = sequence_weighted_pheno_indices,
  #beta = beta,
#  resdf$alpha_div = alpha_div
#  resdf$gamme_div = gamme_div


  # Calculate overall weighted indices across phenos
  #overall_index <- sum(pheno_indices * pheno_weights) / sum(pheno_weights)
  #overall_max_index <- sum(pheno_max_indices * pheno_weights) / sum(pheno_weights)
  #affinity_weighted_index <- sum(sequence_weighted_pheno_indices * pheno_weights) / sum(pheno_weights)

  # Calculate overall normalized index
  # if (overall_max_index > 0) {
  #   overall_normalized_index <- overall_index / overall_max_index
  # } else {
  #   overall_normalized_index <- NA
  # }
  #
  # Return results
  # result <- list(
  #   index = overall_index,
  #   max_index = overall_max_index,
  #   normalized_index = overall_normalized_index,
  #   affinity_weighted_index = affinity_weighted_index
  # )
  #
  # if (return_pheno_values) {
  #   result$pheno_indices <- pheno_indices
  #   result$pheno_max_indices <- pheno_max_indices
  #   result$sequence_weighted_pheno_indices <- sequence_weighted_pheno_indices
  #   result$pheno_weights <- pheno_weights
  # }

  return(resdf)
}










#' Calculate the sum of matrix elements where i < j
#'
#' @description
#' This function computes the sum of elements in the upper triangular portion
#' of a matrix (excluding the diagonal). Mathematically, it calculates:
#' \deqn{\sum_{i=1}^{n} \sum_{j=i+1}^{n} a_{ij}}
#' where a_{ij} is the element at position (i,j) in an n×n matrix.
#'
#' @param matrix A numeric matrix (should be square)
#' @param diag_include Logical. Whether to include the diagonal elements (i=j) in the sum.
#'                    Default is FALSE, which gives the sum where i < j.
#'                    If TRUE, gives the sum where i ≤ j.
#'
#' @return The sum of elements where i < j (or i ≤ j if diag_include=TRUE)
#'
#' @examples
#' # Create a 3×3 matrix
#' m <- matrix(1:9, nrow = 3)
#' upper_triangle_sum(m)  # Sum elements where i < j
#' upper_triangle_sum(m, diag_include = TRUE)  # Sum elements where i ≤ j
#'
#' @export
upper_triangle_sum <- function(matrix, diag_include = FALSE) {
  # Check if input is a matrix
  if (!is.matrix(matrix)) {
    stop("Input must be a matrix")
  }

  # Check if matrix is square
  n_rows <- nrow(matrix)
  n_cols <- ncol(matrix)
  if (n_rows != n_cols) {
    warning("Input is not a square matrix. Proceeding with upper triangle of the provided matrix.")
  }

  # Method 1: Using upper.tri() function
  if (diag_include) {
    # Sum upper triangular part including diagonal (i ≤ j)
    return(sum(matrix[upper.tri(matrix, diag = TRUE)]))
  } else {
    # Sum upper triangular part excluding diagonal (i < j)
    return(sum(matrix[upper.tri(matrix, diag = FALSE)]))
  }

  # Alternative Method 2: Explicit summing using nested loops
  # This is slower but might be clearer for educational purposes
  # result <- 0
  # for (i in 1:n_rows) {
  #   start_j <- if (diag_include) i else i + 1
  #   for (j in start_j:n_cols) {
  #     if (j <= n_cols) {  # Ensure j is within bounds
  #       result <- result + matrix[i, j]
  #     }
  #   }
  # }
  # return(result)
}


#' Normalize the columns of a matrix by dividing by column sums to give per colummn
#' proportions
#'
#' @description
#' This function normalizes each column of a matrix by dividing all elements
#' in that column by the column sum. The result is a matrix where each column
#' sums to 1 (unless the original column sum was 0).
#'
#' @param mat A numeric matrix
#' @param handling_zero_cols How to handle columns with sum 0:
#'        "warning" - Leave as is and issue a warning
#'        "zero" - Leave as all zeros (default)
#'        "equal" - Replace with equal probabilities (1/nrow)
#'        "na" - Replace with NAs
#' @param na.rm Logical. Whether to ignore NA values when calculating column sums. Default is FALSE.
#'
#' @return A matrix with normalized columns
#'
#' @examples
#' # Basic usage
#' m <- matrix(1:6, nrow = 2)
#' normalize_columns(m)
#'
#' # Handling zero columns
#' m_with_zeros <- matrix(c(0, 0, 1, 2, 3, 4), nrow = 2)
#' normalize_columns(m_with_zeros, handling_zero_cols = "equal")
#'
#' @export
normalize_columns <- function(mat, handling_zero_cols = c("zero", "warning", "equal", "na"), na.rm = FALSE) {
  # Validate input
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }

  handling_zero_cols <- match.arg(handling_zero_cols)

  # Calculate column sums
  col_sums <- colSums(mat, na.rm = na.rm)

  # Check for zero column sums
  zero_cols <- which(col_sums == 0)
  if (length(zero_cols) > 0) {
    if (handling_zero_cols == "warning") {
      warning("Some columns have sum 0. These columns will not be normalized.")
    }
  }

  # Create result matrix (copy input to avoid modifying it)
  result <- mat

  # Method 1: Using a loop
  # for (j in 1:ncol(mat)) {
  #   if (col_sums[j] != 0) {
  #     result[, j] <- mat[, j] / col_sums[j]
  #   } else {
  #     # Handle zero column sums based on selected method
  #     if (handling_zero_cols == "equal") {
  #       result[, j] <- rep(1/nrow(mat), nrow(mat))
  #     } else if (handling_zero_cols == "na") {
  #       result[, j] <- rep(NA, nrow(mat))
  #     }
  #     # If "zero" or "warning", leave as is (all zeros)
  #   }
  # }

  # Method 2: Vectorized approach (more efficient)
  # For non-zero columns
  non_zero_cols <- which(col_sums != 0)
  if (length(non_zero_cols) > 0) {
    # Use sweep to divide each column by its sum
    result[, non_zero_cols] <- sweep(mat[, non_zero_cols, drop = FALSE],
                                     2, col_sums[non_zero_cols], "/")
  }

  # Handle zero-sum columns
  if (length(zero_cols) > 0) {
    if (handling_zero_cols == "equal") {
      # Replace with equal probabilities
      equal_value <- 1/nrow(mat)
      for (j in zero_cols) {
        result[, j] <- rep(equal_value, nrow(mat))
      }
    } else if (handling_zero_cols == "na") {
      # Replace with NAs
      for (j in zero_cols) {
        result[, j] <- rep(NA, nrow(mat))
      }
    }
    # If "zero" or "warning", leave as is (all zeros)
  }

  return(result)
}
