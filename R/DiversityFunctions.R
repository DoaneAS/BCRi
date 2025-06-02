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
    message("Metacommunity matrix was normalised to sum to 1.")
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

#' distance-class
#'
#' Container for class \code{distance}.
#'
#' @field distance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as types, and elements containing the pairwise
#' distance of types
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field components list containing the components necessary to calculate
#' similarity. This list is empty when \code{precompute_dist = TRUE} when
#' calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#'
#' @name distance-class
#' @rdname distance-class
#' @exportClass distance
#'
setClass("distance",
         slots = c(distance = "matrix",
                   dat_id = "character",
                   components = "list"))


setOldClass("phylo")
#' metacommunity-class
#'
#' Container for class \code{metacommunity}.
#'
#' @field type_abundance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#' @field similarity two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as types, and elements containing the pairwise
#' similarity of types
#' @field similarity_components list containing the components necessary to
#' calculate similarity. This list is empty when \code{precompute_dist = TRUE}
#' when calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field similarity_parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#' @field ordinariness two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing the
#' ordinariness of types within subcommunities
#' @field subcommunity_weights \code{vector} of mode \code{numeric} containing
#' subcommunity weights
#' @field type_weights two-dimensional \code{matrix} of mode \code{numeric},
#' with rows as types, columns as subcommunities, and elements containing
#' weights of types within a subcommunity
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field raw_abundance [Phylogenetic] two-dimensional \code{matrix} of mode
#' \code{numeric} with rows as types, columns as subcommunities, and elements
#' containing the relative abundance of present day species
#' @field raw_structure [Phylogenetic] two-dimensional \code{matrix} of mode
#' \code{numeric} with rows as historical species, columns as present day
#' species, and elements containing historical species lengths within lineages
#' @field parameters [Phylogenetic] \code{data.frame} containing parameters
#' associated with each historic species in the phylogeny
#'
#' @name metacommunity-class
#' @rdname metacommunity-class
#' @exportClass metacommunity
#'
setClass("metacommunity",
         slots = c(type_abundance = "matrix",
                   similarity = "matrix",
                   similarity_components = "list",
                   similarity_parameters = "list",
                   ordinariness = "matrix",
                   subcommunity_weights = "numeric",
                   type_weights = "matrix",
                   dat_id = "character",
                   raw_abundance = "matrix",
                   raw_structure = "matrix",
                   parameters = "data.frame"))


#' powermean-class
#'
#' Container for class \code{powermean}.
#'
#' @field results \code{data.frame} containing rdiversity output
#' @field measure object of class \code{character} naming the diversity
#' measure being calculated
#' @field type_abundance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#' @field ordinariness two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing the
#' ordinariness of types within subcommunities
#' @field subcommunity_weights \code{vector} of mode \code{numeric} containing
#' subcommunity weights
#' @field type_weights two-dimensional \code{matrix} of mode \code{numeric},
#' with rows as types, columns as subcommunities, and elements containing
#' weights of types within a subcommunity
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field similarity_components list containing the components necessary to
#' calculate similarity. This list is empty when \code{precompute_dist = TRUE}
#' when calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field similarity_parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#'
#' @name powermean-class
#' @rdname powermean-class
#' @exportClass powermean
#'
setClass("powermean", slots = c(results = "matrix",
                                measure = "character",
                                type_abundance = "matrix",
                                ordinariness = "matrix",
                                subcommunity_weights = "vector",
                                type_weights = "matrix",
                                dat_id = "character",
                                similarity_components = "list",
                                similarity_parameters = "list"))


#' relativeentropy-class
#'
#' Container for class \code{relativeentropy}.
#'
#' @field results \code{data.frame} containing rdiversity output
#' @field measure object of class \code{character} naming the diversity
#' measure being calculated
#' @field type_abundance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#' @field ordinariness two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing the
#' ordinariness of types within subcommunities
#' @field subcommunity_weights \code{vector} of mode \code{numeric} containing
#' subcommunity weights
#' @field type_weights two-dimensional \code{matrix} of mode \code{numeric},
#' with rows as types, columns as subcommunities, and elements containing
#' weights of types within a subcommunity
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field similarity_components list containing the components necessary to
#' calculate similarity. This list is empty when \code{precompute_dist = TRUE}
#' when calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field similarity_parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#'
#' @name relativeentropy-class
#' @rdname relativeentropy-class
#' @exportClass relativeentropy
#'
setClass("relativeentropy", slots = c(results = "matrix",
                                      measure = "character",
                                      type_abundance = "matrix",
                                      ordinariness = "matrix",
                                      subcommunity_weights = "vector",
                                      type_weights = "matrix",
                                      dat_id = "character",
                                      similarity_components = "list",
                                      similarity_parameters = "list"))


#' similarity-class
#'
#' Container for class \code{similarity}.
#'
#' @field similarity two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as types, and elements containing the pairwise
#' similarity of types
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field components list containing the components necessary to calculate
#' similarity. This list is empty when \code{precompute_dist = TRUE} when
#' calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#'
#' @name similarity-class
#' @rdname similarity-class
#' @exportClass similarity
#'
setClass("similarity", slots = c(similarity = "matrix",
                                 dat_id = "character",
                                 components = "list",
                                 parameters = "list"))


#' Raw alpha (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating alpha
#' diversity.
#'
#' Values generated from \code{raw_alpha()} may be input into \code{subdiv()} and
#' \code{metadiv()} to calculate raw subcommunity and metacommunity alpha
#' diversity.
#'
#' @param meta object of class \code{metacommunity}
#'
#' @return \code{raw_alpha} returns an object of class \code{powermean}
#' @include metacommunity.R subdiv.R metadiv.R
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
#' # Calculate raw alpha component
#' a <- raw_alpha(meta)
#' subdiv(a, 1)
#' metadiv(a, 1)
#'
raw_alpha <- function(meta) {
  results <- 1 / meta@ordinariness
  powermean(results, meta, "raw alpha")
}


#' Normalised alpha (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating
#' normalised alpha diversity.
#'
#' Values generated from \code{norm_alpha()} may be input into \code{subdiv()}
#' and \code{metadiv()} to calculate normalised subcommunity and metacommunity
#' alpha diversity.
#'
#' @inheritParams raw_alpha
#'
#' @return \code{norm_alpha} returns an object of class \code{powermean}
#' @include metacommunity.R subdiv.R metadiv.R
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
#' # Calculate normalised alpha component
#' a <- norm_alpha(meta)
#' subdiv(a, 1)
#' metadiv(a, 1)
#'
norm_alpha <- function(meta) {
  ordinariness.bar <- sapply(seq_along(meta@subcommunity_weights),
                             function(x) meta@ordinariness[, x] /
                               meta@subcommunity_weights[x])
  if (!is.matrix(ordinariness.bar))
    ordinariness.bar <- as.matrix(t(ordinariness.bar))
  colnames(ordinariness.bar) <- colnames(meta@type_abundance)
  results <- 1 / ordinariness.bar

  powermean(results, meta, "normalised alpha")
}


#' Raw rho (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating raw rho
#' diversity.
#'
#' Values generated from \code{raw_rho()} may be input into \code{subdiv()} and
#' \code{metadiv()} to calculate raw subcommunity and metacommunity rho
#' diversity.
#'
#' @inheritParams raw_alpha
#'
#' @return \code{raw_rho} returns an object of class \code{powermean}
#' @include metacommunity.R subdiv.R metadiv.R
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
#' # Calculate raw rho component
#' r <- raw_rho(meta)
#' subdiv(r, 1)
#' metadiv(r, 1)
#'
raw_rho <- function(meta) {
  results <- rowSums(meta@ordinariness, na.rm = TRUE) / meta@ordinariness
  powermean(results, meta, "raw rho")
}


#' Normalised rho (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating
#' normalised rho diversity.
#'
#' Values generated from \code{norm_rho()} may be input into \code{subdiv()} and
#' \code{metadiv()} to calculate normalised subcommunity and metacommunity rho
#' diversity.
#'
#' @inheritParams raw_alpha
#'
#' @return \code{norm_rho} returns an object of class \code{powermean}
#' @include metacommunity.R subdiv.R metadiv.R
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
#' # Calculate normalised rho component
#' r <- norm_rho(meta)
#' subdiv(r, 1)
#' metadiv(r, 1)
#'
norm_rho <- function(meta) {
  ordinariness.bar <- sapply(seq_along(meta@subcommunity_weights),
                             function(x) meta@ordinariness[, x] /
                               meta@subcommunity_weights[x])
  if (!is.matrix(ordinariness.bar))
    ordinariness.bar <- as.matrix(t(ordinariness.bar))
  colnames(ordinariness.bar) <- colnames(meta@type_abundance)
  results <- rowSums(meta@ordinariness, na.rm = TRUE) / ordinariness.bar
  powermean(results, meta, "normalised rho")
}


#' Raw beta (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating raw beta
#' diversity.
#'
#' Values generated from \code{raw_beta()} may be input into \code{subdiv()} and
#' \code{metadiv()} to calculate raw subcommunity and metacommunity beta
#' diversity.
#'
#' @inheritParams raw_alpha
#'
#' @return \code{raw_beta} returns an object of class \code{relativeentropy}
#' @include metacommunity.R subdiv.R metadiv.R
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
#' # Calculate raw beta component
#' b <- raw_beta(meta)
#' subdiv(b, 1)
#' metadiv(b, 1)
#'
raw_beta <- function(meta) {
  rho <- rowSums(meta@ordinariness, na.rm = TRUE) / meta@ordinariness
  results <- 1 / rho
  relativeentropy(results, meta, "raw beta")
}


#' Normalised beta (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating
#' normalised beta diversity.
#'
#' Values generated from \code{norm_beta()} may be input into \code{subdiv()} and
#' \code{metadiv()} to calculate normalised subcommunity and metacommunity beta
#' diversity.
#'
#' @inheritParams raw_alpha
#'
#' @return \code{norm_beta} returns an object of class \code{relativeentropy}
#' @include metacommunity.R subdiv.R metadiv.R
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
#' # Calculate normalised beta component
#' b <- norm_beta(meta)
#' subdiv(b, 1)
#' metadiv(b, 1)
#'
norm_beta <- function(meta) {
  ordinariness.bar <- sapply(seq_along(meta@subcommunity_weights),
                             function(x) meta@ordinariness[, x] /
                               meta@subcommunity_weights[x])
  if (!is.matrix(ordinariness.bar))
    ordinariness.bar <- as.matrix(t(ordinariness.bar))
  colnames(ordinariness.bar) <- colnames(meta@type_abundance)
  normalised.rho <- rowSums(meta@ordinariness, na.rm = TRUE) / ordinariness.bar
  results <- 1 / normalised.rho
  relativeentropy(results, meta, "normalised beta")
}


#' Gamma (low level diversity component)
#'
#' Calculates the low-level diversity component necessary for calculating gamma
#' diversity.
#'
#' Values generated from \code{raw_gamma()} may be input into \code{subdiv()} and
#' \code{metadiv()} to calculate subcommunity and metacommunity gamma diversity.
#'
#' @inheritParams raw_alpha
#'
#' @return \code{raw_gamma} returns an object of class \code{powermean}
#' @include metacommunity.R subdiv.R metadiv.R
#' @export
#'
#' @references R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' pop <- cbind.data.frame(A = c(1,1), B = c(2,0), C = c(3,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate gamma component
#' g <- raw_gamma(meta)
#' subdiv(g, 1)
#' metadiv(g, 1)
#'
raw_gamma <- function(meta) {
  results <- rowSums(meta@ordinariness, na.rm = TRUE)
  results[results == 0] <- NaN
  results <- 1 / results
  N <- nrow(meta@type_abundance)

  results <- apply(meta@type_abundance, 2, function(x) {
    tmp <- rep(0, N)
    tmp[which(x != 0)] <- results[which(x != 0)]
    tmp
  })
  row.names(results) <- row.names(meta@type_abundance)
  results <- as.matrix(results)

  powermean(results, meta, "gamma")
}

#' Raw subcommunity alpha diversity
#'
#' Calculates similarity sensitive raw subcommunity alpha diversity (an
#' estimate of naive-community metacommunity diversity). This measure may be
#' calculated for a series of orders, represented as a vector of \code{qs}.
#'
#' @param meta object of class \code{metacommunity}
#' @param qs \code{vector} of mode \code{numeric} containing \emph{q} values
#'
#' @return \code{raw_sub_alpha} returns a standard output of class \code{rdiv}
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
#' # Calculate raw subcommunity alpha diversity
#' raw_sub_alpha(meta, 0:2)
#'
raw_sub_alpha <- function(meta, qs)
  subdiv(raw_alpha(meta), qs)


#' Normalised subcommunity alpha diversity
#'
#' Calculates similarity-sensitive normalised subcommunity alpha diversity
#' (the diversity of subcommunity \emph{j} in isolation. This measure may be
#' calculated for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
#'
#' @return \code{norm_sub_alpha} returns a standard output of class \code{rdiv}
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
#' # Calculate normalised subcommunity alpha diversity
#' norm_sub_alpha(meta, 0:2)
#'
norm_sub_alpha <- function(meta, qs)
  subdiv(norm_alpha(meta), qs)


#' Raw subcommunity beta diversity
#'
#' Calculates similarity-sensitive raw subcommunity beta diversity (the
#' distinctiveness of subcommunity \emph{j}). This measure may be calculated
#' for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
#'
#' @return \code{raw_sub_beta} returns a standard output of class \code{rdiv}
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
#' # Calculate raw subcommunity beta diversity
#' raw_sub_beta(meta, 0:2)
#'
raw_sub_beta <- function(meta, qs)
  subdiv(raw_beta(meta), qs)


#' Normalised subcommunity beta diversity
#'
#' Calculates similarity-sensitive normalised subcommunity beta diversity (an
#' estimate of the effective number of distinct subcommunities). This
#' measure may be calculated for a series of orders, represented as a vector
#' of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
#'
#' @return \code{norm_sub_beta} returns a standard output of class \code{rdiv}
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
#' # Calculate normalised subcommunity beta diversity
#' norm_sub_beta(meta, 0:2)
#'
norm_sub_beta <- function(meta, qs)
  subdiv(norm_beta(meta), qs)


#' Raw subcommunity rho diversity
#'
#' Calculates similarity-sensitive raw subcommunity rho diversity (the
#' redundancy of subcommunity \emph{j}. This measure may be calculated for
#' a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
#'
#' @return \code{raw_sub_rho} returns a standard output of class \code{rdiv}
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
#' # Calculate raw subcommunity rho diversity
#' raw_sub_rho(meta, 0:2)
#'
raw_sub_rho <- function(meta, qs)
  subdiv(raw_rho(meta), qs)


#' Normalised subcommunity rho diversity
#'
#' Calculates similarity-sensitive normalised subcommunity rho diversity (the
#' representativeness of subcommunity \emph{j}). This measure may be calculated
#' for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
#'
#' @return \code{norm_sub_rho} returns a standard output of class \code{rdiv}
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
#' # Calculate normalised subcommunity rho diversity
#' norm_sub_rho(meta, 0:2)
#'
norm_sub_rho <- function(meta, qs)
  subdiv(norm_rho(meta), qs)


#' Subcommunity gamma diversity
#'
#' Calculates similarity-sensitive subcommunity gamma diversity (the
#' contribution per individual toward metacommunity diversity). This
#' measure may be calculated for a series of orders, represented as a vector
#' of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
#'
#' @return \code{sub_gamma} returns a standard output of class \code{rdiv}
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
#' # Calculate subcommunity gamma diversity
#' sub_gamma(meta, 0:2)
#'
sub_gamma <- function(meta, qs)
  subdiv(raw_gamma(meta), qs)


#' Raw metacommunity alpha diversity
#'
#' Calculates similarity-sensitive raw metacommunity alpha diversity (the
#' naive-community metacommunity diversity). This measure may be calculated
#' for a series of orders, represented as a vector of \code{qs}.
#'
#' @inheritParams raw_sub_alpha
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
#' @inheritParams raw_sub_alpha
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
#' @inheritParams raw_sub_alpha
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
#' @inheritParams raw_sub_alpha
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
#' @inheritParams raw_sub_alpha
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
#' @inheritParams raw_sub_alpha
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
#' @inheritParams raw_sub_alpha
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


#' Calculate individual-level diversity
#'
#' Generic function for calculating individual-level diversity.
#'
#' \code{data} may be input as three different classes:
#' \itemize{
#' \item{\code{power_mean}: calculates raw and normalised subcommunity alpha, rho
#' or gamma diversity by taking the powermean of diversity components}
#' \item{\code{relativeentropy}: calculates raw or normalised subcommunity beta
#' diversity by taking the relative entropy of diversity components}
#' \item{\code{metacommunity}: calculates all subcommunity measures of diversity}
#' }
#'
#' @inheritParams subdiv
#'
#' @return \code{inddiv()} returns a standard output of class \code{rdiv}
#' @exportMethod inddiv
#'
#' @seealso \code{\link{subdiv}} for subcommunity-level diversity and
#' \code{\link{metadiv}} for metacommunity-level diversity.
#' @references Reeve, R., T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' # Define metacommunity
#' pop <- cbind.data.frame(A = c(1,1), B = c(2,0), C = c(3,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate subcommunity gamma diversity (takes the power mean)
#' g <- raw_gamma(meta)
#' inddiv(g, 0:2)
#'
#' # Calculate subcommunity beta diversity (takes the relative entropy)
#' b <- raw_beta(meta)
#' inddiv(b, 0:2)
#'
#' # Calculate all measures of individual diversity
#' inddiv(meta, 0:2)
#'
setGeneric(name = "inddiv",
           def = function(data, qs) {
             standardGeneric("inddiv")
           } )


#' @rdname inddiv
#'
setMethod(f = "inddiv", signature = "powermean",
          definition = function(data, qs) {
            output <- reshape2::melt(data@results)
            param <- data@similarity_parameters
            cbind.data.frame(measure = data@measure,
                             q = rep(qs, each = nrow(output)),
                             type_level = "type",
                             type_name = output$Var1,
                             partition_level = "subcommunity",
                             partition_name = output$Var2,
                             diversity = output$value,
                             dat_id = data@dat_id,
                             transformation = param$transform,
                             normalised = param$normalise,
                             k = param$k,
                             max_d = param$max_d,
                             stringsAsFactors = FALSE)
          } )


#' @rdname inddiv
#'
setMethod(f = "inddiv", signature = "relativeentropy",
          definition = function(data, qs) {
            output <- reshape2::melt(data@results)
            param <- data@similarity_parameters
            cbind.data.frame(measure = data@measure,
                             q = rep(qs, each = nrow(output)),
                             type_level = "type",
                             type_name = output$Var1,
                             partition_level = "subcommunity",
                             partition_name = output$Var2,
                             diversity = output$value,
                             dat_id = data@dat_id,
                             transformation = param$transform,
                             normalised = param$normalise,
                             k = param$k,
                             max_d = param$max_d,
                             stringsAsFactors = FALSE)
          } )


#' @rdname inddiv
#'
setMethod(f = "inddiv", signature = "metacommunity",
          definition = function(data, qs) {
            # Calculate terms
            div.measures <- list(raw_alpha, norm_alpha,
                                 raw_beta, norm_beta,
                                 raw_rho, norm_rho,
                                 raw_gamma)
            # Calculate subcommunity diversity
            output <- lapply(div.measures, function(x) inddiv(x(data), qs))
            do.call(rbind.data.frame, output)
          } )



#' Metacommunity
#'
#' Functions to generate a \code{metacommunity} object.
#'
#' @param partition two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing
#' the relative abundances of types in subcommunities. For phylogenetic
#' diversity, see \emph{Details}
#' @param similarity (optional) object of class \code{similarity}
#'
#' @field type_abundance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#' @field similarity two-dimensional \code{matrix} of mode \code{numeric} with
#' rows as types, columns as types, and elements containing pairwise
#' similarities between types
#' @field similarity_components list containing the components necessary to
#' calculate similarity. This list is empty when \code{precompute_dist = TRUE}
#' when calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field similarity_parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#' @field ordinariness two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing the
#' ordinariness of types within subcommunities
#' @field subcommunity_weights \code{vector} of mode \code{numeric} containing
#' subcommunity weights
#' @field type_weights two-dimensional \code{matrix} of mode \code{numeric},
#' with rows as types, columns as subcommunities, and elements containing
#' weights of types within a subcommunity
#' @field dat_ID object of class \code{character} denoting the type of diversity
#' being calculated. This can be "naive", "genetic", "taxonomic", and so on
#' @field raw_abundance [Phylogenetic] two-dimensional \code{matrix} of mode
#' \code{numeric} with rows as types, columns as subcommunities, and elements
#' containing the relative abundance of present day species
#' @field raw_structure [Phylogenetic] two-dimensional \code{matrix} of mode
#' \code{numeric} with rows as historical species, columns as present day
#' species, and elements containing historical species lengths within lineages
#' @field parameters [Phylogenetic] \code{data.frame} containing parameters
#' associated with each historic species in the phylogeny
#'
#' @return \code{metacommunity()} returns an object of class
#' \code{metacommunity} (see \emph{Fields}).
#' @name metacommunity
#' @rdname metacommunity-methods
#' @exportMethod metacommunity
#'
#' @seealso \code{\link{metacommunity-class}}
#'
#' @examples
#' # Naive-type
#' partition <- cbind(a = c(1,1,1,0,0), b = c(0,1,0,1,1))
#' row.names(partition) <- paste0("sp", 1:5)
#' partition <- partition / sum(partition)
#' meta <- metacommunity(partition)
#'
setGeneric(name = "metacommunity",
           def = function(partition, similarity) {
             standardGeneric("metacommunity")
           } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,data.frame-method
#'
setMethod(f = "metacommunity",
          signature(partition = "data.frame", similarity = "missing"),
          definition = function(partition) {
            # If similarity is data.frame, convert to matrix
            partition <- as.matrix(partition)

            metacommunity(partition)
          } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,numeric-method
#'
setMethod(f = "metacommunity",
          signature(partition = "numeric", similarity = "missing"),
          definition = function(partition) {
            # If similarity is numeric/vector, convert to matrix
            partition <- as.matrix(partition)

            metacommunity(partition)
          } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,matrix-method
#'
setMethod(f = "metacommunity",
          signature(partition = "matrix", similarity = "missing"),
          definition = function(partition) {
            # If similarity is not input, create identity matrix
            Z <- diag(1, nrow(partition))
            row.names(Z) <- row.names(partition)
            colnames(Z) <- row.names(partition)

            similarity <- new("similarity",
                              similarity = Z,
                              dat_id = "naive",
                              parameters = list(transform = NA,
                                                k = NA,
                                                normalise = NA,
                                                max_d = NA))

            metacommunity(partition, similarity)
          } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,similarity-method
#'
setMethod(f = "metacommunity",
          signature(partition = "missing", similarity = "similarity"),
          definition = function(partition, similarity) {
            # If partition is missing, assume an even distribution
            tips <- similarity$tip.label
            partition <- matrix(rep(1 / length(tips), length(tips)))
            row.names(partition) <- tips
            colnames(partition) <- "sc1"

            metacommunity(partition, similarity)
          } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,similarity-method
#'
setMethod(f = "metacommunity",
          signature(partition = "numeric", similarity = "similarity"),
          definition = function(partition, similarity) {
            partition <- as.matrix(partition)

            metacommunity(partition, similarity)
          } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,similarity-method
#'
setMethod(f = "metacommunity",
          signature(partition = "data.frame", similarity = "similarity"),
          definition = function(partition, similarity) {
            partition <- as.matrix(partition)

            metacommunity(partition, similarity)
          } )


#' @rdname metacommunity-methods
#' @aliases metacommunity,similarity-method
#'
setMethod(f = "metacommunity",
          signature(partition = "matrix", similarity = "similarity"),
          definition = function(partition, similarity) {

            # If a similarity matrix is available (within the similarity
            # object), then generate a metacommunity object in the normal way
            if (length(similarity@similarity) != 0) {

              # Check partition and simliarity matrices
              type_abundance <- check_partition(partition)
              Z <- similarity@similarity
              Z <- check_similarity(Z, partition)

              # Calculate parameters
              subcommunity_weights <- colSums(type_abundance) /
                sum(type_abundance)
              type_weights <- apply(type_abundance, 2, function(x) x / sum(x))
              Zp.j <- Z %*% type_abundance

              # Mark all of the species that have nothing similar as NaNs
              # because diversity of an empty group is undefined
              Zp.j[Zp.j == 0] <- NaN

              if (!is.matrix(type_weights)) {
                type_weights <- t(as.matrix(type_weights))
                row.names(type_weights) <- row.names(type_abundance)
              }

              return(new("metacommunity",
                         type_abundance = type_abundance,
                         similarity = Z,
                         ordinariness = Zp.j,
                         subcommunity_weights = subcommunity_weights,
                         type_weights = type_weights,
                         dat_id = similarity@dat_id,
                         similarity_components = similarity@components,
                         similarity_parameters = similarity@parameters))

              # .. else calculate branch-based phylogenetic similarity and
              # generate a metacommunity object in the normal way
            }else if (similarity@dat_id == "phybranch") {

              components <- similarity@components
              ps <- phy_struct(components$tree, partition)
              return(chainsaw(partition = partition,
                              ps = ps,
                              depth = components$tree_depth))

              # .. otherwise calculate ordinariness line by line and generate
              # a metacommunity object in the normal way
            }else {

              components <- similarity@components

              # Calculate parameters
              type_abundance <- check_partition(partition)
              subcommunity_weights <- colSums(type_abundance) /
                sum(type_abundance)
              type_weights <- apply(type_abundance, 2, function(x) x / sum(x))

              Zp.j <- lapply(seq_len(nrow(type_abundance)), function(x) {
                tmp <- get(components$ordinariness)(similarity, x)
                tmp <- matrix(tmp, nrow = 1)
                tmp %*% type_abundance
              })

              Zp.j <- do.call(rbind.data.frame, Zp.j)
              row.names(Zp.j) <- row.names(partition)

              # Mark all of the species that have nothing similar as NaNs
              # because diversity of an empty group is undefined
              Zp.j[Zp.j == 0] <- NaN
              Zp.j <- as.matrix(Zp.j)

              if (!is.matrix(type_weights)) {
                type_weights <- t(as.matrix(type_weights))
                row.names(type_weights) <- row.names(type_abundance)
              }

              return(new("metacommunity",
                         type_abundance = type_abundance,
                         ordinariness = Zp.j,
                         subcommunity_weights = subcommunity_weights,
                         type_weights = type_weights,
                         dat_id = similarity@dat_id,
                         similarity_components = similarity@components,
                         similarity_parameters = similarity@parameters))
            }
          } )


#' @rdname metacommunity-class
#' @param object object of class \code{metacommunity}
#'
setMethod(f = "show", signature = "metacommunity",
          definition = function(object) {
            cat("Object of class `metacommunity`, containing all of the data required to calculate diversity.")
          } )

#' Metacommunity-level diversity
#'
#' Generic function for calculating metacommunity-level diversity.
#'
#' \code{data} may be input as one of three different classes:
#' \itemize{
#' \item{\code{powermean}: raw or normalised metacommunity alpha, rho or gamma
#' diversity components; will calculate metacommunity-level raw or normalised
#' metacommunity alpha, rho or gamma diversity}
#' \item{\code{relativeentropy}: raw or normalised metacommunity beta
#' diversity components; will calculate metacommunity-level raw or normalised
#' metacommunity beta diversity}
#' \item{\code{metacommunity}: will calculate all metacommunity measures of
#' diversity}
#' }
#'
#' @inheritParams subdiv
#'
#' @return \code{metadiv()} returns a standard output of class \code{rdiv}
#' @exportMethod metadiv
#'
#' @seealso \code{\link{inddiv}} for type-level diversity and
#' \code{\link{subdiv}} for subcommunity-level diversity.
#' @references Reeve, R., T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' # Define metacommunity
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' pop <- pop / sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate metacommunity gamma diversity (takes the power mean)
#' g <- raw_gamma(meta)
#' metadiv(g, 0:2)
#'
#' # Calculate metacommunity beta diversity (takes the relative entropy)
#' b <- raw_beta(meta)
#' metadiv(b, 0:2)
#'
#' # Calculate all measures of metacommunity diversity
#' metadiv(meta, 0:2)
#'
setGeneric(name = "metadiv",
           def = function(data, qs) {
             standardGeneric("metadiv")
           } )


#' @rdname metadiv
#'
setMethod(f = "metadiv", signature = "powermean",
          definition = function(data, qs) {
            # Calculate
            results <- lapply(seq_along(qs), function(x) {
              subdiv <- vapply(seq_len(ncol(data@results)), function(y)
                power_mean(data@results[, y], 1 - qs[x], data@type_weights[, y]),
                numeric(1))
              power_mean(subdiv, 1 - qs[x], data@subcommunity_weights)
            })
            # Tidy up
            output <- do.call(cbind, results)
            row.names(output) <- "metacommunity"
            colnames(output) <- qs
            output <- reshape2::melt(output)
            # Output
            param <- data@similarity_parameters
            cbind.data.frame(measure = data@measure,
                             q = output$Var2,
                             type_level = "types",
                             type_name = "",
                             partition_level = "metacommunity",
                             partition_name = "",
                             diversity = output$value,
                             dat_id = data@dat_id,
                             transformation = param$transform,
                             normalised = param$normalise,
                             k = param$k,
                             max_d = param$max_d,
                             stringsAsFactors = FALSE)
          } )


#' @rdname metadiv
#'
setMethod(f = "metadiv", signature = "relativeentropy",
          definition = function(data, qs) {
            # Calculate
            results <- lapply(seq_along(qs), function(x) {
              subdiv <- vapply(seq_len(ncol(data@results)), function(y)
                power_mean(data@results[, y], qs[x] - 1, data@type_weights[, y]),
                numeric(1))
              power_mean(subdiv, 1 - qs[x], data@subcommunity_weights)
            })
            # Tidy up
            output <- do.call(cbind, results)
            row.names(output) <- "metacommunity"
            colnames(output) <- qs
            output <- reshape2::melt(output)
            # Output
            param <- data@similarity_parameters
            cbind.data.frame(measure = data@measure,
                             q = output$Var2,
                             type_level = "types",
                             type_name = "",
                             partition_level = "metacommunity",
                             partition_name = "",
                             diversity = output$value,
                             dat_id = data@dat_id,
                             transformation = param$transform,
                             normalised = param$normalise,
                             k = param$k,
                             max_d = param$max_d,
                             stringsAsFactors = FALSE)
          } )


#' @rdname metadiv
#'
setMethod(f = "metadiv", signature = "metacommunity",
          definition = function(data, qs) {
            # Calculate terms
            div.measures <- list(raw_alpha, norm_alpha,
                                 raw_beta, norm_beta,
                                 raw_rho, norm_rho,
                                 raw_gamma)
            # Calculate metacommunity diversity
            output <- lapply(div.measures, function(x) metadiv(x(data), qs))
            do.call(rbind.data.frame, output)
          } )


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


#' Calculate power mean
#'
#' Functions to coerce an object into a \code{powermean} (\code{raw_alpha()},
#' \code{norm_alpha()}, \code{raw_rho()}, \code{norm_rho()}, and/or
#' \code{raw_gamma()}).
#'
#' @param results \code{data.frame} containing rdiversity outputs associated
#' with \code{norm_alpha()}, \code{raw_alpha()}, \code{raw_rho()},
#' \code{norm_rho()}, and/or \code{raw_gamma()}
#' @param meta object of class \code{metacommunity} containing the proportional
#' abundance of types, pair-wise similarity, and other associated variables
#' @param tag object of class \code{character} naming the diversity measure
#' being calculated
#'
#' @field results \code{data.frame} containing rdiversity outputs associated
#' with \code{norm_alpha()}, \code{raw_alpha()}, \code{raw_rho()},
#' \code{norm_rho()}, and/or \code{raw_gamma()}
#' @field measure object of class \code{character} naming the diversity
#' measure being calculated
#' @field type_abundance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#' @field ordinariness two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing the
#' ordinariness of types within subcommunities
#' @field subcommunity_weights \code{vector} of mode \code{numeric} containing
#' subcommunity weights
#' @field type_weights two-dimensional \code{matrix} of mode \code{numeric},
#' with rows as types, columns as subcommunities, and elements containing
#' weights of types within a subcommunity
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field similarity_components list containing the components necessary to
#' calculate similarity. This list is empty when \code{precompute_dist = TRUE}
#' when calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field similarity_parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#'
#' @return \code{powermean(x)} returns an object of class \code{powermean}.
#' @include class-powermean.R
#'
#' @noRd
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate subcommunity raw alpha diversity (takes the powermean)
#' a <- raw_alpha(meta)
#' class(a)
#'
powermean <- function(results, meta, tag) {
  new("powermean",
      results = results,
      measure = tag,
      type_abundance = meta@type_abundance,
      ordinariness = meta@ordinariness,
      subcommunity_weights = meta@subcommunity_weights,
      type_weights = meta@type_weights,
      dat_id = meta@dat_id,
      similarity_components = meta@similarity_components,
      similarity_parameters = meta@similarity_parameters)
}



#' @rdname powermean
#' @param object object of class \code{powermean}
#' @return \code{print(x)} prints an object object of class \code{powermean}
#'
#' @noRd
#'
setMethod(f = "show", signature = "powermean",
          definition = function(object) {
            cat("Object of class powermean.")
          } )


#' Calculate relative entropy
#'
#' Functions to coerce an object into a \code{relativeentropy}
#' (\code{raw_beta()} and/or \code{norm_beta()}).
#'
#' @param results \code{data.frame} containing rdiversity outputs associated
#' with \code{raw_beta()} and/or \code{norm_beta()}
#' @param meta object of class \code{metacommunity} containing the proportional
#' abundance of types, pair-wise similarity, and other associated variables
#' @param tag object of class \code{character} naming the diversity measure
#' being calculated
#'
#' @field results \code{data.frame} containing rdiversity outputs associated
#' with \code{raw_beta()} and/or \code{norm_beta()}
#' @field measure object of class \code{character} naming the diversity
#' measure being calculated
#' @field type_abundance two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of historical species, which is
#' calculated from the proportional abundance of terminal taxa
#' @field ordinariness two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types, columns as subcommunities, and elements containing the
#' ordinariness of types within subcommunities
#' @field subcommunity_weights \code{vector} of mode \code{numeric} containing
#' subcommunity weights
#' @field type_weights two-dimensional \code{matrix} of mode \code{numeric},
#' with rows as types, columns as subcommunities, and elements containing
#' weights of types within a subcommunity
#' @field dat_id object of class \code{character} describing the class of
#' distance / similarity being used, e.g. "naive", "taxonomic", and so on
#' @field similarity_components list containing the components necessary to
#' calculate similarity. This list is empty when \code{precompute_dist = TRUE}
#' when calculating distance. When a pairwise distance matrix is too large and
#' \code{precompute_dist = FALSE}, this list contains all the information
#' required to calculate pairwise distance between types
#' @field similarity_parameters list containing parameters associated with
#' converting pairwise distances to similarities (the \code{dist2sim()}
#' arguments)
#'
#' @return object of class \code{relativeentropy}
#' @include class-relativeentropy.R
#'
#' @noRd
#'
#' @examples
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate raw subcommunity beta diversity
#' a <- raw_beta(meta)
#' class(a)
#'
relativeentropy <- function(results, meta, tag) {
  new("relativeentropy",
      results = results,
      measure = tag,
      type_abundance = meta@type_abundance,
      ordinariness = meta@ordinariness,
      subcommunity_weights = meta@subcommunity_weights,
      type_weights = meta@type_weights,
      dat_id = meta@dat_id,
      similarity_components = meta@similarity_components,
      similarity_parameters = meta@similarity_parameters)
}



#' @rdname relativeentropy
#' @param object object of class \code{relativeentropy}
#'
#' @noRd
#'
setMethod(f = "show", signature = "relativeentropy",
          definition = function(object) {
            cat("Object of class relativeentropy, containing:\n")
            cat("@results: inddiv() results\n")
            cat("@measure: measure\n")
            cat("@type_abundance: Matrix of relative abundances (",
                ncol(object@type_abundance), "subcommunities,",
                nrow(object@type_abundance), "types )\n")
            cat("@ordinariness: Matrix of type ordinariness\n")
            cat("@subcommunity_weights: Vector of subcommunity weights\n")
            cat("@type_weights: Vector of type weights\n")
          } )

#' Generate similarity object
#'
#' Container for class \code{similarity}.
#'
#' @param similarity similarity matrix
#' @param dat_id object of class \code{character} denoting the type of diversity
#' being calculated. This can be "naive", "genetic", "taxonomic", and so on
#'
#' @return \code{similarity()} returns an object of class \code{similarity}.
#' @name similarity
#' @rdname similarity-methods
#' @exportMethod similarity
#'
setGeneric(name = "similarity",
           def = function(similarity, dat_id) {
             standardGeneric("similarity")
           } )


#' @rdname similarity-methods
#' @aliases similarity
#'
setMethod(f = "similarity",
          signature(similarity = "matrix", dat_id = "character"),
          definition = function(similarity, dat_id) {
            similarity <- check_similarity(similarity)

            new("similarity",
                similarity = similarity,
                dat_id = dat_id,
                parameters = list(transform = NA,
                                  k = NA,
                                  normalise = NA,
                                  max_d = NA))
          } )


#' @rdname similarity-methods
#' @aliases similarity
#'
setMethod(f = "similarity",
          signature(similarity = "matrix", dat_id = "missing"),
          definition = function(similarity, dat_id) {
            similarity <- check_similarity(similarity)

            new("similarity",
                similarity = similarity,
                dat_id = "UserGenerated",
                parameters = list(transform = NA,
                                  k = NA,
                                  normalise = NA,
                                  max_d = NA))
          } )


#' @rdname similarity-class
#' @param object object of class \code{similarity}
#'
setMethod(f = "show", signature = "similarity",
          definition = function(object) {
            cat("Object of class `similarity`, containing either:\n (1) a similarity matrix; or\n (2) all of the data required to calculate a similarity matrix.")
          } )


#' Calculate subcommunity-level diversity
#'
#' Generic function for calculating subcommunity-level diversity.
#'
#' \code{data} may be input as one of three different classes:
#' \itemize{
#' \item{\code{powermean}: raw or normalised metacommunity alpha, rho or gamma
#' diversity components; will calculate subcommunity-level raw or normalised
#' metacommunity alpha, rho or gamma diversity}
#' \item{\code{relativeentropy}: raw or normalised metacommunity beta
#' diversity components; will calculate subcommunity-level raw or normalised
#' metacommunity beta diversity}
#' \item{\code{metacommunity}: will calculate all subcommunity measures of
#' diversity}
#' }
#'
#' @param data \code{matrix} of mode \code{numeric}; containing diversity
#' components
#' @param qs \code{vector} of mode \code{numeric} containing \emph{q} values
#'
#' @return \code{subdiv()} returns a standard output of class \code{rdiv}
#' @exportMethod subdiv
#'
#' @seealso \code{\link{inddiv}} for type-level diversity and
#' \code{\link{metadiv}} for metacommunity-level diversity.
#' @references Reeve, R., T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
#' S. Mitchell, and L. Matthews. 2016. How to partition diversity.
#' arXiv 1404.6520v3:1–9.
#'
#' @examples
#' # Define metacommunity
#' pop <- data.frame(a = c(1,3), b = c(1,1))
#' row.names(pop) <- paste0("sp", 1:2)
#' pop <- pop/sum(pop)
#' meta <- metacommunity(pop)
#'
#' # Calculate subcommunity gamma diversity (takes the power mean)
#' g <- raw_gamma(meta)
#' subdiv(g, 0:2)
#'
#' # Calculate subcommunity beta diversity (takes the relative entropy)
#' b <- raw_beta(meta)
#' subdiv(b, 0:2)
#'
#' # Calculate all measures of subcommunity diversity
#' subdiv(meta, 0:2)
#'
setGeneric(name = "subdiv",
           def = function(data, qs) {
             standardGeneric("subdiv")
           } )


#' @rdname subdiv
#'
setMethod(f = "subdiv", signature = "powermean",
          definition = function(data, qs) {
            # Calculate
            results <- lapply(seq_along(qs), function(x)
              vapply(seq_len(ncol(data@results)), function(y)
                power_mean(data@results[, y], 1 - qs[x], data@type_weights[, y]),
                numeric(1)))
            # Tidy up
            output <- do.call(cbind, results)
            row.names(output) <- colnames(data@results)
            colnames(output) <- qs
            output <- reshape2::melt(output)
            # Output
            param <- data@similarity_parameters
            cbind.data.frame(measure = data@measure,
                             q = output$Var2,
                             type_level = "types",
                             type_name = "",
                             partition_level = "subcommunity",
                             partition_name = output$Var1,
                             diversity = output$value,
                             dat_id = data@dat_id,
                             transformation = param$transform,
                             normalised = param$normalise,
                             k = param$k,
                             max_d = param$max_d,
                             stringsAsFactors = FALSE)
          } )


#' @rdname subdiv
#'
setMethod(f = "subdiv", signature = "relativeentropy",
          definition = function(data, qs) {
            # Calculate
            results <- lapply(seq_along(qs), function(x)
              vapply(seq_len(ncol(data@results)), function(y)
                power_mean(data@results[, y], qs[x] - 1, data@type_weights[, y]),
                numeric(1)))
            # Tidy up
            output <- do.call(cbind, results)
            row.names(output) <- colnames(data@results)
            colnames(output) <- qs
            output <- reshape2::melt(output)
            # Output
            param <- data@similarity_parameters
            cbind.data.frame(measure = data@measure,
                             q = output$Var2,
                             type_level = "types",
                             type_name = "",
                             partition_level = "subcommunity",
                             partition_name = output$Var1,
                             diversity = output$value,
                             dat_id = data@dat_id,
                             transformation = param$transform,
                             normalised = param$normalise,
                             k = param$k,
                             max_d = param$max_d,
                             stringsAsFactors = FALSE)
          } )


#' @rdname subdiv
#'
setMethod(f = "subdiv", signature = "metacommunity",
          definition = function(data, qs) {
            # Calculate terms
            div_measures <- list(raw_alpha, norm_alpha,
                                 raw_beta, norm_beta,
                                 raw_rho, norm_rho,
                                 raw_gamma)
            # Calculate subcommunity diversity
            output <- lapply(div_measures, function(x) subdiv(x(data), qs))
            do.call(rbind.data.frame, output)
          } )


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







#' Repartition metacommunity
#'
#' Randomly reshuffles the relative abundance of types (\emph{e.g}. species) in
#' a metacommunity (whilst maintaining the relationship between the relative
#' abundance of a particular species across subcommunities). In the case of a
#' phylogenetic metacommunity, the relative abundance of terminal taxa are
#' randomly reshuffled and the relative abundance of types (historical species)
#' are calculated from the resulting partition.
#'
#' @param meta object of class \code{metacommunity}.
#' @param new_partition two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of terminal taxa. If this argument
#' is missing, all species / tips will be shuffled
#'
#' @return \code{repartition()} returns an object of class \code{metacommunity}
#' @export
#'
repartition <- function(meta, new_partition) {

  # Non-phylogenetic metacommunity
  if (isTRUE(all.equal(0, length(meta@raw_structure)))) {

    partition <- meta@type_abundance
    if (missing(new_partition))
      new_partition <- partition[sample(seq_along(row.names(partition))), ,
                                 drop = FALSE]
    new_partition <- check_partition(new_partition)
    # Check
    if (any(dim(partition) != dim(new_partition)))
      stop("Dimensionality has changed during repartition()ing")

    row.names(new_partition) <- row.names(partition)

    s <- similarity(meta@similarity)

    new_meta <- metacommunity(new_partition, s)



    # Phylogenetic metacommunity
  }else {

    raw_abundance <- meta@raw_abundance
    if (missing(new_partition))
      new_partition <- raw_abundance[sample(seq_along(row.names(raw_abundance))),
                                     , drop = FALSE]
    new_partition <- check_partition(new_partition)
    # Check
    if (any(dim(raw_abundance) != dim(new_partition)))
      stop("Dimensionality has changed during repartition()ing")

    row.names(new_partition) <- row.names(raw_abundance)

    hs_abundance <- phy_abundance(new_partition, meta@raw_structure)

    new_similarity <- new("similarity",
                          similarity = meta@similarity * sum(hs_abundance),
                          dat_id = meta@dat_id)

    new_meta <- metacommunity(partition = hs_abundance / sum(hs_abundance),
                              similarity = new_similarity)
    new_meta@raw_abundance <- new_partition
    new_meta@raw_structure <- meta@raw_structure / sum(hs_abundance)
    new_meta@parameters <- meta@parameters
  }

  new_meta
}




#' Repartition metacommunity
#'
#' Randomly reshuffles the relative abundance of types (\emph{e.g}. species) in
#' a metacommunity (whilst maintaining the relationship between the relative
#' abundance of a particular species across subcommunities). In the case of a
#' phylogenetic metacommunity, the relative abundance of terminal taxa are
#' randomly reshuffled and the relative abundance of types (historical species)
#' are calculated from the resulting partition.
#'
#' @param meta object of class \code{metacommunity}.
#' @param new_partition two-dimensional \code{matrix} of mode \code{numeric}
#' with rows as types (species), columns as subcommunities, and each
#' element containing the relative abundance of types in each subcommunity
#' relative to the metacommunity as a whole. In the phylogenetic case, this
#' corresponds to the proportional abundance of terminal taxa. If this argument
#' is missing, all species / tips will be shuffled
#'
#' @return \code{repartition()} returns an object of class \code{metacommunity}
#' @export
#'
repartition_boot <- function(jdmat, aff_mat, nboot, qs) {

  make_p <- function(cids, partition, s) {
    new_partition <- partition[cids, ]
    row.names(new_partition) <- row.names(partition)
    new_meta <- metacommunity(new_partition, s)
    norm_sub_alpha(new_meta, qs=qs)
  }

  # Non-phylogenetic metacommunity

  s <- rdiversity::similarity(aff_mat, "genetic")

  meta <- rdiversity::metacommunity(jdmat,s)

  partition <- meta@type_abundance
  s <- similarity(meta@similarity)
  boots = mosaic::do(nboot)* make_p(mosaic::resample(row.names(partition)), partition, s)
  return(boots)
}


