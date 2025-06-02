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

# Show methods for classes
#' @rdname metacommunity-class
#' @param object object of class \code{metacommunity}
#'
setMethod(f = "show", signature = "metacommunity",
          definition = function(object) {
            cat("Object of class `metacommunity`, containing all of the data required to calculate diversity.")
          } )

#' @rdname powermean-class
#' @param object object of class \code{powermean}
#' @return \code{print(x)} prints an object object of class \code{powermean}
#'
#' @noRd
#'
setMethod(f = "show", signature = "powermean",
          definition = function(object) {
            cat("Object of class powermean.")
          } )

#' @rdname relativeentropy-class
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

#' @rdname similarity-class
#' @param object object of class \code{similarity}
#'
setMethod(f = "show", signature = "similarity",
          definition = function(object) {
            cat("Object of class `similarity`, containing either:\n (1) a similarity matrix; or\n (2) all of the data required to calculate a similarity matrix.")
          } )