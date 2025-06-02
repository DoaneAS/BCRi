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
#' @param data \code{matrix} of mode \code{numeric}; containing diversity
#' components
#' @param qs \code{vector} of mode \code{numeric} containing \emph{q} values
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
#' @inheritParams inddiv
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
#' @inheritParams inddiv
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