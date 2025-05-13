## abundance functions
##
# Helper function for computing alpha diversity from bootrstrap outputs
#
# \code{helperAlpha} divides a set of bootstrapped clones by group annotation,
# and computes the diversity of each set.
#
# @param    boot_output  data.frame from\link{AbundanceCurve} object containing bootstrapped clonal
#                        abundance curves.
# @param    q            vector of Hill Diversity indices to test for diversity calculations.
# @param    clone        name of the \code{boot_output} column containing clone identifiers.
# @param    group        name of the \code{boot_output} column containing grouping information for
#                        diversity calculation.
#
# @return   data.frame containing diversity calculations for each bootstrap iteration.
helperAlphaSS <- function(boot_output, q, clone="clone_id", group=NULL) {
  ## DEBUG
  # abundance <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
  # clone <- abundance@clone_by
  # group <- abundance@group_by

  # Compute diversity from a column of each bootstrap
  output <- boot_output %>%
    dplyr::ungroup() %>%
    dplyr::select(-one_of(c(clone, group))) %>%
    as.matrix() %>%
    apply(2, calcDiversity, q=q) %>%
    data.frame() %>%
    mutate(q=q)

  return(output)
}


# Perform boostrap abundance calculation
#
# @param    x       named vector of observed abundance values.
# @param    n       number of samples to draw from the estimate complete abundance distribution.
# @param    nboot   number of bootstrap realizations.
# @param    method  complete abundance inferrence method.
#                   One of "before", "after" or "none" for complete abundance distribution
#                   inferrence before sampling, after sampling, or uncorrected, respectively.
#
# @return   A matrix of bootstrap results.
bootstrapAbundance <- function(x, n, nboot=200, method="before") {
  ## DEBUG
  # x=abund_obs; method="before"
  # Check argumets
  method <- match.arg(method)

  if (method == "before") {
    # Calculate estimated complete abundance distribution
    p <- inferCompleteAbundance(x)
    # Bootstrap abundance
    boot_mat <- rmultinom(nboot, n, p) / n
  } else if (method == "after") {
    # Calculate estimated complete abundance distribution
    p <- x / sum(x, na.rm=TRUE)
    boot_sam <- rmultinom(nboot, n, p)
    boot_list <- apply(boot_sam, 2, inferCompleteAbundance)

    # Convert to matrix
    boot_names <- unique(unlist(sapply(boot_list, names)))
    boot_mat <- matrix(0, nrow=length(boot_names), ncol=nboot)
    rownames(boot_mat) <- boot_names
    for (i in 1:nboot) {
      boot_mat[names(boot_list[[i]]), i] <- boot_list[[i]]
    }
  } else if (method == "none") {
    # Raw sampling of input
    p <- x / sum(x, na.rm=TRUE)
    boot_sam <- rmultinom(nboot, n, p)
  } else {
    stop("Invalid method: ", method)
  }

  return(boot_mat)
}

#' Estimates the complete clonal relative abundance distribution
#'
#' \code{estimateAbundance} estimates the complete clonal relative abundance distribution
#' and confidence intervals on clone sizes using bootstrapping.
#'
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of
#'                     copy numbers within each clonal group.
#' @param    group     name of the \code{data} column containing group identifiers.
#'                     If \code{NULL} then no grouping is performed and the \code{group}
#'                     column of the output will contain the value \code{NA} for each row.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} then no
#'                     maximum is set.
#' @param    uniform   if \code{TRUE} then uniformly resample each group to the same
#'                     number of observations. If \code{FALSE} then allow each group to
#'                     be resampled to its original size or, if specified, \code{max_size}.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' @param    progress  if \code{TRUE} show a progress bar.
#'
#' @return   A \link{AbundanceCurve} object summarizing the abundances.
#'
#' @references
#' \enumerate{
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population.
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers:
#'            A framework for sampling and estimation in species diversity studies.
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by
#'            generalizing the Good-Turing sample coverage theory.
#'            Ecology. 2015 96, 11891201.
#' }
#'
#' @seealso
#' See \link{plotAbundanceCurve} for plotting of the abundance distribution.
#' See \link{alphaDiversity} for a similar application to clonal diversity.
#'
#' @examples
#' abund <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
#'
#' @export
estimateAbundance <- function(data, clone="clone_id", copy=NULL, group=NULL,
                              min_n=30, max_n=NULL, uniform=TRUE, ci=0.95, nboot=200,
                              progress=FALSE) {

  # TODO:
  # Add alakazam style cellIdColumn=NULL, locusColumn="locus", locusValues=c("IGH")
  # similar to distToNearest
  # filter based on locusValues
  # if cellIdColumn
  #    for rows that have unique cell_id, ok
  #    if rows have cell_id not unique, count only once
  # if not cellIdColumn, count heavy chains (locusValues will be IGH)
  # If mixed bulk and sc do calculation but raise warning because different type of abundances

  ## DEBUG
  # data=ExampleDb; group="sample_id"; clone="clone_id"; copy=NULL; min_n=1; max_n=NULL; ci=0.95; uniform=F; nboot=100
  # copy="duplicate_count"
  # group=NULL

  # Hack for visibility of dplyr variables
  #. <- NULL

  # Check input
  if (!is.data.frame(data)) {
    stop("Input data is not a data.frame")
  }

  # Check columns that are reported are real columns (can be NULL)
  check <- checkColumns(data, c(clone, copy, group))
  if (check != TRUE) { stop(check) }

  # Set confidence interval
  ci_z <- ci + (1 - ci) / 2
  ci_x <- qnorm(ci_z)

  # Tabulate clonal abundance
  count_col <- if (!is.null(copy)) { "copy_count" } else { "seq_count" }
  clone_tab <- countClones(data, copy=copy, clone=clone, groups=group) %>%
    dplyr::mutate(clone_count=!!rlang::sym(count_col))

  # Tabulate group sizes
  if (!is.null(group)) {
    # Summarize groups
    group_tab <- clone_tab %>%
      group_by(!!rlang::sym(group)) %>%
      dplyr::summarize(count=sum(!!rlang::sym("clone_count"), na.rm=TRUE)) %>%
      rename(group=!!rlang::sym(group))
  } else {
    group_tab <- data.frame(v="All", count=sum(clone_tab$clone_count, na.rm=T))
    names(group_tab)[1] <- "group"
  }
  group_all <- as.character(group_tab$group)
  group_tab <- group_tab[group_tab$count >= min_n, ]
  group_keep <- as.character(group_tab$group)

  # Set number of sampled sequence
  if (uniform) {
    nsam <- min(group_tab$count, max_n)
    nsam <- setNames(rep(nsam, length(group_keep)), group_keep)
  } else {
    nsam <- if (is.null(max_n)) { group_tab$count } else { pmin(group_tab$count, max_n) }
    nsam <- setNames(nsam, group_keep)
  }

  # Warn if groups removed
  if (length(group_keep) < length(group_all)) {
    warning("Not all groups passed threshold min_n=", min_n, ".",
            " Excluded: ", paste(setdiff(group_all, group_keep), collapse=", "))
  }

  # Generate abundance bootstrap
  if (progress) {
    pb <- progressBar(length(group_keep))
  }
  boot_list <- list()
  abund_list <- list()
  for (g in group_keep) {
    n <- nsam[g]

    # Extract abundance vector
    if (!is.null(group)) {
      abund_obs <- clone_tab$clone_count[clone_tab[[group]] == g]
      names(abund_obs) <- clone_tab[[clone]][clone_tab[[group]] == g]
    } else {
      # Extract abundance vector
      abund_obs <- clone_tab$clone_count
      names(abund_obs) <- clone_tab[[clone]]
    }

    # Infer complete abundance distribution
    boot_mat <- bootstrapAbundance(abund_obs, n, nboot=nboot, method="before")

    # Assign confidence intervals based on variance of bootstrap realizations
    p_mean <- apply(boot_mat, 1, mean)
    p_sd <- apply(boot_mat, 1, sd)
    p_err <- ci_x * p_sd
    p_lower <- pmax(p_mean - p_err, 0)
    p_upper <- p_mean + p_err

    # Assemble and sort abundance data.frame
    abund_df <- tibble::tibble(!!clone := rownames(boot_mat), p=p_mean, p_sd=p_sd,
                               lower=p_lower, upper=p_upper) %>%
      dplyr::arrange(desc(!!rlang::sym("p"))) %>%
      dplyr::mutate(rank=1:n())

    # Save summary
    abund_list[[g]] <- abund_df

    # Save bootstrap
    boot_list[[g]] <- as.data.frame(boot_mat) %>%
      tibble::rownames_to_column(clone)

    if (progress) { pb$tick() }
  }
  id_col <- "group"
  if (!is.null(group)) { id_col <- group }
  abundance_df <- as.data.frame(bind_rows(abund_list, .id=id_col))
  bootstrap_df <- as.data.frame(bind_rows(boot_list, .id=id_col))

  # Create a new diversity object with bootstrap
  abund_obj <- new("AbundanceCurve",
                   bootstrap=bootstrap_df,
                   abundance=abundance_df,
                   clone_by=clone,
                   group_by=id_col,
                   #groups=if_else(is.null(group), as.character(NA), group_keep),
                   groups=group_keep,
                   n=nsam,
                   nboot=nboot,
                   ci=ci)

  return(abund_obj)
}


### Coverage functions ####

#' Calculate sample coverage
#'
#' \code{calcCoverage} calculates the sample coverage estimate, a measure of sample
#' completeness, for varying orders using the method of Chao et al, 2015, falling back
#' to the Chao1 method in the first order case.
#'
#' @param    x  numeric vector of abundance counts.
#' @param    r  coverage order to calculate.
#'
#' @return   The sample coverage of the given order \code{r}.
#'
#' @references
#' \enumerate{
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population.
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by
#'            generalizing the Good-Turing sample coverage theory.
#'            Ecology. 2015 96, 11891201.
#' }
#'
#' @seealso
#' Used by \link{alphaDiversity}.
#'
#' @examples
#' # Calculate clone sizes
#' clones <- countClones(ExampleDb, groups="sample_id")
#'
#' # Calculate 1first order coverage for a single sample
#' calcCoverage(clones$seq_count[clones$sample_id == "+7d"])
#'
#' @export
calcCoverage <- function(x, r=1) {
  # Use traditional calculation for 1st order coverage
  if (r == 1) { return(calcChao1Coverage(x)) }

  # Use general form for 2nd order and higher coverage
  x <- x[x >= 1]
  n <- sum(x)
  fr <- sum(x == r)
  fs <- sum(x == r + 1)

  if (fr == 0) {
    stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", r, ".")
  }
  if (fs == 0) {
    stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", r + 1, ".")
  }

  a <- factorial(r)*fr / sum(x[x >= r]^r)
  b <- ((n - r)*fr / ((n - r)*fr + (r + 1)*fs))^r
  rC <- 1 - a*b

  return(rC)
}


# Calculate first order coverage
#
# @param    x  a numeric vector of species abundance as counts
#
# @returns  Coverage estimate.
calcChao1Coverage <- function(x) {
  x <- x[x >= 1]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)

  if (f2 > 0) {
    rC1 <- 1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
  } else {
    rC1 <- 1 - (f1 / n) * (((n - 1) * (f1 - 1)) / ((n - 1) * (f1 - 1) + 2))
  }

  return(rC1)
}


# Calculates diversity under rarefaction
#
# Calculates Hill numbers under rarefaction
#
# @param    x  vector of observed abundance counts.
# @param    m  the sequence count to rarefy to.
#
# @return   The first order coverage estimate
inferRarefiedCoverage <- function(x, m) {
  x <- x[x >= 1]
  n <- sum(x)
  if (m > n) {
    stop("m must be <= the total count of observed sequences.")
  }

  # Unrarefied case
  if (m == n) {
    return(calcCoverage(x, r=1))
  }

  # Calculate rarefied coverage
  # TODO: Read up on this and fix
  #rC1 <- iNEXT:::Chat.Ind(x, m)
  y <- x[(n - x) >= m]
  rC1 <- 1 - sum(y/n * exp(lgamma(n - y + 1) - lgamma(n - m - y + 1) - lgamma(n) + lgamma(n - m)))

  return(rC1)
}


#### Abundance functions ####

# Calculate undetected species
#
# Calculates the lower bound of undetected species counts using the Chao1 estimator.
#
# @param    x  vector of observed abundance counts.
#
# @return   The count of undetected species.
inferUnseenCount <- function(x) {
  x <- x[x >= 1]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)

  if (f2 > 0) {
    f0 <- ceiling(((n - 1) * f1^2) / (n * 2 * f2))
  } else {
    f0 <- ceiling(((n - 1) * f1 * (f1 - 1)) / (n * 2))
  }

  return(f0)
}


# Define undetected species relative abundances
#
# @param    x  vector of detected species abundance counts.
#
# @return   An adjusted detected species relative abundance distribution.
inferUnseenAbundance <- function(x) {
  x <- x[x >= 1]

  # Coverage
  rC1 <- calcCoverage(x, r=1)

  # Unseen count
  f0 <- inferUnseenCount(x)

  # Assign unseen relative abundance
  p <- rep((1 - rC1) / f0, f0)

  return(p)
}


# Adjustement to observed relative abundances
#
# @param    x  vector of observed abundance counts
#
# @return   An adjusted observed species relative abundance distribution.
adjustObservedAbundance <- function(x) {
  x <- x[x >= 1]
  n <- sum(x)

  # Coverage
  rC1 <- calcCoverage(x, r=1)

  # Calculate tuning parameter
  lambda <- (1 - rC1) / sum(x/n * exp(-x))

  # Define adjusted relative abundance
  p <- x/n * (1 -  lambda * exp(-x))

  return(p)
}


# Combined unseen inferrence and observed abundance adjustment
#
# @param    x  named vector of observed abundance counts by clone.
#
# @return   A vector containing the complete inferred abundance distribution.
#           Unseen species will be denote by a clone name starting with "U".
inferCompleteAbundance <- function(x) {
  # Infer complete abundance distribution
  p1 <- adjustObservedAbundance(x)
  p2 <- inferUnseenAbundance(x)
  names(p2) <- if (length(p2) > 0) { paste0("U", 1:length(p2)) } else { NULL }

  return(c(p1, p2))
}
