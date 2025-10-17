
#' @title Compute diversity metrics
#' @export
computeDs <- function(jdmat, aff_mat=NULL, qs=0:2){
  if (is.null(aff_mat)) {
   #stop("aff_mat is NULL")
    m <- metacommunity(jdmat)
    d1 = raw_sub_alpha(m, qs=qs)
    d2 = norm_sub_alpha(m, qs=qs)
    d3 = raw_sub_beta(m, qs=qs)
    d4 = norm_sub_beta(m, qs=qs)

    #d5 = raw_meta_alpha(m, qs=0:2)
    #d6 = norm_meta_alpha(m, qs=0:2)
    dout <- rbind(d1, d2, d3, d4)
  } else {
  s <- similarity(aff_mat, "genetic")

  m <- metacommunity(jdmat,s)

  d1 = raw_sub_alpha(m, qs=qs)
  d2 = norm_sub_alpha(m, qs=qs)
  d3 = raw_sub_beta(m, qs=qs)
  d4 = norm_sub_beta(m, qs=qs)
  d5 = norm_sub_rho(m, qs=qs)
  d6 = raw_sub_rho(m, qs=qs)


  dout <- rbind(d1, d2, d3, d4, d5, d6)
  }
  return(dout)
}



#' @title Compute diversity metrics
#' @export
computeMetaDs <- function(jdmat, aff_mat=NULL, qs=0:2){
  if (is.null(aff_mat)) {
   #stop("aff_mat is NULL")
    m <- metacommunity(jdmat)
    d1 <- raw_meta_alpha(m, qs=qs)
    d2 <- norm_meta_alpha(m, qs=qs)
    d3 <- raw_meta_beta(m, qs=qs)
    d4 <- raw_meta_beta(m, qs=qs)

    #d5 = raw_meta_alpha(m, qs=0:2)
    #d6 = norm_meta_alpha(m, qs=0:2)
    dout <- rbind(d1, d2, d3, d4)
  } else {
  s <- similarity(aff_mat, "genetic")

  m <- metacommunity(jdmat,s)

  d1 <- raw_meta_alpha(m, qs=qs)
  d2 <- norm_meta_alpha(m, qs=qs)
  d3 <- raw_meta_beta(m, qs=qs)
  d4 <- raw_meta_beta(m, qs=qs)
  d5 <- raw_meta_rho(m, qs=qs)
  d6 <- norm_meta_rho(m, qs=qs)


  dout <- rbind(d1, d2, d3, d4, d5, d6)
  }
  return(dout)
}




#' @export
computeDs_shuffle <- function(jdmat, aff_mat=NULL, qs=0:2){
  if (is.null(aff_mat)) {
    #stop("aff_mat is NULL")
    m <- metacommunity(jdmat)
    m <- repartition(m)
    d1 = raw_sub_alpha(m, qs=qs)
    d2 = norm_sub_alpha(m, qs=qs)
    d3 = raw_sub_beta(m, qs=qs)
    d4 = norm_sub_beta(m, qs=qs)

    #d5 = raw_meta_alpha(m, qs=0:2)
    #d6 = norm_meta_alpha(m, qs=0:2)
    dout <- rbind(d1, d2, d3, d4)
  } else {
    s <- similarity(aff_mat, "genetic")

    m <- metacommunity(jdmat,s)
    m <- repartition(m)
    d1 = raw_sub_alpha(m, qs=qs)
    d2 = norm_sub_alpha(m, qs=qs)
    d3 = raw_sub_beta(m, qs=qs)
    d4 = norm_sub_beta(m, qs=qs)
    d5 = norm_sub_rho(m, qs=qs)
    d6 = raw_sub_rho(m, qs=qs)


    dout <- rbind(d1, d2, d3, d4, d5, d6)
  }
  return(dout)
}




makeBoot <- function(db, jdmat, aff_mat, group, qs, nboot=100, clone="cid", min_n=10, max_n=NULL, uniform=TRUE) {
  count_col = "seq_count"
  # Tabulate clonal abundance
  clone_tab <- alakazam::countClones(db,  clone=clone, groups=group) %>%
    dplyr::mutate(clone_count=!!rlang::sym(count_col))
  group_tab <- clone_tab %>%
    group_by(!!rlang::sym(group)) %>%
    dplyr::summarize(count=sum(!!rlang::sym("clone_count"), na.rm=TRUE)) %>%
    dplyr::rename(group=!!rlang::sym(group))
  group_all <- as.character(group_tab$group)
  group_tab <- group_tab[group_tab$count >= min_n, ]
  group_keep <- as.character(group_tab$group)
  if (uniform) {
    nsam <- min(group_tab$count, max_n)
    nsam <- setNames(rep(nsam, length(group_keep)), group_keep)
  } else {
    nsam <- if (is.null(max_n)) { group_tab$count } else { pmin(group_tab$count, max_n) }
    nsam <- setNames(nsam, group_keep) }

  boot_list <- list()
  abund_list <- list()

  for (g in group_keep) {
    n <- nsam[g]

    # Tabulate group sizes
    if (!is.null(group)) {
      # Summarize groups
      group_tab <- clone_tab %>%
        group_by(!!rlang::sym(group)) %>%
        dplyr::summarize(count=sum(!!rlang::sym("clone_count"), na.rm=TRUE)) %>%
        dplyr::rename(group=!!rlang::sym(group))
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

    n <- nsam[g]
    abund_obs <- clone_tab$clone_count[clone_tab[[group]] == g]
    names(abund_obs) <- clone_tab[[clone]][clone_tab[[group]] == g]
    boot_list[[g]] <- mosaic::do(100)*bootD(mosaic::resample(names(abund_obs), size=n,replace = TRUE), jdmat, aff_mat, f=g, qs)
  }
  bootdf <- data.table::rbindlist(boot_list, idcol = "phenotype")
  return(bootdf)
}

## function that TAKES the bootstrapped sample
bootD <- function(cids, jdmat, aff_mat, f, qs) {

  s <- similarity(aff_mat[cids,cids], "genetic")

  m <- metacommunity(jdmat[cids,f],s)
  d2 = norm_sub_alpha(m, qs=qs)
  return(d2)
}



