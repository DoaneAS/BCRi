
#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#'
#' @param db AIRR-compliant data.frame
#' @param cloneID ID of clone to analyze
#' @param phenotype_var phenotype variable in db for which to calculate diversity metrics
#' @param cell_id cell id column name in db. Default is "cell_id", which applies to AIRR-compliant data.frames. Use `cell_id` = NULL for bulk repertoire data.
#' @param discreteVar column name of discrete variable in db representing clonal families. Default is "clone_id". 
#' @param similarity logical indicating whether to calculate similarity matrix. Default is TRUE.
#' @param group column name of group variable in db. Default is "subject_id", which applies to AIRR-compliant data.frames. Use `group` = NULL for bulk repertoire data.
#' @param cdr3 logical indicating whether to calculate CDR3 diversity metrics. Default is FALSE.
#' @param useAffinityWeights logical indicating whether to use affinity weights. Default is TRUE.
#' @param qs vector of q values to calculate diversity metrics for. Default is 0:2.
#' @param log_base base of log to use for shannon entropy calculations. Default is exp(1).
#' @param distanceCutoff logical indicating whether to use distance cutoff. Default is FALSE.
#' @param log_base base of log to use for shannon entropy calculations
#' @param germline column name of germline sequence in db. Default is "germline_alignment".
#' @param sequence column name of sequence in db. Default is "sequence_alignment".
#' @param junction column name of junction sequence in db. Default is "junction".
#' @param v_call column name of V-call in db. Default is "v_call".
#' @param j_call column name of J-call in db. Default is "j_call".
#' @param locus column name of locus in db. Default is "locus".
#' @param only_heavy logical indicating whether to only use heavy chain sequences. Default is TRUE.
#' @param split_light logical indicating whether to split light chain sequences. Default is FALSE.
#' @param nboot number of bootstrap replicates to calculate diversity metrics. Default is 100.
#' @param subsample number of cells to subsample from the db. Default is NULL.
#' @details
#' * Clonotype refers to the unique BCR sequence.  Clonal family refers to the set of clonotypes that share the same `discreteVar` value.
#' * Diversities calculated using different approaches and returned in the `diversities` list:
#'  `naive` = BCR clonotype abundances without considering functional similarity; 
#'  `similarity_within_family` = BCR clonotype abundances and functional similarity, where similarity is considered only within clonal families as defined by `discreteVar`; 
#'  `global` = diversity BCR clonotype abundances and functional similarity considering similarity between all clonotypes; 
#'  `discrete_clonotypic` = diversity calculated using abundances that result from aggregating each clonotype by `discreteVar`.
#' @returns a [list()] containing 
#' * `diversities`: see @details 
#' *`$functional_similarity_matrix` (the functional similarity matrix).
#' *`$db_clonotype`: the db aggregated by clonotype as defined by unique BCR sequences);
#' * `P`: the joint distribution of clonotype abundances and phenotype;
#' * `P_family`: the joint distribution of clonal family abundances (as defined by `discreteVar`) and phenotype;
#' * `metaDiversities`: meta-diversity metrics calculated using the joint distributions of clonotype abundances and functional similarity with phenotype;
#' * `diversities_shuffle`: diversity metrics calculated using a shuffle approach to estimate the null distribution of diversity metrics;
#' * `diversities_boot`: diversity metrics calculated using a bootstrap approach to estimate the confidence interval of diversity metrics;
#' @export 
functional_diversity <- function(db, groupID=NULL, phenotype_var="subset", phenotype_reference=NULL,  cell_id="cell_id", similarity=TRUE,
                                  group = "subject_id", cdr3=FALSE, useAffinityWeights=TRUE,qs=0:2, germline = "germline_alignment", 
                                  sequence = "sequence_alignment",
                                  junction = "junction", v_call = "v_call", j_call = "j_call", locus = "locus", 
                                  only_heavy = TRUE, split_light = FALSE,
                                  log_base=exp(1), distanceCutoff=FALSE, discreteVar="clone_id", nboot=100, subsample=NULL) {

  # get clone
  print(groupID)
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  normalize = "len" ## c("len", "none"),
  fields = NULL

  targeting_model = shazam::HH_S5F
  len_limit = NULL
  first = FALSE
  #cdr3 = FALSE
  mod3 = FALSE
  max_n = 0
  threshold = 1
  base_sim = 0.95
  iter_max = 1000
  nstart = 1000
  nproc = 4
  verbose = FALSE
  log = NULL
  summarize_clones = FALSE
  indVar = "ind"

  db_clone <- as.data.frame(db[db[[group]] == groupID, ])

  if(!is.null(subsample)){
    if(nrow(db_clone) > subsample){
      cat("Subsampling ", nrow(db_clone), " cells to ", subsample, " cells... \n")
      ix <- sample(1:nrow(db_clone), size=subsample)
      db_clone <- db_clone[ix,]
    }
  }

  results_prep = prepare_clone(db = db_clone,
                               junction = junction, v_call = v_call, j_call = j_call,
                               first = first, cdr3 = cdr3, fields = fields,
                               cell_id = cell_id, locus = locus, only_heavy = only_heavy,
                               mod3 = mod3, max_n = max_n)
  dbfull = data.table::copy(db)
  db <- results_prep$db
  n_rmv_mod3 <- results_prep$n_rmv_mod3
  n_rmv_cdr3 <- results_prep$n_rmv_cdr3
  n_rmv_N <- results_prep$n_rmv_N
  junction_l <- results_prep$junction_l
  cdr3_col <-  results_prep$cdr3_col

  db_l <- db[db[[locus]] %in% c("IGK", "IGL", "TRA", "TRG"), , drop=F]
  db <- db[db[[locus]] %in% c("IGH", "TRB", "TRD"), , drop=F]
  db_gp <- db

  mutabs <- shazam::HH_S5F@mutability
  # Generated by using Rcpp::compileAttributes() -> do not edit by hand
  # Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

  #Rcpp::sourceCpp("~/Projects/code/scratch/scoper/src/RcppMutation.cpp") # scoper/src/RcppMutation.cpp")
  #Rcpp::sourceCpp("./src/RcppMutation.cpp")

  n <- nrow(db_gp)


  germs <- db_gp[[germline]]
  seqs <- db_gp[[sequence]]
  ## accomdate different length sequences for hamming distance
  juncs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
  #juncs <- as.character( set_eq_seqDistance(juncs))
  junc_length <- unique(stringi::stri_length(juncs))
  if (length(junc_length) > 1) {
    juncs <- as.character( set_eq_seqDistance(juncs))
    junc_length <- unique(stringi::stri_length(juncs))
  }


  #juncs <- as.character( set_eq_seqDistance(juncs))
  # find unique seqs
  #.I = NULL
  seqs <- paste(seqs, juncs, germs, sep = "|")
  df <- data.table::as.data.table(seqs)[, list(list(.I)), by=seqs] %>%
    tidyr::separate(col = seqs, into = c("seqs_unq", "juncs_unq", "germs_unq"), sep = "\\|")
  n_unq <- nrow(df)
  ind_unq <- df$V1

  ## make result df here
  dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))
  dgc[[phenotype_var]] = unique(db_gp[[phenotype_var]])
  #dgc[[phenotype_var]] = unique(db_gp[[phenotype_var]])
  dgc$rgsw = 0
  dgc$rgsw_norm = 0
  dgc$clone_id = groupID
  dgc$clone_size = n
  dgc$richness <- n_unq
  dgc$type = "limited"
  dgc$wmax = 0
  dgc$beta1 = 0
  dgc$clonotypic_entropy_sim = 0
  dgc$clonotypic_entropy_sim_mnorm = 0
  dgc$clonotypic_entropy_sim_lnorm = 0
  dgc$clonotypic_entropy_sim_mlnorm = 0
  dgc$clonotypic_entropy <- 0
  dgc$clonotypic_entropy_mnorm <- 0
  dgc$clonotypic_entropy_lnorm <- 0
  dgc$clonotypic_entropy_mlnorm <- 0
  dgc$clonotypic_diversity <- 0
  dgc$clonotypic_diversity_sim <- 0



  if (n_unq <= 2) {

    return(dgc)


    return_list <- list("affinity_mat" = NULL,
                        "disim_mtx" = NULL,
                        "db_clone" = NULL,
                        "jd" = NULL,
                        "weighted_GS" =  NULL,
                        "db_pheno"= dgc)

    return(return_list)

  } else{
    dgc$type = "full"

    # find corresponding unique germs and junctions
    seqs_unq <- df$seqs_unq
    germs_unq <- df$germs_unq
    juncs_unq <- df$juncs_unq
    cat("Calculating distance matrix...\n")
    # calculate unique junctions distance matrix
    dist_mtx <- alakazam::pairwiseDist(seq = juncs_unq,
                                       dist_mat = getDNAMatrix(gap = 0))
    # count mutations from unique sequence imgt
    results <- pairwiseMutions(germ_imgt = germs_unq,
                               seq_imgt = seqs_unq,
                               junc_length = junc_length,
                               len_limit = len_limit,
                               cdr3 = cdr3,
                               mutabs = mutabs)
    tot_mtx <- results$pairWiseTotalMut
    sh_mtx <- results$pairWiseSharedMut
    mutab_mtx <- results$pairWiseMutability
    # calculate likelihhod matrix
    lkl_mtx <- likelihoods(tot_mtx = tot_mtx,
                           sh_mtx = sh_mtx,
                           mutab_mtx = mutab_mtx)
    # calculate weighted matrix
    disim_mtx <- dist_mtx * (1.0 - lkl_mtx)


    mtx = disim_mtx

    nearest_dist <- apply(mtx, 2,  function(x) {
      gt0 <- which(x > 0)
      if (length(gt0) != 0) { min(x[gt0]) } else { NA }
    })


    krnl_mtx <- krnlMtxGenerator(mtx = mtx)


    if (distanceCutoff) {
      aff_mtx <- makeAffinity(mtx_o = mtx,
                              mtx_k = krnl_mtx,
                              thd = max(nearest_dist, na.rm = T))
      #aff_mtx[aff_mtx > max(nearest_dist, na.rm = T)] <- 0
    } else {
      aff_mtx <- krnl_mtx
    }

    threshold = NULL
    base_sim = 0.95
    iter_max = 1000
    nstart = 1000
    #{
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])


    aff_mtx[is.na(aff_mtx)] <- 0

    ## return clone db with unique seq identifier
    db_gp$ind <- 0

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }


    ## option to zero out across discrete categories
    if (!is.null(discreteVar)) {
      cids = list()
      for (i in 1:n_unq) {
        cids[[i]] = db_gp[[discreteVar]][db_gp[["ind"]]==i][1]

      }


      cids <- unlist(cids)
      cmat = outer(cids , cids, "==")
      mode(cmat) = "numeric"

      d_mtx = aff_mtx * cmat
    } else {
      null_mtx =  aff_mtx
      null_mtx[] <- 0
      diag(null_mtx) <- 1
      d_mtx <- null_mtx
      cmat <- null_mtx
    }


    if (!similarity){
      aff_mtx = d_mtx
    }




    dbj = bcrCounts(db_gp, inds = "ind")
    # prob_mat =  aff_mtx
    # prob_mat[] <- 0
    #
    # for (i in 1:nrow(dbj)) {
    #   for (j in 1:nrow(dbj)) {
    #     prob_mat[i, j] = dbj$p[i] * dbj$p[j]
    #   }
    # }

    dbjp = bcrCounts_pheno(db_gp, inds = "ind", pheno=phenotype_var)

    dbj$w = 0
    dbj$wn = 0
    dbj$wnv = 0
    db_gp$w = 0
    db_gp$wn = 0
    dbjp$w = 0
    dbjp$wn = 0


    kdmatrix = 1 - aff_mtx

    for (i in dbj$ind) {
      #dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
      #dbj$wn[dbj$ind==i] = sum(aff_mtx[i, ])/sum(aff_mtx) # s
      dbj$w[dbj$ind==i] = sum(aff_mtx[i, ])  # sum of affinities for each cell
      dbj$wn[dbj$ind==i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell * richness

      dbjp$w[dbjp$ind==i] = sum(aff_mtx[i, ])
      dbjp$wn[dbjp$ind==i] = sum(aff_mtx[i, ]) * nrow(dbjp)
      db_gp$w[db_gp$ind == i] = sum(aff_mtx[i, ])  # sum of affinities for each cell
      db_gp$wn[db_gp$ind == i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
      #db_gp$wn[db_gp$ind == i] = sum(aff_mtx[i, ]) / sum(aff_mtx)
    }

    wijmat = matrix(0, nrow = nrow(dbj), ncol = nrow(dbj))

    dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
    dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)
    #dbp$wn <- dbp$w / sum(dbp$w)
    #
    dbp = db_gp %>% dplyr::group_by(.data[[indVar]], .data[["w"]],  .data[["wn"]], .data[[phenotype_var]]) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::spread(key = {{phenotype_var}}, value = n, fill=0)



    jdmat <-  db_gp %>% dplyr::group_by(.data[[indVar]], .data[[phenotype_var]]) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::spread(key = {{phenotype_var}}, value = n, fill=0) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data[[indVar]]) %>% as.matrix()

    jdmat_c <- db_gp %>% dplyr::group_by(.data[[discreteVar]],  .data[[phenotype_var]]) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::spread(key = {{phenotype_var}}, value = n, fill=0) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data[[discreteVar]]) %>% as.matrix()

    dbp_p = dbp

    ###proportion private to phenotype
    for (f in unique(db_gp[[phenotype_var]])) {
      dbp_p[[f]] = dbp_p[[f]] / sum(dbp_p[[f]])
    }
    ###

    for (f in unique(db_gp[[phenotype_var]])) {
      rgslist=list()
      wlist = list()
      #dlist = list()
      #dwlist = list()
      ix = dbp_p[[indVar]]
      if (sum(dbp_p[[f]]) == 0) {
        next
      }

      for (i in ix) {
        #jx = setdiff(ix, i)
        wn = dbp_p$wn[i]
        #d = dbp_p$d[i]
        #dw = dbp_p$dw[i]
        #pj = sum(dbp_p[[f]][jx])
        rgslist[i] =  dbp_p[[f]][i] * (1-dbp_p[[f]][i]) * wn
        wlist[i] = wn
        #dlist[i] =  dbp_p[[f]][i] * (1-dbp_p[[f]][i]) * d
        # dwlist[i] =  dbp_p[[f]][i] * (1-dbp_p[[f]][i]) * dw
      }

      if (sum(unlist(rgslist)) == 0) {
        next
      }

      dgc$rgsw[dgc$pheno==f] =  sum(unlist(rgslist))
      dgc$wmax[dgc$pheno==f] = max(unlist(wlist))
    }

    dgc$beta1 = dgc$wmax * (1 - (1/n_unq))
    dgc$rgsw_norm = dgc$rgsw / dgc$beta1
 
    cat("calculating diversities...")
    probs_matrix <- normalize_columns(jdmat)
    totp = rowSums(jdmat)
    totP = totp / sum(totp)
    cell_n = colSums(jdmat)
    probs_matrixT <- cbind(totP, probs_matrix)
    colnames(probs_matrixT)[1] = "Mixed"
    rownames(aff_mtx) <- paste0("cid_", dbp$ind)
    colnames(aff_mtx) <- rownames(aff_mtx)
    rownames(jdmat) <- paste0("cid_", dbp$ind)
    db_gp$cid <- paste0("cid_", db_gp$ind)


    ddf = computeDs(jdmat, aff_mtx, qs=qs)
    nmat <- aff_mtx
    nmat[] <- 0
    diag(nmat) <- 1
    ddf_n = computeDs(jdmat, aff_mat=NULL, qs=qs)
    ddf_i = computeDs(jdmat, aff_mat=d_mtx, qs=qs)

    ddf_clone <- computeDs(jdmat_c,aff_mat=NULL, qs=qs)
    ddf_n$sim = "naive"
    ddf_i$sim = "similarity_within_family"
    ddf$sim = "global"
    ddf_clone$sim <- "discrete_clonotypic"
    ddf_shuffle <- computeDs_shuffle(jdmat,aff_mtx, qs=qs)
    ddf_shuffle$sim <- "shuffle"

    ddf_boot = repartition_boot(jdmat = jdmat, aff_mat = aff_mtx, nboot = nboot, qs=qs)
    ddf_boot$groupID <- groupID

    ddfm <- computeMetaDs(jdmat, aff_mat=aff_mtx, qs=qs)
    ddfm$groupID <- groupID

    db_gp$cid <- paste0("cid_", db_gp$ind)
    db_gp$clone <- db_gp$cid
    bres <- makeBoot(db_gp, jdmat=jdmat, aff_mat=aff_mtx, group=phenotype_var, qs=qs, nboot=nboot, clone="cid", min_n=2, max_n=NULL)

    bres$groupID <- groupID
    ddf <- rbind(ddf, ddf_n, ddf_i, ddf_clone)
    ddf$groupID <- groupID

    cat("completed\n")

    return_list <- list("functional_similarity_matrix" = aff_mtx,
                        "disim_mtx" = disim_mtx,
                        "db_clonotype" = db_gp,
                        "P" = jdmat,
                        "P_family" = jdmat_c,
                        "diversities" = ddf,
                        ##"diversities_corrected" = ddf_a,
                        "metaDiversities" = ddfm,
                        "diversity_shuffle" = ddf_boot,
                        "diversities_bootstrap" = bres
    )

    return(return_list)
  }
}

#' Run functional_diversity using metadata stored in a Seurat object
#'
#' Extracts the cell-level metadata from a Seurat object, optionally remaps
#' column names to the defaults expected by `functional_diversity`, and then
#' forwards the augmented metadata to `functional_diversity`.
#'
#' @param seurat_obj A Seurat or SeuratObject containing the required metadata.
#' @rdname functional_diversity
#' @param groupID Identifier (or identifiers) of the group in the column supplied via `group`.
#' @param phenotype_var Name of the metadata column describing phenotypes; passed to `functional_diversity`.
#' @param phenotype_reference Optional phenotype label to use as reference; passed to `functional_diversity`.
#' @param group Name of the metadata column describing the grouping variable; passed to `functional_diversity`.
#' @param column_map Named list mapping canonical column names used by `functional_diversity`
#'   (`clone_id`, `germline_alignment`, `sequence_alignment`, `junction`, `locus`) to the
#'   corresponding column names present in the Seurat metadata. Defaults assume metadata already
#'   uses the canonical names.
#' @param metadata Optional data.frame to use instead of extracting metadata from `seurat_obj`.
#' @param ... Additional arguments forwarded to `functional_diversity`.
#'
#' @return The result of calling `functional_diversity` on the extracted metadata.
#' @export
functional_diversity_from_seurat <- function(seurat_obj,
                                             groupID,
                                             phenotype_var = "subset",
                                             phenotype_reference = NULL,
                                             group = "subject_id",
                                             column_map = list(
                                               clone_id = "clone_id",
                                               germline_alignment = "germline_alignment",
                                               sequence_alignment = "sequence_alignment",
                                               junction = "junction",
                                               locus = "locus"
                                             ),
                                             metadata = NULL,
                                             ...) {

  dots <- list(...)

  if (!inherits(seurat_obj, c("Seurat", "SeuratObject"))) {
    stop("`seurat_obj` must inherit from `Seurat` or `SeuratObject`.")
  }

  if (is.null(metadata)) {
    metadata <- tryCatch(seurat_obj[[]], error = function(e) NULL)
    if (is.null(metadata)) {
      metadata <- tryCatch(seurat_obj@meta.data, error = function(e) NULL)
    }
  }

  if (is.null(metadata) || !is.data.frame(metadata)) {
    stop("Unable to extract metadata from `seurat_obj`; supply it via the `metadata` argument.")
  }

  db <- as.data.frame(metadata, stringsAsFactors = FALSE)

  required_map_names <- c("clone_id", "germline_alignment", "sequence_alignment", "junction", "locus")
  missing_map_names <- setdiff(required_map_names, names(column_map))
  if (length(missing_map_names) > 0) {
    stop("`column_map` must provide entries for: ", paste(missing_map_names, collapse = ", "), ".")
  }

  required_cols <- unique(c(group, phenotype_var, unlist(column_map)))
  missing_cols <- setdiff(required_cols, colnames(db))
  if (length(missing_cols) > 0) {
    stop("The following required metadata columns are missing: ", paste(missing_cols, collapse = ", "), ".")
  }

  for (canonical_name in names(column_map)) {
    source_col <- column_map[[canonical_name]]
    if (!identical(canonical_name, source_col) || !canonical_name %in% colnames(db)) {
      db[[canonical_name]] <- db[[source_col]]
    }
  }

  if (!"cell_id" %in% colnames(db)) {
    db$cell_id <- rownames(db)
  }

  if (is.null(groupID) && !"groupID" %in% names(dots)) {
    stop("`groupID` must be supplied to identify the group to analyse.")
  }

  default_args <- list(
    db = db,
    groupID = groupID,
    phenotype_var = phenotype_var,
    phenotype_reference = phenotype_reference,
    group = group
  )

  if (!"cell_id" %in% names(dots)) {
    default_args$cell_id <- "cell_id"
  }

  if (!"discreteVar" %in% names(dots) && "clone_id" %in% names(db)) {
    default_args$discreteVar <- "clone_id"
  }

  args <- utils::modifyList(default_args, dots)

  do.call(functional_diversity, args)
}

