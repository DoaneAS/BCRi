#' Helper Function for pairwise mutations, adapted from scoper
#'
#'
pairwiseMutions <- function(germ_imgt,
                            seq_imgt,
                            junc_length,
                            len_limit = NULL,
                            cdr3 = FALSE,
                            mutabs = NULL,
                            norm_fact = TRUE) {

  ##### get number of seqs
  n <- unique(c(length(seq_imgt), length(germ_imgt)))
  ##### check number of sequences
  if (length(n) > 1) stop("germ_imgt and seq_imgt number should be the same")
  if (n == 1) stop("there should be at least two seqs")
  # check consensus length
  if (!is.null(len_limit)) {
    lenConsensus <- len_limit@seqLength
    seq_imgt  <- substr(seq_imgt, start = 1, stop = lenConsensus)
    germ_imgt <- substr(germ_imgt, start = 1, stop = lenConsensus)
    eff_germ <- ifelse(length(unique(germ_imgt)) == 1,
                       unique(germ_imgt),
                       shazam::consensusSequence(sequences = unique(germ_imgt),
                                                 muFreqColumn = NULL,
                                                 lenLimit = lenConsensus,
                                                 method = "catchAll",
                                                 minFreq = NULL,
                                                 includeAmbiguous = FALSE,
                                                 breakTiesStochastic = FALSE,
                                                 breakTiesByColumns = NULL,
                                                 db = NULL)$cons)
  } else {
    ##### constants
    lv <- ifelse(cdr3, shazam::IMGT_V@seqLength, shazam::IMGT_V@seqLength - 3)
    trim_l <- junc_length
    ##### trim out junction/cdr3 segments from seq_imgt
    seq_imgt <- sapply(1:length(seq_imgt), function(i){
      x <- strsplit(seq_imgt[i], split="")[[1]]
      x[(lv+1):(lv+trim_l)] <- ""   # x[(lv+1):(lv+trim_l[i])] <- ""
      return(paste(x, collapse=""))
    })
    ##### Pads ragged ends
    l <- unique(stringi::stri_length(seq_imgt))
    if (length(l) > 1) {
      seq_imgt <- padSeqEnds(seq = seq_imgt, len = NULL, start = FALSE, pad_char = "N")
    }
    ##### trim out junction/cdr3 segments from germ_imgt
    germ_imgt <- sapply(1:length(germ_imgt), function(i){
      x <- strsplit(germ_imgt[i], split="")[[1]]
      x[(lv+1):(lv+trim_l)] <- ""  # x[(lv+1):(lv+trim_l[i])] <- ""
      return(paste(x, collapse=""))
    })
    ##### Pads ragged ends
    l <- unique(stringi::stri_length(germ_imgt))
    if (length(l) > 1) {
      germ_imgt <- padSeqEnds(seq = germ_imgt, len = NULL, start = FALSE, pad_char = "N")
    }
    ##### find consensus germline (allel level grouping)
    # see arg "method" from shazam::collapseClones function
    eff_germ <- ifelse(length(unique(germ_imgt)) == 1,
                       unique(germ_imgt),
                       consensusSequence(sequences = unique(germ_imgt),
                                         muFreqColumn = NULL,
                                         lenLimit = NULL,
                                         method = "catchAll",
                                         minFreq = NULL,
                                         includeAmbiguous = FALSE,
                                         breakTiesStochastic = FALSE,
                                         breakTiesByColumns = NULL,
                                         db = NULL)$cons)
    ##### check germ and seqs lengths
    seq_imgt_lent <- unique(stringi::stri_length(seq_imgt))
    germ_imgt_lent <- unique(stringi::stri_length(germ_imgt))
    eff_germ_lent <- stringi::stri_length(eff_germ)
    lenConsensus <- min(seq_imgt_lent, germ_imgt_lent, eff_germ_lent)
    ##### trim extra characters
    if ( seq_imgt_lent > lenConsensus)  { seq_imgt <- substr( seq_imgt, start = 1, stop = lenConsensus) }
    if (germ_imgt_lent > lenConsensus) { germ_imgt <- substr(germ_imgt, start = 1, stop = lenConsensus) }
    if ( eff_germ_lent > lenConsensus)  { eff_germ <- substr( eff_germ, start = 1, stop = lenConsensus) }
  }
  ##### count informative positions
  if (norm_fact) {
    informative_pos <- sapply(1:n, function(x){ sum(stringi::stri_count(seq_imgt[x], fixed = c("A","C","G","T"))) })
  } else {
    informative_pos <- rep(1, n)
  }
  ##### convert eff_germ and seq_imgt to matrices
  seqsMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
  effMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
  for (i in 1:n) {
    seqsMtx[i, ] <- strsplit(seq_imgt[i], split = "")[[1]][1:lenConsensus]
    effMtx[i, ] <- strsplit(eff_germ, split = "")[[1]][1:lenConsensus]
  }
  ##### make a distance matrix
  dnaMtx <- getDNAMatrix(gap = 0)
  mutMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
  for (i in 1:n) {
    mutMtx[i, ] <- sapply(1:lenConsensus, function(j) {
      return(dnaMtx[effMtx[i,j], seqsMtx[i,j]])
    })
  }
  ##### make a mutation matrix
  mutMtx <- matrix(paste0(effMtx, mutMtx, seqsMtx), nrow=n, ncol=lenConsensus)
  ##### clean non-mutated elements
  mutMtx[grepl(pattern="0", mutMtx)] <- NA
  ##### check mutabilities
  ##### make a motif matrix
  motifMtx <- matrix(0, nrow=n, ncol=lenConsensus)
  if (!is.null(mutabs)) {
    for (i in 1:n) {
      for (j in 3:(lenConsensus-2)) {
        motifMtx[i, j] <- mutabs[substr(germ_imgt[i], start = j-2, stop = j+2)]
      }
    }
    motifMtx[is.na(motifMtx)] <- 0
  }
  ##### calculate mutation matrix
  results <- pairwiseMutMatrix(informative_pos = informative_pos,
                               mutMtx = mutMtx,
                               motifMtx = motifMtx)
  sh_mtx <- results$sh_mtx
  tot_mtx <- results$tot_mtx
  mutab_mtx <- results$mutab_mtx
  ##### make symmetric matrix
  sh_mtx[lower.tri(sh_mtx)] <- t(sh_mtx)[lower.tri(sh_mtx)]
  tot_mtx[lower.tri(tot_mtx)] <- t(tot_mtx)[lower.tri(tot_mtx)]
  mutab_mtx[lower.tri(mutab_mtx)] <- t(mutab_mtx)[lower.tri(mutab_mtx)]
  # return results
  return_list <- list("pairWiseSharedMut" = sh_mtx,
                      "pairWiseTotalMut" = tot_mtx,
                      "pairWiseMutability" = mutab_mtx)
  return(return_list)
}
# *****************************************************************************

pairwiseMutMatrix <- function(informative_pos, mutMtx, motifMtx) {
  pairwiseMutMatrixRcpp(informative_pos, mutMtx, motifMtx)
}
# ***




#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#'
#' @param db db
#' @param cloneID specific clone
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name for clone
#' @importFrom igraph diversity degree graph_from_adjacency_matrix simplify
#' @returns list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using Shannon diversity
#' @export
intraclonal_shannon <- function(db,
                                cloneID,
                                phenotype_var="subset",
                                cell_id=NULL,
                                clone = "clone_id",
                                callClones = FALSE,
                                normalize="len",
                                germline = "germline_alignment",
                                sequence = "sequence_alignment",
                                junction = "junction",
                                v_call = "v_call",
                                j_call = "j_call",
                                fields = NULL,
                                locus = "locus",
                                only_heavy = TRUE,
                                split_light = FALSE,
                                targeting_model = NULL,
                                len_limit = NULL,
                                first = FALSE,
                                cdr3 = FALSE,
                                mod3 = FALSE,
                                max_n = 0,
                                threshold = NA,
                                base_sim = 0.95,
                                iter_max = 1000,
                                nstart = 1000,
                                nproc = 4,
                                verbose = FALSE,
                                log = NULL,
                                summarize_clones = FALSE) {
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  #normalize = "len" ## c("len", "none"),
  #germline = "germline_alignment_d_mask"
  #sequence = "sequence_alignment"
  #junction = "junction"
  #v_call = "v_call_genotyped"
  #j_call = "j_call"
  #fields = NULL
  #locus = "locus"
  #only_heavy = TRUE
  #split_light = FALSE
  if (!is.null(targeting_model)) {
    if (targeting_model == "none") {
      targeting_model = NULL
    } else {
      mutabs <- targeting_model@mutability
      #targeting_model = targeting_model@mutability
    }
  } else {
    targeting_model = shazam::HH_S5F
    mutabs <- targeting_model@mutability
  }

  # if callClones {
  #   scoper::spectralClones(db = db,
  #                          model = model, method = method, linkage = linkage,
  #                          normalize = normalize, germline = germline,
  #                          sequence = sequence, junction = junction,
  #                          v_call = v_call, j_call = j_call,
  #                          fields = fields, locus = locus,
  #                          only_heavy = only_heavy, split_light = split_light,
  #                          targeting_model = targeting_model,
  #                          len_limit = len_limit, first = first,
  #                          cdr3 = cdr3, mod3 = mod3, max_n = max_n,
  #                          threshold = threshold, base_sim = base_sim,
  #                          iter_max = iter_max, nstart = nstart,
  #                          nproc = nproc, verbose = verbose,
  #                          log=log, summarize_clones=summarize_clones)
  #
  # }
  #targeting_model = HH_S5F
  #len_limit = NULL
  #first = FALSE
  #cdr3 = FALSE
  #mod3 = FALSE
  #max_n = 0
  #threshold = 1
  #base_sim = 0.95
  #iter_max = 1000
  #nstart = 1000
  #nproc = 4
  #verbose = FALSE
  #log = NULL
  #summarize_clones = FALSE

  # get clone


  print(cloneID)
  db_clone <- as.data.frame(db[db$clone_id == cloneID, ])
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
  #groupBy <- c("vjl_group")

  ### summary of the groups
  # vjl_gps <- db %>%
  #   dplyr::group_by(!!!rlang::syms(groupBy)) %>%
  #   dplyr::summarise(group_v_call = stringi::stri_join(unique(!!rlang::sym(v_call)), collapse=","),
  #                    group_j_call = stringi::stri_join(unique(!!rlang::sym(j_call)), collapse=","),
  #                    group_junction_length = unique(!!rlang::sym(junction_l)),
  #                    group_size = n())
  # vjl_gps$group_v_call <- sapply(1:nrow(vjl_gps),
  #                                function(i){ stringi::stri_join(unique(stringi::stri_split_fixed(vjl_gps$group_v_call[i], ",")[[1]]), collapse=",") })
  # vjl_gps$group_j_call <- sapply(1:nrow(vjl_gps),
  #                                function(i){ stringi::stri_join(unique(stringi::stri_split_fixed(vjl_gps$group_j_call[i], ",")[[1]]), collapse=",") })
  # n_groups <- nrow(vjl_gps)

  ## a single clone will have only 1 vjl group
  # gp <- vjl_gps$vjl_group[1]
  #
  # len_limit = NULL
  # vjl_gp <- vjl_gps$vjl_group[gp]
  # gp_vcall <- vjl_gps$group_v_call[gp]
  # gp_jcall <- vjl_gps$group_j_call[gp]
  # gp_lent <- vjl_gps$group_junction_length[gp]
  # gp_size <- vjl_gps$group_size[gp]
  # db_gp <- dplyr::filter(db, !!rlang::sym("vjl_group") == vjl_gp)
  #
  mutabs <- shazam::HH_S5F@mutability
  mutabs <- targeting_model@mutability
  # Generated by using Rcpp::compileAttributes() -> do not edit by hand
  # Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

  #Rcpp::sourceCpp("./src/RcppMutation.cpp")

  n <- nrow(db_gp)

  ### cloning
  # if (method == "vj") {
  #   ### check targeting model
  #   if (!is.null(mutabs)) {
  #     mutabs <- mutabs@mutability
  #   } else {
  #     mutabs <- NULL
  #   }
  # get required info based on the method
  germs <- db_gp[[germline]]
  seqs <- db_gp[[sequence]]
  juncs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
  junc_length <- unique(stringi::stri_length(juncs))
  # find unique seqs
  seqs <- paste(seqs, juncs, germs, sep = "|")
  df <- data.table::as.data.table(seqs)[, list(list(.I)), by=seqs] %>%
    tidyr::separate(col = seqs, into = c("seqs_unq", "juncs_unq", "germs_unq"), sep = "\\|")
  n_unq <- nrow(df)
  ind_unq <- df$V1
  if (n_unq == 1) {
    db_gp$ind <- 0
    db_gp$Hi = 0
    db_gp$vertex_diversity = 0
    db_gp$gll = 0

    ind = "ind"

    dgc = db_gp %>% dplyr::group_by({{ ind }}, {{ phenotype_var }}) %>%
      dplyr::summarise(n = n())

    return_list = list("affinity_mat" = NULL,
                       "db_clone" = db_gp,
                       "g" = NULL,
                       "db_pheno"= dgc)

    return(return_list)
  }
  # find corresponding unique germs and junctions
  seqs_unq <- df$seqs_unq
  germs_unq <- df$germs_unq
  juncs_unq <- df$juncs_unq
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
  threshold = NULL
  base_sim = 0.95
  iter_max = 1000
  nstart = 1000
  mtx = disim_mtx
  #{
  ### constants
  n <- nrow(mtx)
  bs <- (1 - base_sim)*junc_length
  off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])

  krnl_mtx <- krnlMtxGenerator(mtx = mtx)
  aff_mtx <- krnl_mtx

  ## return clone db with unique seq identifier
  db_gp$ind <- 0

  for (i in 1:n_unq) {
    #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
    db_gp$ind[ind_unq[[i]]] <- i
  }

  g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
  g =  igraph::simplify(g, remove.loops = TRUE, remove.multiple = FALSE)
  igraph::V(g)$k = igraph::degree(g)

  #v = igraph::eigen_centrality(g)$vector
  #db_gp$ev = v[ db_gp$ind]

  cent.div <- igraph::diversity(g)


  igraph::V(g)$Di <-  igraph::diversity(g)
  igraph::V(g)$Hi =  igraph::V(g)$Di * log10(igraph::V(g)$k)
  db_gp$Di = igraph::V(g)$Di[ db_gp$ind]
  db_gp$Hi = igraph::V(g)$Hi[ db_gp$ind]

  # for (i in 1:length(g)) {
  #   #spij = sum(E(g)[.from(i)]$weight / log(V(g)$k[i]))
  #   #spij = sum(g[i,] / log(V(g)$k[i]))
  #   V(g)$sumPij[i] = sum( g[i,] / log(V(g)$k[i]))
  #   #set_vertex_attr(g, index = V(g)[i], name = "sumPij", value = spij)
  #   #V(g)$sumPij[i] =  sum(g[i,] / log(V(g)$k[i]))
  # }


  #gll = list()
  for (i in 1:length(g)) {
    ## proportional weight
    pij = g[i,] / sum(g[i,])
    pij = pij[pij > 0]
    Hiw = (-1 * sum(pij * log10(pij)))
    Diw = Hiw / log10(igraph::V(g)$k[i])
    igraph::V(g)$Diw[i] = Diw
    igraph::V(g)$Hiw[i] = Hiw
    igraph::V(g)$sumw[i] = sum(pij * log10(pij))

    #pijg = g[i,] / sum(g[])
    # V(g)$pijg[i] = sum(pijg)

    ### different norms for edge weights
    # pijk = g[i,] / V(g)$k[i]
    # pijk = pijk[pijk > 0]
    # Hik = -1 * sum(pijk * log(pijk))
    # Dik = Hik / log10(V(g)$k[i])
    # V(g)$sumPijk[i] = sum(pijk)
    #
    # V(g)$Hik[i] = Hik
    # V(g)$Dik[i] = Dik
    #
    #
    # piju = g[i,]
    # piju = piju[piju > 0]
    #
    # Hiu = -1 * sum(piju * log(piju))
    # Diu = Hiu / log10(V(g)$k[i])
    # V(g)$sumPiju[i] = sum(piju)
    #
    # V(g)$Hiu[i] = Hiu
    # V(g)$Diu[i] = Diu

    #gll[i] = sum(pij) * log(sum(pij))
  }


  #db_gp$sumPij = V(g)$sumPij[ db_gp$ind]
  db_gp$Diw = igraph::V(g)$Diw[ db_gp$ind]
  db_gp$Hiw = igraph::V(g)$Hiw[ db_gp$ind]
  #db_gp$Dik = V(g)$Dik[ db_gp$ind]
  #db_gp$Hik = V(g)$Hik[ db_gp$ind]
  #db_gp$sumPijk = V(g)$sumPijk[ db_gp$ind]
  #db_gp$sumw = V(g)$sumw[ db_gp$ind]
  #db_gp$pijg = V(g)$pijg[ db_gp$ind]

  #db_gp$Diu = V(g)$Diu[ db_gp$ind]
  #db_gp$Hiu = V(g)$Hiu[ db_gp$ind]
  #db_gp$sumPiju = V(g)$sumPiju[ db_gp$ind]
  db_gp$n_clone = n_unq
  ind = "ind"
  #  grps = c(ind, phenotype_var)
  #dgc = db_gp %>% dplyr::group_by({{phenotype_var}}) %>%
  #  dplyr::summarise(n = n())
  #dgc$intra_clonotypic = 0

  dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))
  #dgc[[phenotype_var]] = unique(db_gp[[phenotype_var]])
  dgc$intra_clonotypic = 0
  dgc$inter_clonotypic_norm = 0
  dgc$sumPij = 0
  ## now get joint phenotype - entropy
  for (f in unique(db_gp[[phenotype_var]])) {
    fnodes = db_gp[db_gp[[phenotype_var]] == f,]
    fnodesindex = unique(fnodes$ind)

    pij = g[fnodesindex,] / sum(g[fnodesindex,])
    pij = pij[pij > 0]
    k = length(pij)
    Hiw = (-1 * sum(pij * log10(pij)))
    Diw = Hiw / log10(k)

    dgc$intra_clonotypic_entropy[dgc$pheno == f] = Hiw
    dgc$intra_clonotypic_diversity[dgc$pheno == f] = Diw

    sumpij =  sum(g[fnodesindex,])
    dgc$sumPij[dgc$pheno==f] = sumpij

    dgc$intra_cl[dgc$pheno == f] = sum(g[fnodesindex,]) * log10(sum(g[fnodesindex,]))


    ## joint distribution ####
    ## not needed- very few exact matches
    # dbj = db_gp %>% dplyr::group_by(ind, subset) %>%
    #   dplyr::summarise(n = n()) %>% tidyr::spread(key = subset, value = n, fill=0)
    #

    # each pij for nodes that are phenotype == f
    #pijg = db_gp[db_gp[[phenotype_var]] == f,]$pijg
    #-1 * sum(pijg * log10(pijg)


    #for (i in fnodesindex) {
    # sum(g[i,]) * log10(sum(g[fnodesindex,]))
    #}
    #db_gp$gll[i] = sum(g[i,]) * log10(sum(g[i,]))
  }


  return_list <- list("affinity_mat" = aff_mtx,
                      "db_clone" = db_gp,
                      #"g" = g,
                      "db_pheno"= dgc)

  return(return_list)
}



#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#'
#' @param db db
#' @param cloneID ID of clone to analyze
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name of clone variable in db
#' @export
#' @returns data.frame of diversity metrics and phenotype ##list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using Simposon diversity
intraclonal_simpson <- function(db, cloneID=NULL, phenotype_var="subset", cell_id=NULL, clone = "clone_id", cdr3=FALSE) {
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  normalize = "len" ## c("len", "none"),
  germline = "germline_alignment"
  sequence = "sequence_alignment"
  junction = "junction"
  v_call = "v_call"
  j_call = "j_call"
  fields = NULL
  locus = "locus"
  only_heavy = TRUE
  split_light = FALSE
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

  # get clone
  print(cloneID)
  db_clone <- as.data.frame(db[db$clone_id == cloneID, ])
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
  juncs <- as.character( set_eq_seqDistance(juncs))
  junc_length <- unique(stringi::stri_length(juncs))
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
  dgc$intra_clonotypic_entropy_sim = 0
  dgc$intra_clonotypic_entropy_sim_cells = 0
  dgc$intra_clonotypic_entropy_sim_seq = 0
  dgc$simpson_index = 0
  dgc$clone_id = cloneID
  dgc$clone_size = n

  if (n_unq <= 2) {
    db_gp$ind <- 0
    db_gp$n_clone = n_unq

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }


    #db_gp$vertex_diversity = 0
    #db_gp$gll = 0

    ind = "ind"

    #dgc = db_gp %>% dplyr::group_by({{ ind }}, {{ phenotype_var }}) %>%
    #  dplyr::summarise(n = n())

    #dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))


    #dgc$intra_clonotypic_entropy_sim = 0

    return_list = list("affinity_mat" = NULL,
                       "db_clone" = db_gp,
                       # "g" = NULL,
                       "db_pheno"= dgc)

    return(dgc)
  } else{

    # find corresponding unique germs and junctions
    seqs_unq <- df$seqs_unq
    germs_unq <- df$germs_unq
    juncs_unq <- df$juncs_unq
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

    if (all(disim_mtx == dist_mtx) | any(rowSums(disim_mtx) == 0)) {
      # get required info based on the method
      seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
      ## set seq lengths with pad
      seqs <- as.character( set_eq_seqDistance(seqs))
      junc_length <- unique(stringi::stri_length(seqs))
      # find unique seqs
      df <- data.table::as.data.table(seqs)[, list(list(.I)), by = seqs]
      n_unq <- nrow(df)
      ind_unq <- df$V1
      seqs_unq <- df$seqs
      if (n_unq == 1) {
        return(dgc)
        #return(list("idCluster" = rep(1, n),
        #            "n_cluster" = 1,
        #            "eigen_vals" = rep(0, n)))
      }
      # calculate unique seuences distance matrix
      disim_mtx <- alakazam::pairwiseDist(seq = seqs_unq,
                                          dist_mat = getDNAMatrix(gap = 0))
    }
    mtx = disim_mtx
    threshold = NULL
    base_sim = 0.95
    iter_max = 1000
    nstart = 1000
    mtx = disim_mtx
    #{
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])

    krnl_mtx <- krnlMtxGenerator(mtx = mtx)
    aff_mtx <- krnl_mtx

    aff_mtx[is.na(aff_mtx)] <- 0

    ## return clone db with unique seq identifier
    db_gp$ind <- 0

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }

    # aff_mtx_c <- aff_mtx[db_gp$ind,db_gp$ind]
    #
    #
    # #diag(aff_mtx) <- 0
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE, diag = FALSE)
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx[db_gp$ind,db_gp$ind], mode = "undirected", weighted = TRUE)
    # g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
    #
    #
    # #g =  igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
    # g =  igraph::simplify(g)
    #
    # igraph::V(g)$k = igraph::degree(g)
    #
    # #v = igraph::eigen_centrality(g)$vector
    # #db_gp$ev = v[ db_gp$ind]
    #
    # igraph::V(g)$Di <-  igraph::diversity(g)
    # igraph::V(g)$Hi =  igraph::V(g)$Di * log10(igraph::V(g)$k)
    #

    # for (i in 1:length(g)) {
    #   #spij = sum(E(g)[.from(i)]$weight / log(V(g)$k[i]))
    #   #spij = sum(g[i,] / log(V(g)$k[i]))
    #   V(g)$sumPij[i] = sum( g[i,] / log(V(g)$k[i]))
    #   #set_vertex_attr(g, index = V(g)[i], name = "sumPij", value = spij)
    #   #V(g)$sumPij[i] =  sum(g[i,] / log(V(g)$k[i]))
    # }

    dbj = bcrCounts(db_gp, inds = "ind")
    prob_mat =  aff_mtx
    prob_mat[] <- 0

    for (i in 1:nrow(dbj)) {
      for (j in 1:nrow(dbj)) {
        prob_mat[i, j] = dbj$p[i] * dbj$p[j]
      }
    }


    diag(aff_mtx) = 1

    ## Pi * Pj * Sij, where p is the probability of a cell with BCR sequence x having that sequence
    ## Useful to correct for cases where more than 1 cell has the same BCR sequence
    smat = prob_mat * aff_mtx

    dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
    dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)




    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp$ind[ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #     #dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = 0
    #   }
    #   dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = sum(smat[ix,])
    # }


    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp[[indVar]][ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   for (i in ix) {
    #     # for (j in 1:nrow(dbj)) {
    #     pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ])
    #   }
    #   dgc$intra_clonotypic_entropy_sim2[dgc$pheno==f] =  sum(unlist(pilist))
    # }


    for (f in unique(db_gp[[phenotype_var]])) {
      ix = dbp[[indVar]][ dbp[[f]] >=1]
      if (length(ix) == 0) {
        next
      }
      pilist = list()
      pcslist = list()
      pslist = list()
      for (i in ix) {
        #for (j in 1:nrow(dbj)) {
        #dbp[[f]][ix] * dbj$p[j]

        pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ]) # pi*pj*wij over i and all js
        pslist[i] = sum(aff_mtx[i,])
        pclist = sum(prob_mat[i,])
      }
      dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  1 - sum(unlist(pilist))
      dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
      dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
      dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
      dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
      dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
      dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
      dgc$simpson_index[dgc$pheno==f] = 1 - sum(dpp[[f]]^2)

    }


    return_list <- list("affinity_mat" = aff_mtx,
                        "db_clone" = db_gp,
                        #"g" = g,
                        "db_pheno"= dgc)

    return(dgc)
  }
}





#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#'
#' @param db db
#' @param cloneID ID of clone to analyze
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name of clone variable in db
#' @export
#' @returns data.frame of diversity metrics and phenotype ##list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using Simposon diversity
weighted_simpson <- function(db, cloneID=NULL, phenotype_var="subset", cell_id=NULL, clone = "clone_id", cdr3=FALSE, useAffinityWeights=TRUE) {
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  normalize = "len" ## c("len", "none"),
  germline = "germline_alignment"
  sequence = "sequence_alignment"
  junction = "junction"
  v_call = "v_call"
  j_call = "j_call"
  fields = NULL
  locus = "locus"
  only_heavy = TRUE
  split_light = FALSE
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

  # get clone
  print(cloneID)
  db_clone <- as.data.frame(db[db$clone_id == cloneID, ])
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
  dgc$clone_id = cloneID
  dgc$clone_size = n
  dgc$richness <- n_unq
  dgc$type = "limited"
  dgc$wmax = 0
  dgc$beta1 = 0

  if (n_unq <= 2) {

    return(dgc)

    # db_gp$ind <- 0
    # db_gp$n_clone = n_unq
    #
    # for (i in 1:n_unq) {
    #   #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
    #   db_gp$ind[ind_unq[[i]]] <- i
    # }
    #
    # seqs_unq <- df$seqs_unq
    # germs_unq <- df$germs_unq
    # juncs_unq <- df$juncs_unq
    #
    #
    #
    # ind = "ind"
    #
    # dbj = bcrCounts(db_gp, inds = "ind")
    # dbjp = bcrCounts_pheno(db_gp, inds = "ind", pheno=phenotype_var)
    # aff_mtx = matrix(1, nrow = n_unq, ncol = n_unq)
    # diag(aff_mtx) = 0
    #
    # dbj$w = 0
    # dbj$wn = 0
    # dbj$wnv = 0
    # db_gp$w = 0
    # db_gp$wn = 0
    # db_gp$wnv = 0
    # dbjp$w = 0
    # dbjp$wn = 0
    # dbjp$wnv = 0
    #
    # for (i in dbj$ind) {
    #   dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
    #   dbj$wn[dbj$ind==i] = sum(aff_mtx[i, ])/sum(aff_mtx) # s
    #   dbj$wnv[dbj$ind==i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
    #
    #   dbjp$w[dbjp$ind==i] = sum(aff_mtx[i, ])
    #   dbjp$wn[dbjp$ind==i] = sum(aff_mtx[i, ]) / sum(aff_mtx) # sum of affinities for each cell)
    #   db_gp$w[db_gp$ind == i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
    #   db_gp$wn[db_gp$ind == i] = sum(aff_mtx[i, ]) / sum(aff_mtx)
    # }
    # dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
    # dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)
    # dbp$wn <- dbp$w / sum(dbp$w)
    # dbp_p = dbp
    #
    # for (f in unique(db_gp[[phenotype_var]])) {
    #   dbp_p[[f]] = dbp_p[[f]] / sum(dbp[[f]])
    # }
    #
    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbjp[[indVar]][ dbjp[[phenotype_var]] == f]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   rgslist=list()
    #   wnvlist=list()
    #   pilist = list()
    #   pilist2 = list()
    #   ilist = list()
    #   pnlist = list()
    #   pcslist = list()
    #   pslist = list()
    #   wlist = list()
    #   wnlist = list()
    #   pplist = list()
    #   for (i in ix) {
    #     #for (j in 1:nrow(dbj)) {
    #     #dbp[[f]][ix] * dbj$p[j]
    #     w = dbj$w[dbj$ind==i]
    #     wn = dbj$wn[dbj$ind==i]
    #     wnv = dbj$wnv[dbj$ind==i]
    #
    #     ixx = ((dbjp[[indVar]]==i) & (dbjp[[phenotype_var]]==f))
    #     pii = dbjp$p[ixx] * sum(dbjp$p[!ixx]) ## pi * pj
    #
    #     ## pi and pj within phenotype.
    #     piii = dbp_p[[f]][dbp_p$ind==i] * sum(dbp_p[[f]][dbp_p$ind!=i]) ## pi * pj
    #
    #     #ixp = ((dbjp[[indVar]]!=i) & (dbjp[[phenotype_var]]==f)) ##
    #
    #     # sum of affinities for each cell
    #     #pii = dpp[dpp[[indVar]]==i,][[f]] * sum(dpp[dpp[[indVar]]!=i,][[f]])
    #     if (w > 0) {
    #       pilist[i] = pii * w
    #       ilist[i] = pii * (1/wn)
    #       pnlist[i] = pii * wn
    #       wnlist[[i]] = wn
    #       wlist[[i]] = w
    #       pplist[i] = piii * wn
    #       wnvlist[[i]] = wnv
    #       rgslist[i] = pii * wnv
    #     } else if (w == 0) {
    #       pilist[i] = pii
    #       ilist[i] = pii
    #       pnlist[i] = pii
    #       pplist[i] = piii
    #       wnvlist[[i]] = 0
    #       rgslist[i] = 0
    #     }
    #     #pilist[i] = pii * (1 / dbjp$w[ixx])  # pi*pj*wij over i and all js
    #     #pilist[i] = sum( prob_mat[i, ]  * sum(aff_mtx[i, ])) # pi*pj*wij over i and all js
    #     # pslist[i] = dbjp$w[ixx]
    #     pcslist[i] = pii
    #   }
    #   dgc$rgsw[dgc$pheno==f] =  sum(unlist(rgslist))
    #   dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
    #   dgc$intra_clonotypic_entropy_simn[dgc$pheno==f] =  sum(unlist(pnlist))
    #   dgc$intra_clonotypic_entropy_simi[dgc$pheno==f] =  sum(unlist(ilist))
    #   dgc$intra_clonotypic_entropy_seqwn[dgc$pheno==f] =  sum(unlist(wnlist))
    #   dgc$intra_clonotypic_entropy_seqw[dgc$pheno==f] =  sum(unlist(wlist))
    #   dgc$intra_clonotypic_entropy_simn_private[dgc$pheno==f] =  sum(unlist(pplist))
    #   #dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
    #   #dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
    #   dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
    #   dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
    #   dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
    #   dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
    #   dgc$simpson_index[dgc$pheno==f] = sum(unlist(pcslist))
    #   dgc$type[dgc$pheno==f] = "limited"
    #   dgc$wmax[dgc$pheno==f] = max(unlist(wnvlist))
    #
    # }
    # ix0 = dgc$wmax==0
    # dgc$rgsw[ix0] = 0
    # dgc$rgsw_norm[ix0] = 0
    # dgc$beta1[ix0] = 0
    # dgc$beta1[!ix0] = dgc$wmax[!ix0] * (1 - (1/n_unq))
    # dgc$rgsw_norm[!ix0] = dgc$rgsw[!ix0] / dgc$beta1[!ix0]
    #
    # dgc$beta1[!ix0] = dgc$wmax[!ix0] * (1 - (1/n_unq))
    # dgc$rgsw_norm[!ix0] = dgc$rgsw[!ixo] / dgc$beta1[!ix0]
    #
    # #dgc$intra_clonotypic_entropy_sim = 0
    #
    # return_list = list("affinity_mat" = NULL,
    #                    "db_clone" = db_gp,
    #                    # "g" = NULL,
    #                    "db_pheno"= dgc)
    #
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
    aff_mtx <- krnl_mtx
   ### use only for global
    # aff_mtx <- makeAffinity(mtx_o = mtx,
     #                       mtx_k = krnl_mtx,
      #                      thd = max(nearest_dist, na.rm = T))



    threshold = NULL
    base_sim = 0.95
    iter_max = 1000
    nstart = 1000
    #{
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])


    #aff_mtx <- krnl_mtx

    aff_mtx[is.na(aff_mtx)] <- 0

    ## return clone db with unique seq identifier
    db_gp$ind <- 0

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }

    diag(aff_mtx) = 1

    dbj = bcrCounts(db_gp, inds = "ind")
    prob_mat =  aff_mtx
    prob_mat[] <- 0

    for (i in 1:nrow(dbj)) {
      for (j in 1:nrow(dbj)) {
        prob_mat[i, j] = dbj$p[i] * dbj$p[j]
      }
    }

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


    ## normalize 0-1

    #dbj$w <- rangeAtoB(dbj$w, 0, 1)
    #dbj$wn <- dbj$w / sum(dbj$w)
    #dbjp$wn <- dbjp$w / sum(dbjp$w)

    #db_gp$w <- rangeAtoB(db_gp$w, 0, 1)
    #dbjp$w <- rangeAtoB(dbjp$w, 0, 1)
    #diag(aff_mtx) = 1


    ## Pi * Pj * Sij, where p is the probability of a cell with BCR sequence x having that sequence
    ## Useful to correct for cases where more than 1 cell has the same BCR sequence
    smat = prob_mat * aff_mtx

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


    regw <- weighted_rich_gini_simpson_affinity_pairs(jdmat,
                                                      sequence_weights = kdmatrix,
                                                    #sequence_weights = disim_mtx,
                                                    affinity_weights = dbp$w,
                                                    #useRichness= TRUE,
                                                    #normalize_phenos = TRUE,
                                                    pheno_weights = NULL,
                                                    zero_handling = "ignore")



    rgw <- weighted_rich_gini_simpson_affinity(jdmat,
                                                      affinity_weights = dbp$w,
                                                      #useRichness= TRUE,
                                                      #normalize_phenos = TRUE,
                                                      #pheno_weights = NULL,
                                                      zero_handling = "ignore")


    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp$ind[ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #     #dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = 0
    #   }
    #   dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = sum(smat[ix,])
    # }


    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp[[indVar]][ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   for (i in ix) {
    #     # for (j in 1:nrow(dbj)) {
    #     pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ])
    #   }
    #   dgc$intra_clonotypic_entropy_sim2[dgc$pheno==f] =  sum(unlist(pilist))
    # }

    # browser()

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

      #dgc$rgsw_d[dgc$pheno==f] =  sum(unlist(dlist))
     # dgc$rgsw_dw[dgc$pheno==f] =  sum(unlist(dwlist))
    }

    dgc$beta1 = dgc$wmax * (1 - (1/n_unq))
    #dgc$dmax = max(dbp_p$d)
    #dgc$dwmax = max(dbp_p$dw)
    dgc$rgsw_norm = dgc$rgsw / dgc$beta1

    #dgc$rgsw_d_norm = dgc$rgsw_d / dgc$beta1
    #dgc$rgsw_dw_norm = dgc$rgsw_dw / (dgc$dwmax * (1 - (1/n_unq)))
    #dgc$rgsw_d_norm = dgc$rgsw_dw / (dgc$dmax * (1 - (1/n_unq)))



    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   pilist2 = list()
    #   ilist = list()
    #   pnlist = list()
    #   pcslist = list()
    #   pslist = list()
    #   wlist = list()
    #   wnlist = list()
    #   pplist = list()
    # }
    #
    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbjp[[indVar]][ dbjp[[phenotype_var]] == f]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   rgslist=list()
    #   wnvlist=list()
    #   pilist = list()
    #   pilist2 = list()
    #   ilist = list()
    #   pnlist = list()
    #   pcslist = list()
    #   pslist = list()
    #   wlist = list()
    #   wnlist = list()
    #   pplist = list()
    #   for (i in ix) {
    #     #for (j in 1:nrow(dbj)) {
    #     #dbp[[f]][ix] * dbj$p[j]
    #     w = dbj$w[dbj$ind==i]
    #     #wn = dbj$wn[dbj$ind==i]
    #     #wnv = dbj$wnv[dbj$ind==i]
    #
    #     ixx = ((dbjp[[indVar]]==i) & (dbjp[[phenotype_var]]==f))
    #     pii = dbjp$p[ixx] * sum(dbjp$p[!ixx]) ## pi * pj
    #
    #     ## pi and pj within phenotype.
    #     piii = dbp_p[[f]][dbp_p$ind==i] * sum(dbp_p[[f]][dbp_p$ind!=i]) ## pi * pj
    #
    #     #ixp = ((dbjp[[indVar]]!=i) & (dbjp[[phenotype_var]]==f)) ##
    #
    #     # sum of affinities for each cell
    #     #pii = dpp[dpp[[indVar]]==i,][[f]] * sum(dpp[dpp[[indVar]]!=i,][[f]])
    #     if (w > 0) {
    #       pilist[i] = pii * w
    #       ilist[i] = pii * (1/wn)
    #       pnlist[i] = pii * wn
    #       wnlist[[i]] = wn
    #       wlist[[i]] = w
    #       pplist[i] = piii * wn
    #       wnvlist[[i]] = wnv
    #       rgslist[i] = pii * wnv
    #     } else if (w == 0) {
    #       pilist[i] = pii
    #       ilist[i] = pii
    #       pnlist[i] = pii
    #       pplist[i] = piii
    #       wnvlist[[i]] = 0
    #     }
    #
    #     #pilist[i] = pii * (1 / dbjp$w[ixx])  # pi*pj*wij over i and all js
    #     #pilist[i] = sum( prob_mat[i, ]  * sum(aff_mtx[i, ])) # pi*pj*wij over i and all js
    #     # pslist[i] = dbjp$w[ixx]
    #     pcslist[i] = pii
    #   }
      # dgc$rgsw[dgc$pheno==f] =  sum(unlist(rgslist))
      # dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
      # dgc$intra_clonotypic_entropy_simn[dgc$pheno==f] =  sum(unlist(pnlist))
      # dgc$intra_clonotypic_entropy_simi[dgc$pheno==f] =  sum(unlist(ilist))
      # dgc$intra_clonotypic_entropy_seqwn[dgc$pheno==f] =  sum(unlist(wnlist))
      # dgc$intra_clonotypic_entropy_seqw[dgc$pheno==f] =  sum(unlist(wlist))
      # dgc$intra_clonotypic_entropy_simn_private[dgc$pheno==f] =  sum(unlist(pplist))
      # #dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
      # #dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
      # dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
      # dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
      # dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
      # dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
      # dgc$simpson_index[dgc$pheno==f] = sum(unlist(pcslist))
      # dgc$type[dgc$pheno==f] = "full"
      # dgc$wmax[dgc$pheno==f] = max(unlist(wnvlist))

    #}



    setDT(dgc)
    dgc <- dgc[rgsw >0,]

    return_list <- list("affinity_mat" = aff_mtx,
                        "disim_mtx" = disim_mtx,
                         "db_clone" = db_gp,
                        "jd" = dbp,
                        "weighted_RGS_pairs" =  regw,
                        "weighted_RGS" = rgw,
                         "db_pheno"= dgc)

    #setDT(dgc)
    #dgc <- dgc[rgsw >0,]
    return(return_list)
  }
}




#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#'
#' @param db db
#' @param cloneID ID of clone to analyze
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name of clone variable in db
#' @export
#' @returns data.frame of diversity metrics and phenotype ##list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using Simposon diversity
weighted_simpson_global <- function(db, cloneID=NULL, phenotype_var="subset", cell_id=NULL, clone = "clone_id", cdr3=FALSE) {
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  normalize = "len" ## c("len", "none"),
  germline = "germline_alignment"
  sequence = "sequence_alignment"
  junction = "junction"
  v_call = "v_call"
  j_call = "j_call"
  fields = NULL
  locus = "locus"
  only_heavy = TRUE
  split_light = FALSE
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

  # get clone
  print(cloneID)

  results_prep = prepare_clone(db = db,
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
  dgc$intra_clonotypic_entropy_sim = 0
  dgc$intra_clonotypic_entropy_simi = 0
  dgc$intra_clonotypic_entropy_simn = NA
  dgc$intra_clonotypic_entropy_simn_private = 0
  dgc$intra_clonotypic_entropy_sim_cells = 0
  dgc$intra_clonotypic_entropy_sim_seq = 0
  dgc$intra_clonotypic_entropy_seqw = 0
  dgc$intra_clonotypic_entropy_seqwn = NA
  dgc$simpson_index = 0
  dgc$clone_id = cloneID
  dgc$clone_size = n
  dgc$type = "limited"

  #if (n_unq <= 2) {
  #  db_gp$ind <- 0
  #  db_gp$n_clone = n_unq
  #
  #  for (i in 1:n_unq) {
  #    #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
  #    db_gp$ind[ind_unq[[i]]] <- i
  #  }
  #
  #  seqs_unq <- df$seqs_unq
  #  germs_unq <- df$germs_unq
  #  juncs_unq <- df$juncs_unq
  #
  #
  #
  #  ind = "ind"
  #  dbjp = bcrCounts_pheno(db_gp, inds = "ind", pheno=phenotype_var)
  #
  #
  #  dbj = bcrCounts(db_gp, inds = "ind")
  #
  #  dbj$w = 1
  #  dbj$wn = dbj$w / sum(dbj$w)
  #  db_gp$w = 1
  #  db_gp$wn = 1
  #  dbjp$w = 1
  #  dbjp$wn = dbjp$w / sum(dbjp$w)
  #
  #
  #
  #
  #  dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
  #  dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)
  #
  #  dbp$wn <- dbp$w / sum(dbp$w)
  #
  #  dbp_p = dbp
  #
  #
  #
  #
  #  for (f in unique(db_gp[[phenotype_var]])) {
  #    dbp_p[[f]] = dbp_p[[f]] / sum(dbp[[f]])
  #  }
  #
  #  for (f in unique(db_gp[[phenotype_var]])) {
  #    ix = dbjp[[indVar]][ dbjp[[phenotype_var]] == f]
  #    if (length(ix) == 0) {
  #      next
  #    }
  #    pilist = list()
  #    pilist2 = list()
  #    ilist = list()
  #    pnlist = list()
  #    pcslist = list()
  #    pslist = list()
  #    wlist = list()
  #    wnlist = list()
  #    pplist = list()
  #    for (i in ix) {
  #      w = dbj$w[dbj$ind==i]
  #      wn = dbj$wn[dbj$ind==i]
  #
  #      ixx = ((dbjp[[indVar]]==i) & (dbjp[[phenotype_var]]==f))
  #      pii = dbjp$p[ixx] * sum(dbjp$p[!ixx]) ## pi * pj
  #
  #      ## pi and pj within phenotype.
  #      piii = dbp_p[[f]][dbp_p$ind==i] * sum(dbp_p[[f]][dbp_p$ind!=i]) ## pi * pj
  #
  #      #ixp = ((dbjp[[indVar]]!=i) & (dbjp[[phenotype_var]]==f)) ##
  #
  #      # sum of affinities for each cell
  #      #pii = dpp[dpp[[indVar]]==i,][[f]] * sum(dpp[dpp[[indVar]]!=i,][[f]])
  #      if (w > 0) {
  #        pilist[i] = pii * w
  #        ilist[i] = pii * (1/wn)
  #        pnlist[i] = pii * wn
  #        wnlist[[i]] = wn
  #        wlist[[i]] = w
  #        pplist[i] = piii * wn
  #      } else if (w == 0) {
  #        pilist[i] = pii
  #        ilist[i] = pii
  #        pnlist[i] = pii
  #        pplist[i] = piii
  #      }
  #
  #      #pilist[i] = pii * (1 / dbjp$w[ixx])  # pi*pj*wij over i and all js
  #      #pilist[i] = sum( prob_mat[i, ]  * sum(aff_mtx[i, ])) # pi*pj*wij over i and all js
  #      # pslist[i] = dbjp$w[ixx]
  #      pcslist[i] = pii
  #    }
  #    dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
  #    dgc$intra_clonotypic_entropy_simn[dgc$pheno==f] =  sum(unlist(pnlist))
  #    dgc$intra_clonotypic_entropy_simi[dgc$pheno==f] =  sum(unlist(ilist))
  #    dgc$intra_clonotypic_entropy_simn_private[dgc$pheno==f] =  sum(unlist(pplist))
  #    #dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
  #    #dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
  #    dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
  #    dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
  #    dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
  #    dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
  #    dgc$simpson_index[dgc$pheno==f] = sum(unlist(pcslist))
  #
  #  }
  #  #dgc = db_gp %>% dplyr::group_by({{ ind }}, {{ phenotype_var }}) %>%
  #  #  dplyr::summarise(n = n())
  #
  #  #dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))
  #
  #
  #  #dgc$intra_clonotypic_entropy_sim = 0
  #
  #  return_list = list("affinity_mat" = NULL,
  #                     "db_clone" = db_gp,
  #                     # "g" = NULL,
  #                     "db_pheno"= dgc)
  #
  #  return(dgc)
  #}
  #else{

    # find corresponding unique germs and junctions
    seqs_unq <- df$seqs_unq
    germs_unq <- df$germs_unq
    juncs_unq <- df$juncs_unq
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

    if (all(disim_mtx == dist_mtx) | any(rowSums(disim_mtx) == 0)) {
      # get required info based on the method

      db_gp$ind <- 0
      db_gp$n_clone = n_unq

      for (i in 1:n_unq) {
        #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
        db_gp$ind[ind_unq[[i]]] <- i
      }

      ind = "ind"
      dbjp = bcrCounts_pheno(db_gp, inds = "ind", pheno=phenotype_var)


      dbj = bcrCounts(db_gp, inds = "ind")

      dbj$w = 1
      dbj$wn = dbj$w / sum(dbj$w)
      db_gp$w = 1
      db_gp$wn = 1
      dbjp$w = 1
      dbjp$wn = dbjp$w / sum(dbjp$w)




      dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
      dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)

      dbp$wn <- dbp$w / sum(dbp$w)

      dbp_p = dbp




      for (f in unique(db_gp[[phenotype_var]])) {
        dbp_p[[f]] = dbp_p[[f]] / sum(dbp[[f]])
      }

      for (f in unique(db_gp[[phenotype_var]])) {
        ix = dbjp[[indVar]][ dbjp[[phenotype_var]] == f]
        if (length(ix) == 0) {
          next
        }
        pilist = list()
        pilist2 = list()
        ilist = list()
        pnlist = list()
        pcslist = list()
        pslist = list()
        wlist = list()
        wnlist = list()
        pplist = list()
        for (i in ix) {
          w = dbj$w[dbj$ind==i]
          wn = dbj$wn[dbj$ind==i]

          ixx = ((dbjp[[indVar]]==i) & (dbjp[[phenotype_var]]==f))
          pii = dbjp$p[ixx] * sum(dbjp$p[!ixx]) ## pi * pj

          ## pi and pj within phenotype.
          piii = dbp_p[[f]][dbp_p$ind==i] * sum(dbp_p[[f]][dbp_p$ind!=i]) ## pi * pj

          #ixp = ((dbjp[[indVar]]!=i) & (dbjp[[phenotype_var]]==f)) ##

          # sum of affinities for each cell
          #pii = dpp[dpp[[indVar]]==i,][[f]] * sum(dpp[dpp[[indVar]]!=i,][[f]])
          if (w > 0) {
            pilist[i] = pii * w
            ilist[i] = pii * (1/wn)
            pnlist[i] = pii * wn
            wnlist[[i]] = wn
            wlist[[i]] = w
            pplist[i] = piii * wn
          } else if (w == 0) {
            pilist[i] = pii
            ilist[i] = pii
            pnlist[i] = pii
            pplist[i] = piii
          }

          #pilist[i] = pii * (1 / dbjp$w[ixx])  # pi*pj*wij over i and all js
          #pilist[i] = sum( prob_mat[i, ]  * sum(aff_mtx[i, ])) # pi*pj*wij over i and all js
          # pslist[i] = dbjp$w[ixx]
          pcslist[i] = pii
        }
        dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
        dgc$intra_clonotypic_entropy_simn[dgc$pheno==f] =  sum(unlist(pnlist))
        dgc$intra_clonotypic_entropy_simi[dgc$pheno==f] =  sum(unlist(ilist))
        dgc$intra_clonotypic_entropy_simn_private[dgc$pheno==f] =  sum(unlist(pplist))
        #dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
        #dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
        dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
        dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
        dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
        dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
        dgc$simpson_index[dgc$pheno==f] = sum(unlist(pcslist))

      }
      #dgc = db_gp %>% dplyr::group_by({{ ind }}, {{ phenotype_var }}) %>%
      #  dplyr::summarise(n = n())

      #dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))


      #dgc$intra_clonotypic_entropy_sim = 0

      return_list = list("affinity_mat" = NULL,
                         "db_clone" = db_gp,
                         # "g" = NULL,
                         "db_pheno"= dgc)

      return(dgc)

    }

    mtx = disim_mtx

    nearest_dist <- apply(mtx, 2,  function(x) {
      gt0 <- which(x > 0)
      if (length(gt0) != 0) { min(x[gt0]) } else { NA }
    })


    krnl_mtx <- krnlMtxGenerator(mtx = mtx)

    aff_mtx <- makeAffinity(mtx_o = mtx,
                            mtx_k = krnl_mtx,
                            thd = max(nearest_dist, na.rm = T))



    threshold = NULL
    base_sim = 0.95
    iter_max = 1000
    nstart = 1000
    #{
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])


    #aff_mtx <- krnl_mtx

    aff_mtx[is.na(aff_mtx)] <- 0

    ## return clone db with unique seq identifier
    db_gp$ind <- 0

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }

    # aff_mtx_c <- aff_mtx[db_gp$ind,db_gp$ind]
    #
    #
    # #diag(aff_mtx) <- 0
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE, diag = FALSE)
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx[db_gp$ind,db_gp$ind], mode = "undirected", weighted = TRUE)
    # g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
    #
    #
    # #g =  igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
    # g =  igraph::simplify(g)
    #
    # igraph::V(g)$k = igraph::degree(g)
    #
    # #v = igraph::eigen_centrality(g)$vector
    # #db_gp$ev = v[ db_gp$ind]
    #
    # igraph::V(g)$Di <-  igraph::diversity(g)
    # igraph::V(g)$Hi =  igraph::V(g)$Di * log10(igraph::V(g)$k)
    #

    # for (i in 1:length(g)) {
    #   #spij = sum(E(g)[.from(i)]$weight / log(V(g)$k[i]))
    #   #spij = sum(g[i,] / log(V(g)$k[i]))
    #   V(g)$sumPij[i] = sum( g[i,] / log(V(g)$k[i]))
    #   #set_vertex_attr(g, index = V(g)[i], name = "sumPij", value = spij)
    #   #V(g)$sumPij[i] =  sum(g[i,] / log(V(g)$k[i]))
    # }

    # browser()
    diag(aff_mtx) = 0

    dbj = bcrCounts(db_gp, inds = "ind")

    aff_mtx = makeCloneMat(aff_mtx, db_gp)

    dbjp = bcrCounts_pheno(db_gp, inds = "ind", pheno=phenotype_var)

    dbj$w = 0
    dbj$wn = 0
    db_gp$w = 0
    db_gp$wn = 0
    dbjp$w = 0
    dbjp$wn = 0


    for (i in dbj$ind) {
      dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
      dbj$wn[dbj$ind==i] = sum(aff_mtx[i, ])/sum(aff_mtx) # s
      dbjp$w[dbjp$ind==i] = sum(aff_mtx[i, ])
      dbjp$wn[dbjp$ind==i] = sum(aff_mtx[i, ]) / sum(aff_mtx) # sum of affinities for each cell)
      db_gp$w[db_gp$ind == i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
      db_gp$wn[db_gp$ind == i] = sum(aff_mtx[i, ]) / sum(aff_mtx)
    }


    ## normalize 0-1

    #dbj$w <- rangeAtoB(dbj$w, 0, 1)
    #dbj$wn <- dbj$w / sum(dbj$w)
    #dbjp$wn <- dbjp$w / sum(dbjp$w)

    #db_gp$w <- rangeAtoB(db_gp$w, 0, 1)
    #dbjp$w <- rangeAtoB(dbjp$w, 0, 1)
    #diag(aff_mtx) = 1


    ## Pi * Pj * Sij, where p is the probability of a cell with BCR sequence x having that sequence
    ## Useful to correct for cases where more than 1 cell has the same BCR sequence
    smat = prob_mat * aff_mtx

    dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
    dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)
    dbp$wn <- dbp$w / sum(dbp$w)





    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp$ind[ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #     #dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = 0
    #   }
    #   dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = sum(smat[ix,])
    # }


    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp[[indVar]][ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   for (i in ix) {
    #     # for (j in 1:nrow(dbj)) {
    #     pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ])
    #   }
    #   dgc$intra_clonotypic_entropy_sim2[dgc$pheno==f] =  sum(unlist(pilist))
    # }

    # browser()

    dbp_p = dbp


    ###proportion private to phenotype
    for (f in unique(db_gp[[phenotype_var]])) {
      dbp_p[[f]] = dbp_p[[f]] / sum(dbp[[f]])
    }
    ###

    for (f in unique(db_gp[[phenotype_var]])) {
      ix = dbjp[[indVar]][ dbjp[[phenotype_var]] == f]
      if (length(ix) == 0) {
        next
      }
      pilist = list()
      pilist2 = list()
      ilist = list()
      pnlist = list()
      pcslist = list()
      pslist = list()
      wlist = list()
      wnlist = list()
      pplist = list()
      for (i in ix) {
        #for (j in 1:nrow(dbj)) {
        #dbp[[f]][ix] * dbj$p[j]
        w = dbj$w[dbj$ind==i]
        wn = dbj$wn[dbj$ind==i]

        ixx = ((dbjp[[indVar]]==i) & (dbjp[[phenotype_var]]==f))
        pii = dbjp$p[ixx] * sum(dbjp$p[!ixx]) ## pi * pj

        ## pi and pj within phenotype.
        piii = dbp_p[[f]][dbp_p$ind==i] * sum(dbp_p[[f]][dbp_p$ind!=i]) ## pi * pj

        #ixp = ((dbjp[[indVar]]!=i) & (dbjp[[phenotype_var]]==f)) ##

        # sum of affinities for each cell
        #pii = dpp[dpp[[indVar]]==i,][[f]] * sum(dpp[dpp[[indVar]]!=i,][[f]])
        if (w > 0) {
          pilist[i] = pii * w
          ilist[i] = pii * (1/wn)
          pnlist[i] = pii * wn
          wnlist[[i]] = wn
          wlist[[i]] = w
          pplist[i] = piii * wn
        } else if (w == 0) {
          pilist[i] = pii
          ilist[i] = pii
          pnlist[i] = pii
          pplist[i] = piii
        }

        #pilist[i] = pii * (1 / dbjp$w[ixx])  # pi*pj*wij over i and all js
        #pilist[i] = sum( prob_mat[i, ]  * sum(aff_mtx[i, ])) # pi*pj*wij over i and all js
        # pslist[i] = dbjp$w[ixx]
        pcslist[i] = pii
      }
      dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
      dgc$intra_clonotypic_entropy_simn[dgc$pheno==f] =  sum(unlist(pnlist))
      dgc$intra_clonotypic_entropy_simi[dgc$pheno==f] =  sum(unlist(ilist))
      dgc$intra_clonotypic_entropy_seqwn[dgc$pheno==f] =  sum(unlist(wnlist))
      dgc$intra_clonotypic_entropy_seqw[dgc$pheno==f] =  sum(unlist(wlist))
      dgc$intra_clonotypic_entropy_simn_private[dgc$pheno==f] =  sum(unlist(pplist))
      #dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
      #dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
      dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
      dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
      dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
      dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
      dgc$simpson_index[dgc$pheno==f] = sum(unlist(pcslist))
      dgc$type[dgc$pheno==f] = "full"

    }


    return_list <- list("affinity_mat" = aff_mtx,
                        "db_clone" = db_gp,
                        #"g" = g,
                        "db_pheno"= dgc)

    return(dgc)
}









#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#'
#' @param db db
#' @param cloneID ID of clone to analyze
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name of clone variable in db
#' @export
#' @returns data.frame of diversity metrics and phenotype ##list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using Simposon diversity
getRGSw <- function(db, cloneID=NULL, phenotype_var="subset", cell_id=NULL, clone = "clone_id", cdr3=FALSE) {
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  normalize = "len" ## c("len", "none"),
  germline = "germline_alignment"
  sequence = "sequence_alignment"
  junction = "junction"
  v_call = "v_call"
  j_call = "j_call"
  fields = NULL
  locus = "locus"
  only_heavy = TRUE
  split_light = FALSE
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

  # get clone
  print(cloneID)
  db_clone <- as.data.frame(db[db$clone_id == cloneID, ])
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
  dgc$rgsw_d = 0
  dgc$rgsw_dw = 0
  dgc$rgsw_norm = 0
  dgc$rgsw_d_norm = 0
  dgc$rgsw_dw_norm = 0
  dgc$type = "limited"
  dgc$wmax = 0
  dgc$beta1 = 0

  if (n_unq <= 2) {

    return(dgc)

  } else{
    dgc$type = "full"

    # find corresponding unique germs and junctions
    seqs_unq <- df$seqs_unq
    germs_unq <- df$germs_unq
    juncs_unq <- df$juncs_unq
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
    aff_mtx <- krnl_mtx
    ### use only for global
    # aff_mtx <- makeAffinity(mtx_o = mtx,
    #                       mtx_k = krnl_mtx,
    #                      thd = max(nearest_dist, na.rm = T))



    threshold = NULL
    base_sim = 0.95
    iter_max = 1000
    nstart = 1000
    #{
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])


    #aff_mtx <- krnl_mtx

    aff_mtx[is.na(aff_mtx)] <- 0

    ## return clone db with unique seq identifier
    db_gp$ind <- 0

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }

    diag(aff_mtx) = 1

    dbj = bcrCounts(db_gp, inds = "ind")
    prob_mat =  aff_mtx
    prob_mat[] <- 0


    dbjp = bcrCounts_pheno(db_gp, inds = "ind", pheno=phenotype_var)

    dbj$w = 0
    dbj$d = 0
    dbj$dw = 0
    dbj$wnv = 0
    db_gp$w = 0
    db_gp$d = 0
    db_gp$dw = 0
    dbjp$w = 0
    dbjp$d = 0
    dbjp$dw = 0


    kdmatrix = 1 - aff_mtx

    for (i in dbj$ind) {
      #dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
      #dbj$wn[dbj$ind==i] = sum(aff_mtx[i, ])/sum(aff_mtx) # s
      dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
      dbj$d[dbj$ind==i] = sum(mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
      dbj$dw[dbj$ind==i] = sum(kdmatrix[i, ]) * nrow(dbj) # sum of affinities for each cell

      dbjp$w[dbjp$ind==i] = sum(aff_mtx[i, ]) * nrow(dbjp)
      dbjp$d[dbjp$ind==i] = sum(mtx[i, ]) * nrow(dbjp) # sum of affinities for each cell
      dbjp$dw[dbjp$ind==i] = sum(kdmatrix[i, ]) * nrow(dbj) # sum of affinities for each cell
      #dbjp$wn[dbjp$ind==i] = sum(aff_mtx[i, ]) / sum(aff_mtx) # sum of affinities for each cell)
      db_gp$w[db_gp$ind == i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
      db_gp$d[db_gp$ind == i] = sum(mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
      db_gp$dw[db_gp$ind == i] = sum(kdmatrix[i, ]) * nrow(dbj) # sum of affinities for each cell
      #db_gp$wn[db_gp$ind == i] = sum(aff_mtx[i, ]) / sum(aff_mtx)
    }


    wijmat = matrix(0, nrow = nrow(dbj), ncol = nrow(dbj))

    for (i in dbj$ind) {
      for (j in dbj$ind) {
        wijmat[i, j] = dbj$p[i] * dbj$p[j]
      }
#dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) # sum of affinities for each cell
#dbj$wn[dbj$ind==i] = sum(aff_mtx[i, ])/sum(aff_mtx) # s
dbj$w[dbj$ind==i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
dbj$d[dbj$ind==i] = sum(mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
dbj$dw[dbj$ind==i] = sum(kdmatrix[i, ]) * nrow(dbj) # sum of affinities for each cell

dbjp$w[dbjp$ind==i] = sum(aff_mtx[i, ]) * nrow(dbjp)
dbjp$d[dbjp$ind==i] = sum(mtx[i, ]) * nrow(dbjp) # sum of affinities for each cell
dbjp$dw[dbjp$ind==i] = sum(kdmatrix[i, ]) * nrow(dbj) # sum of affinities for each cell
#dbjp$wn[dbjp$ind==i] = sum(aff_mtx[i, ]) / sum(aff_mtx) # sum of affinities for each cell)
db_gp$w[db_gp$ind == i] = sum(aff_mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
db_gp$d[db_gp$ind == i] = sum(mtx[i, ]) * nrow(dbj) # sum of affinities for each cell
db_gp$dw[db_gp$ind == i] = sum(kdmatrix[i, ]) * nrow(dbj) # sum of affinities for each cell
#db_gp$wn[db_gp$ind == i] = sum(aff_mtx[i, ]) / sum(aff_mtx)
  }


    ## normalize 0-1

    #dbj$w <- rangeAtoB(dbj$w, 0, 1)
    #dbj$wn <- dbj$w / sum(dbj$w)
    #dbjp$wn <- dbjp$w / sum(dbjp$w)

    #db_gp$w <- rangeAtoB(db_gp$w, 0, 1)
    #dbjp$w <- rangeAtoB(dbjp$w, 0, 1)
    #diag(aff_mtx) = 1


    ## Pi * Pj * Sij, where p is the probability of a cell with BCR sequence x having that sequence
    ## Useful to correct for cases where more than 1 cell has the same BCR sequence
    smat = prob_mat * aff_mtx

    dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
    dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)
    #dbp$wn <- dbp$w / sum(dbp$w)
    #
    dbp = db_gp %>% dplyr::group_by(.data[[indVar]], .data[["w"]], .data[["d"]], .data[["dw"]],  .data[[phenotype_var]]) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::spread(key = {{phenotype_var}}, value = n, fill=0)

    jdmat <-  db_gp %>% dplyr::group_by(.data[[indVar]], .data[["w"]], .data[[phenotype_var]]) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::spread(key = {{phenotype_var}}, value = n, fill=0) %>%
      dplyr::select(-.data[[indVar]]) %>% as.martix()


    regw <- weighted_rich_gini_simpson_conservation(jdmat,
                                                    sequence_weights = disim_mtx,
                                                    conservation_weights = dbp$w,
                                                    #useRichness= TRUE,
                                                    #normalize_phenos = TRUE,
                                                    #pheno_weights = NULL,
                                                    zero_handling = "ignore")



    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp$ind[ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #     #dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = 0
    #   }
    #   dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = sum(smat[ix,])
    # }


    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp[[indVar]][ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   for (i in ix) {
    #     # for (j in 1:nrow(dbj)) {
    #     pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ])
    #   }
    #   dgc$intra_clonotypic_entropy_sim2[dgc$pheno==f] =  sum(unlist(pilist))
    # }

    # browser()

    dbp_p = dbp


    ###proportion private to phenotype
    for (f in unique(db_gp[[phenotype_var]])) {
      dbp_p[[f]] = dbp_p[[f]] / sum(dbp_p[[f]])
    }
    ###

    for (f in unique(db_gp[[phenotype_var]])) {
      rgslist=list()
      wlist = list()
      dlist = list()
      dwlist = list()
      ix = dbp_p[[indVar]]
      if (sum(dbp_p[[f]]) == 0) {
        next
      }

      for (i in ix) {
        jx = setdiff(ix, i)
        w = dbp_p$w[i]
        d = dbp_p$d[i]
        dw = dbp_p$dw[i]
        pj = sum(dbp_p[[f]][jx])
        rgslist[i] =  dbp_p[[f]][i] * (1-dbp_p[[f]][i]) * w
        wlist[i] = w
        dlist[i] =  dbp_p[[f]][i] * (1-dbp_p[[f]][i]) * d
        dwlist[i] =  dbp_p[[f]][i] * (1-dbp_p[[f]][i]) * dw
      }

      if (sum(unlist(rgslist)) == 0) {
        next
      }

      dgc$rgsw[dgc$pheno==f] =  sum(unlist(rgslist))
      dgc$wmax[dgc$pheno==f] = max(unlist(wlist))

      dgc$rgsw_d[dgc$pheno==f] =  sum(unlist(dlist))
      dgc$rgsw_dw[dgc$pheno==f] =  sum(unlist(dwlist))
    }

    dgc$beta1 = dgc$wmax * (1 - (1/n_unq))
    dgc$dmax = max(dbp_p$d)
    dgc$dwmax = max(dbp_p$dw)
    dgc$rgsw_norm = dgc$rgsw / dgc$beta1

    dgc$rgsw_d_norm = dgc$rgsw_d / dgc$beta1
    dgc$rgsw_dw_norm = dgc$rgsw_dw / (dgc$dwmax * (1 - (1/n_unq)))
    dgc$rgsw_d_norm = dgc$rgsw_dw / (dgc$dmax * (1 - (1/n_unq)))



    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   pilist2 = list()
    #   ilist = list()
    #   pnlist = list()
    #   pcslist = list()
    #   pslist = list()
    #   wlist = list()
    #   wnlist = list()
    #   pplist = list()
    # }
    #
    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbjp[[indVar]][ dbjp[[phenotype_var]] == f]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   rgslist=list()
    #   wnvlist=list()
    #   pilist = list()
    #   pilist2 = list()
    #   ilist = list()
    #   pnlist = list()
    #   pcslist = list()
    #   pslist = list()
    #   wlist = list()
    #   wnlist = list()
    #   pplist = list()
    #   for (i in ix) {
    #     #for (j in 1:nrow(dbj)) {
    #     #dbp[[f]][ix] * dbj$p[j]
    #     w = dbj$w[dbj$ind==i]
    #     #wn = dbj$wn[dbj$ind==i]
    #     #wnv = dbj$wnv[dbj$ind==i]
    #
    #     ixx = ((dbjp[[indVar]]==i) & (dbjp[[phenotype_var]]==f))
    #     pii = dbjp$p[ixx] * sum(dbjp$p[!ixx]) ## pi * pj
    #
    #     ## pi and pj within phenotype.
    #     piii = dbp_p[[f]][dbp_p$ind==i] * sum(dbp_p[[f]][dbp_p$ind!=i]) ## pi * pj
    #
    #     #ixp = ((dbjp[[indVar]]!=i) & (dbjp[[phenotype_var]]==f)) ##
    #
    #     # sum of affinities for each cell
    #     #pii = dpp[dpp[[indVar]]==i,][[f]] * sum(dpp[dpp[[indVar]]!=i,][[f]])
    #     if (w > 0) {
    #       pilist[i] = pii * w
    #       ilist[i] = pii * (1/wn)
    #       pnlist[i] = pii * wn
    #       wnlist[[i]] = wn
    #       wlist[[i]] = w
    #       pplist[i] = piii * wn
    #       wnvlist[[i]] = wnv
    #       rgslist[i] = pii * wnv
    #     } else if (w == 0) {
    #       pilist[i] = pii
    #       ilist[i] = pii
    #       pnlist[i] = pii
    #       pplist[i] = piii
    #       wnvlist[[i]] = 0
    #     }
    #
    #     #pilist[i] = pii * (1 / dbjp$w[ixx])  # pi*pj*wij over i and all js
    #     #pilist[i] = sum( prob_mat[i, ]  * sum(aff_mtx[i, ])) # pi*pj*wij over i and all js
    #     # pslist[i] = dbjp$w[ixx]
    #     pcslist[i] = pii
    #   }
    # dgc$rgsw[dgc$pheno==f] =  sum(unlist(rgslist))
    # dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
    # dgc$intra_clonotypic_entropy_simn[dgc$pheno==f] =  sum(unlist(pnlist))
    # dgc$intra_clonotypic_entropy_simi[dgc$pheno==f] =  sum(unlist(ilist))
    # dgc$intra_clonotypic_entropy_seqwn[dgc$pheno==f] =  sum(unlist(wnlist))
    # dgc$intra_clonotypic_entropy_seqw[dgc$pheno==f] =  sum(unlist(wlist))
    # dgc$intra_clonotypic_entropy_simn_private[dgc$pheno==f] =  sum(unlist(pplist))
    # #dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  1 - sum(unlist(pclist))
    # #dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
    # dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
    # dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
    # dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
    # dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
    # dgc$simpson_index[dgc$pheno==f] = sum(unlist(pcslist))
    # dgc$type[dgc$pheno==f] = "full"
    # dgc$wmax[dgc$pheno==f] = max(unlist(wnvlist))

    #}


    # return_list <- list("affinity_mat" = aff_mtx,
    #                     "db_clone" = db_gp,
    #                     #"g" = g,
    #                     "db_pheno"= dgc)

    setDT(dgc)
    dgc <- dgc[rgsw >0,]
    return(dgc)
}
  }










#' Takes a db, a cloneID, and the name of a phenotype variable and returns the affinity matrix, db, and data.frame of per-phenotype diversity metrics for that clone
#' This variant looks at the global abundance of the clone
#' @param db db
#' @param cloneID ID of clone to analyze
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name of clone variable in db
#' @export
#' @returns list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using Simposon diversity
intraclonal_simpson_global_p <- function(db, cloneID=NULL, phenotype_var="subset", cell_id=NULL, clone = "clone_id", cdr3=FALSE) {
  model = "spectral"
  method = "vj"
  linkage = c("single", "average", "complete")
  normalize = "len" ## c("len", "none"),
  germline = "germline_alignment"
  sequence = "sequence_alignment"
  junction = "junction"
  v_call = "v_call"
  j_call = "j_call"
  fields = NULL
  locus = "locus"
  only_heavy = TRUE
  split_light = FALSE
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

  # get clone
  print(cloneID)
  db_clone <- as.data.frame(db[db$clone_id == cloneID, ])
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
  setDT(db_gp)

  mutabs <- shazam::HH_S5F@mutability
  # Generated by using Rcpp::compileAttributes() -> do not edit by hand
  # Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

  #Rcpp::sourceCpp("~/Projects/code/scratch/scoper/src/RcppMutation.cpp") # scoper/src/RcppMutation.cpp")
  #Rcpp::sourceCpp("./src/RcppMutation.cpp")

  n <- nrow(db_gp)

  ### cloning
  # if (method == "vj") {
  #   ### check targeting model
  #   if (!is.null(mutabs)) {
  #     mutabs <- mutabs@mutability
  #   } else {
  #     mutabs <- NULL
  #   }
  # get required info based on the method
  germs <- db_gp[[germline]]
  seqs <- db_gp[[sequence]]
  ## accomdate different length sequences for hamming distance
  juncs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
  juncs <- as.character( set_eq_seqDistance(juncs))
  junc_length <- unique(stringi::stri_length(juncs))
  # find unique seqs
  seqs <- paste(seqs, juncs, germs, sep = "|")
  df <- data.table::as.data.table(seqs)[, list(list(.I)), by=seqs] %>%
    tidyr::separate(col = seqs, into = c("seqs_unq", "juncs_unq", "germs_unq"), sep = "\\|")
  n_unq <- nrow(df)
  ind_unq <- df$V1

  ## make result df here
  dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))
  dgc[[phenotype_var]] = unique(db_gp[[phenotype_var]])
  #dgc[[phenotype_var]] = unique(db_gp[[phenotype_var]])
  dgc$intra_clonotypic_entropy_sim = 0
  dgc$intra_clonotypic_entropy_sim_cells = 0
  dgc$intra_clonotypic_entropy_sim_seq = 0
  dgc$simpson_index = 0
  dgc$clone_id = cloneID
  dgc$clone_size = n

  if (n_unq <= 2) {
    db_gp$ind <- 0
    db_gp$n_clone = n_unq

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }


    #db_gp$vertex_diversity = 0
    #db_gp$gll = 0

    ind = "ind"

    #dgc = db_gp %>% dplyr::group_by({{ ind }}, {{ phenotype_var }}) %>%
    #  dplyr::summarise(n = n())

    #dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))


    #dgc$intra_clonotypic_entropy_sim = 0

    return_list = list("affinity_mat" = NULL,
                       "db_clone" = db_gp,
                       # "g" = NULL,
                       "db_pheno"= dgc)

    return(dgc)
  } else{

    # find corresponding unique germs and junctions
    seqs_unq <- df$seqs_unq
    germs_unq <- df$germs_unq
    juncs_unq <- df$juncs_unq
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

    if (all(disim_mtx == dist_mtx) | any(rowSums(disim_mtx) == 0)) {
      # get required info based on the method
      seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction)]]
      junc_length <- unique(stringi::stri_length(seqs))
      # find unique seqs
      df <- data.table::as.data.table(seqs)[, list(list(.I)), by = seqs]
      n_unq <- nrow(df)
      ind_unq <- df$V1
      seqs_unq <- df$seqs
      if (n_unq == 1) {
        return(gdc)
        #return(list("idCluster" = rep(1, n),
        #            "n_cluster" = 1,
        #            "eigen_vals" = rep(0, n)))
      }
      # calculate unique seuences distance matrix
      disim_mtx <- alakazam::pairwiseDist(seq = seqs_unq,
                                          dist_mat = getDNAMatrix(gap = 0))
    }
    mtx = disim_mtx
    threshold = NULL
    base_sim = 0.95
    iter_max = 1000
    nstart = 1000
    mtx = disim_mtx
    #{
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])

    krnl_mtx <- krnlMtxGenerator(mtx = mtx)
    aff_mtx <- krnl_mtx

    aff_mtx[is.na(aff_mtx)] <- 0

    ## return clone db with unique seq identifier
    db_gp$ind <- 0

    for (i in 1:n_unq) {
      #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
      db_gp$ind[ind_unq[[i]]] <- i
    }

    # aff_mtx_c <- aff_mtx[db_gp$ind,db_gp$ind]
    #
    #
    # #diag(aff_mtx) <- 0
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE, diag = FALSE)
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx[db_gp$ind,db_gp$ind], mode = "undirected", weighted = TRUE)
    # g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
    # #g = igraph::graph_from_adjacency_matrix(aff_mtx, mode = "undirected", weighted = TRUE)
    #
    #
    # #g =  igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
    # g =  igraph::simplify(g)
    #
    # igraph::V(g)$k = igraph::degree(g)
    #
    # #v = igraph::eigen_centrality(g)$vector
    # #db_gp$ev = v[ db_gp$ind]
    #
    # igraph::V(g)$Di <-  igraph::diversity(g)
    # igraph::V(g)$Hi =  igraph::V(g)$Di * log10(igraph::V(g)$k)
    #

    # for (i in 1:length(g)) {
    #   #spij = sum(E(g)[.from(i)]$weight / log(V(g)$k[i]))
    #   #spij = sum(g[i,] / log(V(g)$k[i]))
    #   V(g)$sumPij[i] = sum( g[i,] / log(V(g)$k[i]))
    #   #set_vertex_attr(g, index = V(g)[i], name = "sumPij", value = spij)
    #   #V(g)$sumPij[i] =  sum(g[i,] / log(V(g)$k[i]))
    # }

    dbj = bcrCounts(db_gp, inds = "ind")
    prob_mat =  aff_mtx
    prob_mat[] <- 0

    for (i in 1:nrow(dbj)) {
      for (j in 1:nrow(dbj)) {
        prob_mat[i, j] = dbj$p[i] * dbj$p[j]
      }
    }


    diag(aff_mtx) = 0

    ## Pi * Pj * Sij, where p is the probability of a cell with BCR sequence x having that sequence
    ## Useful to correct for cases where more than 1 cell has the same BCR sequence
    smat = prob_mat * aff_mtx

    dbp =  get_jd(db_gp, pheno = phenotype_var,indVar = indVar)
    dpp =  get_jd_p(db_gp, pheno = phenotype_var,indVar = indVar)




    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp$ind[ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #     #dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = 0
    #   }
    #   dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] = sum(smat[ix,])
    # }


    # for (f in unique(db_gp[[phenotype_var]])) {
    #   ix = dbp[[indVar]][ dbp[[f]] >= 1]
    #   if (length(ix) == 0) {
    #     next
    #   }
    #   pilist = list()
    #   for (i in ix) {
    #     # for (j in 1:nrow(dbj)) {
    #     pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ])
    #   }
    #   dgc$intra_clonotypic_entropy_sim2[dgc$pheno==f] =  sum(unlist(pilist))
    # }


    for (f in unique(db_gp[[phenotype_var]])) {
      ix = dbp[[indVar]][ dbp[[f]] >=1]
      if (length(ix) == 0) {
        next
      }
      pilist = list()
      pcslist = list()
      pslist = list()
      for (i in ix) {
        #for (j in 1:nrow(dbj)) {
        #dbp[[f]][ix] * dbj$p[j]

        pilist[i] = sum( prob_mat[i, ]  * aff_mtx[i, ])
        pslist[i] = sum(aff_mtx[i,])
        pclist = sum(prob_mat[i,])
      }
      dgc$intra_clonotypic_entropy_sim[dgc$pheno==f] =  sum(unlist(pilist))
      dgc$intra_clonotypic_entropy_sim_cells[dgc$pheno==f] =  sum(unlist(pclist))
      dgc$intra_clonotypic_entropy_sim_seq = sum(unlist(pslist))
      dgc$ncells[dgc$pheno==f] =  sum(dbp[[f]])
      dgc$pcells[dgc$pheno==f] =  sum(dpp[[f]])
      dgc$nseqs[dgc$pheno==f] =  nrow(dbp[ix,])
      dgc$pseqs[dgc$pheno==f] =  nrow(dbp[ix,]) / nrow(dbp)
      dgc$simpson_index[dgc$pheno==f] = sum(dpp[[f]]^2)

    }


    return_list <- list("affinity_mat" = aff_mtx,
                        "db_clone" = db_gp,
                        #"g" = g,
                        "db_pheno"= dgc)

    return(dgc)
  }
}






#' function returns intraclonal diversity metrics using simpson (default) or shannon diversities
#' calculated using a BCR sequence similarity matrix
#' @param db db
#' @param phenotype_var phenotype variable
#' @param cell_id cell id
#' @param clone column name of clone id in db
#' @param cloneID clone to analyze
#' @export
#' @returns list with affinity matrix, db for clone, and data.frame of per-phenotype
#' list with affinity matrix, db for clone, and data.frame of per-phenotype diversity metrics using using Simposon (default) or Shannon diversity calculated using a BCR sequence similarity matrix
intraclonal_diversity <- function(db, cloneID, use_clones=TRUE, method = c("simpson", "shannon"), phenotype_var="subset", cell_id=NULL, clone = "clone_id") {
  ## determine if using provided clone calls or VJL groups
  if (use_clones) {
    #cloneID <- db[db[[clone_id]] == clone, clone_id]
  } else {
    db <- prepare_clone(db = db,
                        junction = "junction",
                        v_call = "v_call",
                        j_call = "j_call",
                        first = FALSE,
                        cdr3 = FALSE,
                        fields = NULL,
                        cell_id = cell_id,
                        locus = "locus",
                        only_heavy = TRUE,
                        mod3 = FALSE,
                        max_n = 0)

    #cloneID <- db[db[[clone_id]] == clone, "vjl_group"]
  }

  if (method == "simpson") {
    return(intraclonal_simposon(db, cloneID, phenotype_var, cell_id, clone))
  } else if (method == "shannon") {
    return(intraclonal_shannon(db, cloneID, phenotype_var, cell_id, clone))
  }
}




#####
#' function to caclulate levenshtein distance matrix from a repertoire db
#' @importFrom Biostrings stringDist
#' @return list of distance matrix and db with indicator variable for unique sequences
getLevD <- function(db, region=c("sequence", "junction"), indVar="ind", germline = "germline_alignment",
                    sequence = "sequence_alignment", junction = "junction",
                    cdr3 = FALSE) {
  germs <- db[[germline]]
  seqs <- db[[sequence]]
  juncs <- db[[ifelse(cdr3, cdr3_col, junction)]]
  junc_length <- unique(stringi::stri_length(juncs))
  # find unique seqs
  seqs <- paste(seqs, juncs, germs, sep = "|")
  df <- data.table::as.data.table(seqs)[, list(list(.I)), by=seqs] %>%
    tidyr::separate(col = seqs, into = c("seqs_unq", "juncs_unq", "germs_unq"), sep = "\\|")
  n_unq <- nrow(df)

  if (region == "junction") {
    seq_set <- df$juncs_unq
  } else if (region == "sequence") {
    seq_set <- df$seqs_unq
  } else {
    stop("region must be either junction or sequence")
  }

  n_unq <- nrow(df)
  ind_unq <- df$V1

  db[[indVar]] <- 0

  for (i in 1:n_unq) {
    #idCluster[ind_unq[[i]]] <- idCluster_unq[i]
    db[[indVar]][ind_unq[[i]]] <- i
  }
  seq_d <- as.matrix(stringDist(seq_set, method="levenshtein"))

  return_list = list("distance_mat" = seq_d,
                     "db" = db)

  return(return_list)
}



#' function calculates the global entropy of the sequences in the db pver marginal phenotypes
global_entropies <- function(db, region="sequence", indVar="ind", phenotype_var = "predictions_scanvi",
                             germline = "germline_alignment", sequence = "sequence_alignment",
                             junction = "junction", v_call = "v_call", j_call = "j_call",
                             fields = NULL, locus = "locus", only_heavy = TRUE,
                             split_light = FALSE, targeting_model = shazam::HH_S5F,
                             len_limit = NULL, first = FALSE, cdr3 = FALSE,
                             mod3 = FALSE, max_n = 0, threshold = 1,
                             base_sim = 0.95, iter_max = 1000,
                             nstart = 1000, nproc = 4, verbose = FALSE,
                             log = NULL) {

  res <- getLevD(db, region=region, indVar=indVar, germline = germline,
                 sequence = sequence, junction = junction,
                 cdr3 = cdr3)

  seq_m <- res$distance_mat
  db_gp <- res$db

  seq_sim <-  krnlMtxGenerator(seq_m)


  #jd = get_jd(db_gp, pheno="predictions_scanvi", indVar=indVar)

  dbj = bcrCounts(db_gp, inds = indVar)



  #
  #
  # pijmat =  seq_sim
  # pijmat[] <- 0
  #
  # for (i in 1:nrow(dpp)) {
  #   for (j in 1:nrow(dbj)) {
  #     pijmat[i, j] = dbj$p[i] * dbj$p[j]
  #   }
  # }
  #
  #
  prob_mat =  seq_sim
  prob_mat[] <- 0

  for (i in 1:nrow(dbj)) {
    for (j in 1:nrow(dbj)) {
      prob_mat[i, j] = dbj$p[i] * dbj$p[j]
    }
  }


  # prob_mat2 =  seq_sim
  # prob_mat2[] <- 0
  # r = nrow(seq_sim)
  #   for (i in 1:nrow(dbj)) {
  #   for (j in 1:nrow(dbj)) {
  #     prob_mat2[i, j] = dbj$p[i] * (1 - dbj$p[j])^r
  #   }
  # }


  smat = prob_mat * seq_sim
  #smat2 = prob_mat2 * seq_sim
  dbp =  phenoprops(db_gp, pheno = phenotype_var,inds = indVar)

  #dpp <- get_jd_p(db_gp, pheno=phenotype_var, indVar = indVar)


  dgc = data.frame(pheno = unique(db_gp[[phenotype_var]]))
  #dgc[[phenotype_var]] = unique(db_gp[[phenotype_var]])
  dgc$global_divertsity = 0
  #dgc$phenotype_probability = 0
  #dgc$global_entropy_gen_simpson = 0


  for (f in unique(db_gp[[phenotype_var]])) {
    ix = dbp[[indVar]][ dbp[[f]] >= 1]
    if (length(ix) == 0) {
      next
    }
    pilist = list()
    for (i in ix) {
      # for (j in 1:nrow(dbj)) {
      pilist[i] = sum( prob_mat[i, ]  * seq_sim[i, ])
    }
    dgc$global_divertsity[dgc$pheno==f] =  sum(unlist(pilist))
  }
  #}

  # for (f in unique(db_gp[[phenotype_var]])) {
  #   ix = dbp[[indVar]][ dbp[[f]] >= 1]
  #   if (length(ix) == 0) {
  #     next
  #   }
  #   dgc$global_divertsity[dgc$pheno==f] = sum(smat[ix,])
  #   dgc$phenotype_diversity[dgc$pheno==f] = sum(prob_mat[ix,])
  #   dgc$global_entropy_gen_simpson[dgc$pheno==f] = sum(smat2[ix,])
  # }
  return(dgc)
}




#' Sequence set distance
#'
#' Calculates the distance between pairs of BCRs based on their aligned sequences.
#' The function assume the sequences are at an even length.
#' If not the function will pad the sequences with to the longest sequence length with Ns.
#'
#' @param    seqs          A character list of the sequences for which the distance is to be calculated.
#' @param    AA                    Logical (FALSE by default). If to calculate the distance based on the amino acid sequences.
#'
#' @return
#' A \code{list} containing a  \code{matrix} of the computed distances between the alleles pairs and a vector with length adjusted sequences.
#'
#' @export
seqDistance <- function(seqs, AA=FALSE) {
  ## check if the input is list.

  if (!is.character(seqs))
    stop("The input germline set is not in a character class.")

  ## check sequences length. If not even pad the sequences to the max length
  seqs <-
    gsub("\\s", "N", format(seqs, width = max(nchar(seqs))))

  ## get the distance matrix
  #seqs = seqs[order(names(seqs))]
  #### change gaps from '.' to '-'
  if(AA){
    seqs <- gsub("X", "-", seqs)
    #### create a dna string set
    seqs <- Biostrings::AAStringSet(seqs)
  }else{
    seqs <- gsub("[.]", "-", seqs)
    #### create a dna string set
    seqs <- Biostrings::DNAStringSet(seqs)
  }

  #### compute the distance between pairs. penalize for gaps
  seq_distance <-
    DECIPHER::DistanceMatrix(
      seqs,
      includeTerminalGaps = FALSE,
      penalizeGapLetterMatches = TRUE,
      verbose = FALSE
    )

  return_list = list("distance_mat" = seq_distance,
                     "sequences" = seqs)
  return(return_list)
}




#' Pads the sequences to the longest sequence length with Ns
#'
#'
#' @param    seqs          A character list of the sequences for which the distance is to be calculated.
#'
#' @return
#' A \code{list} containing a  \code{matrix} of the computed distances between the alleles pairs and a vector with length adjusted sequences.
#'
#' @export
set_eq_seqDistance <- function(seqs) {
  ## check if the input is list.

  if (!is.character(seqs))
    stop("The input germline set is not in a character class.")

  ## check sequences length. If not even pad the sequences to the max length
  seqs <-
    gsub("\\s", "N", format(seqs, width = max(nchar(seqs))))

  ## get the distance matrix
  #seqs = seqs[order(names(seqs))]
  #### change gaps from '.' to '-'
  seqs <- gsub("[.]", "-", seqs)
  #### create a dna string set
  seqs <- Biostrings::DNAStringSet(seqs)

  #### compute the distance between pairs. penalize for gaps

  return(seqs)
}



