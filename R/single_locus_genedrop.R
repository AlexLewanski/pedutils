############################
############################
### GENE DROP SIMULATION ###
############################
############################

#library(pedigree)
#library(tidyverse) #specifically tibble


#################
### FUNCTIONS ###
#################

#' Single locus gene drop on a pedigree
#'
#' @param ped pedigree on which to perform gene drop
#' @param id_col name or index of id column
#' @param sire_col name or index of sire column
#' @param dam_col name or index of dam column
#' @param sex_col name or index of sex column
#' @param sims number of gene drop simulations
#' @param pop_freq the frequency of the allele in the starting population
#' @param fixed_founder_genotypes should founder genotypes be fixed?
#' @param calc_contribution should the ancestry contributions of each founder be recorded?
#' @param report_progress should the simulation progress be shown?
#'
#' @return a list of simulation output
#' @export
#'
simple_gene_drop <- function(ped,
                             id_col,
                             sire_col,
                             dam_col,
                             sex_col,
                             sims = 10,
                             pop_freq = 0.5,
                             fixed_founder_genotypes = FALSE,
                             calc_contribution = FALSE,
                             report_progress = TRUE) {

  #==========================
  #=== INITIAL PROCESSING ===
  #==========================

  ### CHECK INPUTS ###
  report_progress <- isTRUE(report_progress)
  fixed_founder_genotypes <- isTRUE(fixed_founder_genotypes)
  #OTHER CHECKS

  if (report_progress) message('Initial pedigree processing and simulation preparation.')

  ### INITIAL PROCESSING OF PEDIGREE ###
  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)

  # (1) subset down to id, sire, and dam cols
  # (2) rename the columns
  # (3) reorder pedigree so that parents appear before offspring
  # (4) record size of pedigree
  ped_reorder <- ped[,c(id_col, sire_col, dam_col, sex_col)]
  colnames(ped_reorder) <- c('id', 'sire', 'dam', 'sex')

  #ordered_ped <- ped_reorder[order(pedigree::orderPed(ped_reorder)),]
  ordered_ped <- reorder_ped(ped_reorder,
                             id_col = 'id',
                             sire_col = 'sire',
                             dam_col = 'dam')
  ped_size <- nrow(ped_reorder)

  #record the index of sire and dams in new cols (this will make retrieval
  #of sire and dam genotype and genotype origin information easier and faster)
  ordered_ped$sire_id <- match(ordered_ped$sire, ordered_ped$id)
  ordered_ped$dam_id <- match(ordered_ped$dam, ordered_ped$id)


  ### CREATE GENOTYPE AND GENOTYPE ORIGIN MATRICES ###
  geno_mat <- matrix(data = NA, nrow = ped_size, ncol = 2,
                     dimnames = list(ordered_ped$id, c('sire_geno', 'dam_geno')))

  geno_orig_mat <- matrix(data = NA, nrow = ped_size, ncol = 2,
                          dimnames = list(ordered_ped$id, c('sire_geno_origin', 'dam_geno_origin')))

  if (isTRUE(calc_contribution)) {
    contr_mat <- matrix(data = 0, nrow = ped_size, ncol = 2,
                        dimnames = list(ordered_ped$id, c('sire_contr', 'dam_contr')))

    contr_list <- list()
  }


  ### CREATE GENOTYPE INFORMATION FOR FOUNDERS###
  #*** THIS IS AN OBVIOUS PLACE FOR IMPROVEMENT. HOW DO OTHERS DO THIS?

  #indices of founders
  founder_indices <- which(apply(ordered_ped[, c('sire', 'dam')], 1, function(x) all(x == '0')))

  #IF A SINGLE FIXED SET OF FOUNDER ALLELES ACROSS SIMULATIONS
  if (fixed_founder_genotypes) {
    #add genotype information for founders by two independent draws from the binomial
    #geno_mat[founder_indices,] <- t(replicate(length(founder_indices), rbinom(n = 2, size = 1, prob = 0.5), simplify = TRUE))
    geno_mat[founder_indices,] <- matrix(stats::rbinom(n = length(founder_indices)*2, size = 1, prob = pop_freq), ncol = 2)
  }

  #for both alleles in each founder, record the origin as the founder id (not sure if this is right)
  geno_orig_mat[founder_indices,] <- t(sapply(rownames(geno_orig_mat)[founder_indices], function(x) paste(x, 1:2, sep = "_"))) #separate allele identities for each founder
  #geno_orig_mat[founder_indices,] <- rownames(geno_orig_mat)[founder_indices] #if you want to have the same identity for both alleles of each founder

  #=============================
  #=== GENE DROP SIMULATIONS ===
  #=============================

  #initiate progress bar
  if (report_progress) {
    message('Starting gene drop simulation', ifelse(sims > 1, 's.', '.'))
    prog_bar <- utils::txtProgressBar(min = 0, max = sims, initial = 0, char = "*", style = 3)
  }

  gene_drop_list <- list() #list to store each gene drop simulation

  for (i in seq_len(sims)) {

    #IF FOUNDERS SHOULD BE RANDOMLY DRAWN FOR EACH SIMULATION (i.e. if fixed_founder_genotypes is FALSE)
    if (!fixed_founder_genotypes) {
      #add genotype information for founders by two independent draws from the binomial
      geno_mat[founder_indices,] <- matrix(stats::rbinom(n = length(founder_indices)*2,
                                                         size = 1,
                                                         prob = pop_freq),
                                           ncol = 2)
    }

    #transmission probabilities for each allele from each individual's parents
    #for example, the ith element in sire_draw_vec indicates whether the ith individual will inherit
    #allele 1 (sire allele) or allele 2 from its father.
    sire_draw_vec <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L
    dam_draw_vec <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L

    f_vec_mat <- matrix(data = 0, nrow = ped_size, ncol = 2,
                        dimnames = list(ordered_ped$id, c('sire_ibd_count', 'dam_ibd_count')))

    #perform gene drop
    gene_drop_list[[paste0('sim', i)]] <- single_gene_drop(ped_size = ped_size,
                                                           pedigree = ordered_ped,
                                                           sire_draw_vec = sire_draw_vec,
                                                           dam_draw_vec = dam_draw_vec,
                                                           geno_mat = geno_mat,
                                                           geno_orig_mat = geno_orig_mat,
                                                           f_vec_mat = f_vec_mat)

    if (isTRUE(calc_contribution)) {
      contr_list[[paste0('sim', i)]] <- single_genetic_contr(ped_size = ped_size,
                                                             pedigree = ordered_ped,
                                                             sire_draw_vec = sire_draw_vec,
                                                             dam_draw_vec = dam_draw_vec,
                                                             contr_mat = contr_mat)
    }

    if (report_progress) utils::setTxtProgressBar(prog_bar, i)
  }



  #=====================================
  #=== PROCESS AND RETURN SIM OUTPUT ===
  #=====================================
  #THIS CAN BE PRETTY SLOW FOR LARGE NUMBERS OF SIMULATIONS
  #THINK MORE ABOUT HOW TO SPEED THIS UP

  if (report_progress) message('\nProcessing simulation output.')
  output_list <- list()

  output_list[['geno']] <- lapply(gene_drop_list, function(x) {
    tibble::rownames_to_column(as.data.frame(x[['geno_mat']]), var = "id")
  } ) %>%
    dplyr::bind_rows(.id = 'sim')
    #dplyr::bind_rows(., .id = 'sim')

  output_list[['geno_origin']] <- lapply(gene_drop_list, function(x) {
    tibble::rownames_to_column(as.data.frame(x[['geno_orig_mat']]), var = "id")
  }) %>%
    dplyr::bind_rows(.id = 'sim')
    #dplyr::bind_rows(., .id = 'sim')

  output_list[['f_count']] <- lapply(gene_drop_list, function(x) {
    tibble::rownames_to_column(as.data.frame(x[['f_vec_mat']]), var = "id")
  }) %>%
    dplyr::bind_rows(.id = 'sim')
    #dplyr::bind_rows(., .id = 'sim')

  if (isTRUE(calc_contribution)) {
    output_list[['contr_list']] <- contr_list %>%
      dplyr::bind_rows(.id = 'sim')
      #dplyr::bind_rows(., .id = 'sim')
  }

  return(output_list)

}



single_gene_drop <- function(ped_size,
                             pedigree,
                             sire_draw_vec,
                             dam_draw_vec,
                             geno_mat,
                             geno_orig_mat,
                             f_vec_mat) {

  #iterate through every individual and perform the following tasks:
  # (1) randomly draw an allele from each parent
  # (2) copy the origin of the chosen allele

  #notes:
  # - for purposes of speed, the choice of which allele to transmit from
  #   each parent is pre-calculated and fed into the function (sire_draw_vec,
  #   dam_draw_vec)
  # -

  for (i in seq_len(ped_size)) {

    #if an individual doesn't have a sire or dam id, it is a founder
    #founders should already have genotypes so skip to next individual
    #THIS COULD BE IMPROVED (OR CONSIDERED MORE THOUGHTFULLY). THIS
    #WILL PRODUCE WEIRDNESS IF AN INDIVIDUAL HAS A SINGLE KNOWN PARENT
    if (any(is.na(pedigree[i,c('sire_id', 'dam_id')]))) next

    #retrieve sire and dam index for focal individual
    sire_index <- pedigree[i, 'sire_id']
    dam_index <- pedigree[i, 'dam_id']

    #get random index (1 or 2) for which sire and dam ID to transmit
    #sire_draw <- rbinom(n = 1, size = 1, prob = 0.5) + 1L #sum(runif(1, 0, 1) > 0.5) + 1L
    #dam_draw <- rbinom(n = 1, size = 1, prob = 0.5) + 1L #sum(runif(1, 0, 1) > 0.5) + 1L
    sire_draw <- sire_draw_vec[i]
    dam_draw <- dam_draw_vec[i]

    #add chosen sire and dam alleles for the individual's genotype
    geno_mat[i,1] <- geno_mat[sire_index, sire_draw]
    geno_mat[i,2] <- geno_mat[dam_index, dam_draw]

    #record the origin of each allele for the focal individual
    geno_orig_mat[i,1] <- geno_orig_mat[sire_index, sire_draw]
    geno_orig_mat[i,2] <- geno_orig_mat[dam_index, dam_draw]

    sire_ibd <- as.integer(length(unique(geno_orig_mat[sire_index,])) == 1)
    dam_ibd <- as.integer(length(unique(geno_orig_mat[dam_index,])) == 1)
    f_vec_mat[i,1] <- f_vec_mat[sire_index, sire_draw] + sire_ibd
    f_vec_mat[i,2] <- f_vec_mat[dam_index, dam_draw] + dam_ibd

  }

  return(list(geno_mat = geno_mat,
              geno_orig_mat = geno_orig_mat,
              f_vec_mat = f_vec_mat))

}




single_genetic_contr <- function(ped_size,
                                 pedigree,
                                 sire_draw_vec,
                                 dam_draw_vec,
                                 contr_mat) {

  contribute_list <- list()

  for (i in seq_len(ped_size)) {
    contribute_list[[pedigree$id[i]]] <- contr_mat
    contribute_list[[i]][i,] <- 1

    if (i == ped_size) {
      contribute_list[[i]] <- contribute_list[[i]][apply(contribute_list[[i]], 1, function(x) sum(x) != 0),,drop = FALSE] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = 'id')
      break
    }

    for (x in (i + 1):ped_size) {

      if (any(is.na(pedigree[x,c('sire_id', 'dam_id')]))) next

      #retrieve sire and dam index for focal individual
      sire_index <- pedigree[x, 'sire_id']
      dam_index <- pedigree[x, 'dam_id']

      sire_draw <- sire_draw_vec[i]
      dam_draw <- dam_draw_vec[i]

      contribute_list[[i]][x,1] <- contribute_list[[i]][sire_index, sire_draw]
      contribute_list[[i]][x,2] <- contribute_list[[i]][dam_index, dam_draw]
    }

    contribute_list[[i]] <- contribute_list[[i]][apply(contribute_list[[i]], 1, function(x) !all(x == 0)),,drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = 'id')
  }

  return(
    contribute_list %>%
      #dplyr::bind_rows(., .id = 'focal_ind') %>%
      dplyr::bind_rows(.id = 'focal_ind') %>%
      dplyr::mutate(geno_sum = .data$sire_contr + .data$dam_contr) %>%
      dplyr::select(.data$focal_ind, .data$id, .data$geno_sum) %>%
      `rownames<-`( NULL )
  )
}





### ALTERNATIVE VERSION OF GENE DROP THAT USES MATRICES ###


#' Single locus gene drop on a pedigree (using matrices to record the alleles)
#'
#' @param ped pedigree (stored in dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column 3 --> dam. Founder parents should be coded as 0s.
#' @param sims the number of gene drops to perform
#' @param report_progress should the simulation progress be reported?
#'
#' @return a list containing the following items: matrix of alleles derived from the sire (rows = individuals, columns = sims); matrix of alleles derived from the dam (rows = individuals, columns = sims); reorded pedigree; indexed version of reorded pedigree
#' @export
#'
gene_drop_matrix <- function(ped, sims = 10, report_progress = TRUE) {
  ped_reorder <- reorder_ped(ped, 1, 2, 3)
  ped_reorder_index <- index_pedigree(ped_reorder)

  ped_size <- nrow(ped_reorder_index)

  init_allele_mat <- matrix(data = 0, nrow = ped_size, ncol = 2)

  #allele_mat_sire <- matrix(data = 0, nrow = ped_size, ncol = sims)
  #allele_mat_dam <- matrix(data = 0, nrow = ped_size, ncol = sims)

  allele_mat_list <- stats::setNames(replicate(2,
                                               matrix(data = 0,
                                                      nrow = ped_size,
                                                      ncol = sims,
                                                      dimnames = list(ped_reorder[,1,drop=TRUE], NULL)),
                                               simplify = FALSE),
                                     nm = c('sire', 'dam'))

  sire_founder_index <- which(ped_reorder_index[,2,drop=TRUE] == 0)
  dam_founder_index <- which(ped_reorder_index[,3,drop=TRUE] == 0)

  init_allele_mat[sire_founder_index,1] <- sire_founder_index
  init_allele_mat[dam_founder_index,2] <- -dam_founder_index

  #initiate progress bar
  if (report_progress) {
    message('Starting gene drop simulation', ifelse(sims > 1, 's.', '.'))
    prog_bar <- utils::txtProgressBar(min = 0, max = sims, initial = 0, char = "*", style = 3)
  }

  for (SIM in seq_len(sims)) {

    completed_mat <- single_drop_matrix(ped = ped_reorder_index,
                                        ped_size = ped_size,
                                        allele_mat = init_allele_mat)

    for (i in seq_len(2)) allele_mat_list[[i]][,SIM] <- completed_mat[,i]
    #allele_mat_sire[,SIM] <- completed_mat[,1]
    #allele_mat_dam[,SIM] <- completed_mat[,2]

    if (report_progress) utils::setTxtProgressBar(prog_bar, SIM)

  }

  allele_mat_list[['reordered_pedigree']] <- ped_reorder
  allele_mat_list[['indexed_pedigree']] <- ped_reorder_index
  allele_mat_list[['founder_alleles']] <- data.frame(id = c(ped_reorder[sire_founder_index,1,drop=TRUE],
                                                            ped_reorder[dam_founder_index,1,drop=TRUE]),
                                                     allele = c(sire_founder_index,
                                                                -dam_founder_index))
  allele_mat_list[['sim_count']] <- sims

  return(allele_mat_list)
}


single_drop_matrix <- function(ped,
                               ped_size,
                               allele_mat) {

  fid_draw <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L
  mid_draw <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L

  for (i in seq_len(ped_size)) {
    #sire
    if (ped[i,2,drop=TRUE] != 0) {
      allele_mat[i,1] <- allele_mat[ped[i,2,drop=TRUE], fid_draw[i]]
    }

    #dam
    if (ped[i,3,drop=TRUE] != 0) {
      allele_mat[i,2] <- allele_mat[ped[i,3,drop=TRUE], mid_draw[i]]
    }
  }

  return(allele_mat)
}


## @import dplyr
## @import tibble
## @importFrom rlang .data
