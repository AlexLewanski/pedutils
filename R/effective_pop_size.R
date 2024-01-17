
calc_Ne_f <- function(f, eqg) {
  delF <- 1 - ((1 - f)^(1/(eqg - 1)))
  Ne <- 1/(2*mean(delF))

  return(
    list(Ne = Ne,
         stderr = (2/sqrt(length(f))) * (Ne^2) * stats::sd(delF)
    )
  )
}

wrapper_Ne_f <- function(id, ped, kin_mat, f_info = NULL, eqg_info) {

  if (is.null(f_info))
    f_info <- inbr_from_kinmat(ped = ped,
                               kin_mat = kin_mat,
                               id_col = 1, sire_col = 2, dam_col = 3)

  return(
    calc_Ne_f(f = f_info$f_ped[match(id, f_info$id)],
              eqg = eqg_info$equiv_gen[match(id, eqg_info$id)])
  )
}


calc_Ne_c <- function(coanc, mean_eqg) {
  delCoanc <- 1 - ((1 - coanc)^(1/mean_eqg))
  Ne <- 1/(2*mean(delCoanc))
  return(list(
    Ne = Ne,
    stderr = (2/sqrt(length(coanc))) * (Ne^2) * stats::sd(delCoanc)
  ))
}


wrapper_Ne_c <- function(id, coanc_mat, eqg_info, max_compar = 100000) {

  id_pairs <- t(utils::combn(id, 2))

  #THIS SAFEGUARD COULD DEFINITELY BE IMPROVED AND MADE MORE EFFICIENT
  #currently all combinations of individuals are generated and then subsampled
  #perhaps a better way to do this would be to take the coancestry matrix, look
  #at the number of values, and then directly subsample values from that?
  if (nrow(id_pairs) > max_compar)
    id_pairs <- id_pairs[sample(nrow(id_pairs), max_compar, replace = FALSE),]

  coanc_vec <- coanc_mat[id_pairs]
  mean_eqg <- (eqg_info$equiv_gen[match(id_pairs[,1], eqg_info$id)] + eqg_info$equiv_gen[match(id_pairs[,2], eqg_info$id)])/2

  return(calc_Ne_c(coanc = coanc_vec, mean_eqg = mean_eqg))
}




#' Calculate effective population size from a pedigree
#'
#' @param ped A pedigree dataframe in the format output by the process_ped function.
#' @param Ne The estimators under which to calculate Ne. Currently, two options are implemented---one based on inbreeding (Ne_f) and another based on coancestry (Ne_c)
#' @param id The individuals to be considered in the Ne calculation. The option "all" (the default) uses all individuals in the pedigree in the calculation.
#' @param coanc_mat The coancestry (kinship) matrix for the individuals in the pedigree. If this is not supplied, the coancestry matrix will be calculated internally.
#' @param inbr A dataframe containing information on inbreeding. The dataframe must include a column named "id" with id of each individual and a column named "f_ped," which records the inbreeding value for each individual. If this is not supplied, the inbreeding information will be calculated internally.
#' @param eqg A dataframe containing information on effective number of generations for each individual. The dataframe must include a column named "id" with id of each individual and a column named "equiv_gen," which records the equivalent number of generations for each individual. If this is not supplied, the number of equivalent generations will be calculated internally.
#' @param max_compar The number coancestry values to be considered in the calculation of Ne. If the number of values is higher than max_compar, then the coancestry values will be randomly subsampled down to max_compar. This argument is only considered for the "Ne_c" estimator of Ne.
#'
#' @return a list containing the Ne estimates and standard errors
#' @export
calc_ped_Ne <- function(ped,
                        Ne = c('Ne_f', 'Ne_c'),
                        id = 'all',
                        coanc_mat = NULL,
                        inbr = NULL,
                        eqg = NULL,
                        max_compar = 100000) {

  #=====================================
  #=== INITIAL PROCESSING AND CHECKS ===
  #=====================================
  if (is.null(coanc_mat))
    coanc_mat <- calc_pedmat(ped = ped, type = 'kinship')

  if (is.null(eqg)) {
    eqg <- ped_summary_stats(ped = ped, summary = 'equiv_gen')
  } else {
    stopifnot(is.data.frame(eqg))
    if (!'equiv_gen' %in% colnames(eqg))
      stop('The equivalent generation info needs to be contained in a column labelled "equiv_gen"')
  }

  if (identical(id, 'all'))
    id <- ped[,'id',drop = TRUE]


  #=========================
  #=== CALCULATION of Ne ===
  #=========================
  ne_list <- list()

  for (i in Ne) {

    ne_list[[i]] <- switch(i,
                           Ne_f = {
                             wrapper_Ne_f(id = id,
                                          ped = ped,
                                          kin_mat = coanc_mat,
                                          eqg_info = eqg)
                           },
                           Ne_c = {
                             wrapper_Ne_c(id = id,
                                          coanc_mat = coanc_mat,
                                          eqg_info = eqg,
                                          max_compar = max_compar)
                           })

  }

  return(ne_list)
}


