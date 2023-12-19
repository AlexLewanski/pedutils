#########################
### COANCESTRY MATRIX ###
#########################
#' Calculate kinship/coancestry matrix or additive relationship matrix from a pedigree
#'
#' @param ped the pedigree from which to calculate the kinship matrix
#' @param id_col the name or index of the id column
#' @param sire_col the name or index of the sire column
#' @param dam_col the name or index of the dam column
#' @param type Whether the kinship matrix or additive relationship matrix should be calculated
#'
#' @return the kinship or additive relationship matrix
#' @export
calc_pedmat <- function(ped,
                        id_col = 1,
                        sire_col = 2,
                        dam_col = 3,
                        type = c('kinship', 'add_rel')) {
  #resources:
  #https://doi.org/10.3389/fgene.2021.655638
  #https://github.com/mayoverse/kinship2/blob/master/R/kinship.R

  #========================
  #== PREPARATION STEPS ===
  #========================

  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  #perhaps this isn't needed anymore because a custom function is now used...
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)

  #reorder pedigree so parents occur before offspring
  #ped_reorder <- ped[order(pedigree::orderPed(ped[,c(id_col, sire_col, dam_col)])),]
  ped_reorder <- reorder_ped(ped,
                             id_col = id_col,
                             sire_col = sire_col,
                             dam_col = dam_col)

  id_vec <- ped_reorder[,id_col, drop = TRUE] #id vector
  n_plus_one <- length(id_vec) + 1 #number of individuals in pedigree plus 1

  #vectors of parents: vectors hold the index of each parent in the id_vec to
  #facilitate easy retrieval of the parents' kinship values from the kinmat
  #if the parent isn't in id_vec, it is a founder so assign it a value of n_plus_one
  sire_vec <- match(ped_reorder[, sire_col, drop = TRUE], id_vec, nomatch = n_plus_one)
  dam_vec <- match(ped_reorder[, dam_col, drop = TRUE], id_vec, nomatch = n_plus_one)

  #empty matrix to store kinship values
  kinmat <- matrix(data = 0,
                   nrow = n_plus_one, ncol = n_plus_one,
                   dimnames = list(c(id_vec, n_plus_one), c(id_vec, n_plus_one))
  )


  #===========================
  #== KINSHIP CALCULATIONS ===
  #===========================
  for (x in seq_along(id_vec)) {
    kinmat[,x] <- kinmat[x,] <- (kinmat[dam_vec[x],] + kinmat[sire_vec[x],])/2
    kinmat[x,x] <- (1 + kinmat[dam_vec[x], sire_vec[x]])/2
  }

  #return the kinship matrix (without the extra founder column and row)
  return(
    switch(type,
           kinship = {kinmat[-(n_plus_one), -(n_plus_one)]},
           add_rel = {2*kinmat[-(n_plus_one), -(n_plus_one)]})
    )
}





#' Calculate partial kinship matrix for a given pedigree and focal founder
#'
#' @param ped the pedigree from which to calculate the partial kinship matrix
#' @param id_col the name or index of the id column
#' @param sire_col the name or index of the sire column
#' @param dam_col the name or index of the dam column
#' @param founder_val the value in the sire and dam columns that represent founder
#' @param focal_founder the founder for which to calculate the partial kinship matrix
#'
#' @return a partial kinship matrix for the focal founder
#' @export
#'
calc_partial_kinmat <- function(ped,
                                id_col = 1,
                                sire_col = 2,
                                dam_col = 3,
                                founder_val = 0,
                                focal_founder) {
  #resources:
  #https://doi.org/10.3389/fgene.2021.655638
  #https://github.com/mayoverse/kinship2/blob/master/R/kinship.R
  #https://doi.org/10.1016/j.livsci.2006.04.007

  #========================
  #== PREPARATION STEPS ===
  #========================

  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)

  #reorder pedigree so parents occur before offspring
  #ped_reorder <- ped[order(pedigree::orderPed(ped[,c(id_col, sire_col, dam_col)])),] #12/13 --> check that this works
  ped_reorder <- reorder_ped(ped,
                             id_col = id_col,
                             sire_col = sire_col,
                             dam_col = dam_col)

  founders <- ped_reorder[ped_reorder[,sire_col] == founder_val & ped_reorder[,dam_col] == founder_val, id_col, drop = TRUE]
  non_focal_founders <- founders[founders != focal_founder]

  id_vec <- ped_reorder[,id_col] #id vector
  n_plus_one <- length(id_vec) + 1 #number of individuals in pedigree plus 1

  #vectors of parents: vectors hold the index of each parent in the id_vec to
  #facilitate easy retrieval of the parents' kinship values from the kinmat
  #if the parent isn't in id_vec, it is a founder so assign it a value of n_plus_one
  sire_vec <- match(ped_reorder[, sire_col], id_vec, nomatch = n_plus_one)
  dam_vec <- match(ped_reorder[, dam_col], id_vec, nomatch = n_plus_one)
  founder_index <- which(id_vec == focal_founder)

  #matrix of 0s to store kinship values
  partial_kinmat <- matrix(data = 0,
                           nrow = n_plus_one, ncol = n_plus_one,
                           dimnames = list(c(id_vec, n_plus_one), c(id_vec, n_plus_one))
  )


  #===================================
  #== PARTIAL KINSHIP CALCULATIONS ===
  #===================================
  for (x in seq_along(id_vec)) {
    if (id_vec[x] %in% non_focal_founders) {
      partial_kinmat[,x] <- partial_kinmat[x,] <- 0
    } else {
      partial_kinmat[,x] <- partial_kinmat[x,] <- (partial_kinmat[dam_vec[x],] + partial_kinmat[sire_vec[x],])/2

      if (x == founder_index) {
        partial_kinmat[x,x] <- 0.5
      } else {
        partial_kinmat[x,x] <- partial_kinmat[x,founder_index] + 0.5*partial_kinmat[dam_vec[x], sire_vec[x]]
      }
    }
  }

  return(partial_kinmat[-(n_plus_one), -(n_plus_one)])
}






#############################################
### INBREEDING BASED ON COANCESTRY MATRIX ###
#############################################

#' Calculate inbreeding from a kinship/coancestry matrix
#'
#' @param ped the pedigree from which to calculate inbreeding
#' @param kin_mat the kinship matrix from to which to calculate inbreeding. If this is not provided, the matrix is calculated internally.
#' @param id_col the name or index of the id column
#' @param sire_col the name or index of the sire column
#' @param dam_col the name or index of the dam column
#'
#' @return a dataframe with inbreeding values for each individual
#' @export
#'
inbr_from_kinmat <- function(ped,
                             kin_mat = NULL,
                             id_col = 1,
                             sire_col = 2,
                             dam_col = 3) {

  #On my computer, feeding a tibble to orderPed results in the R session being
  #aborted. This doesn't happen if ped is a dataframe
  #MAYBE THIS ISN'T NECESSARY ANYMORE
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)

  #if a kinship matrix isn't provided, calculate one from the input pedigree
  if (is.null(kin_mat)) {
    kin_mat <- calc_pedmat(ped = ped,
                           id_col = id_col,
                           sire_col = sire_col,
                           dam_col = dam_col,
                           type = 'kinship')
  }

  if (!identical(colnames(kin_mat), rownames(kin_mat)))
    stop('The column and row names for kin_mat must be identical (including identical order)')

  ### extract the indices in kin_mat for sire and dam
  #(the dam and sire indices for founders will return NA)
  dam_vec <- match(ped[,dam_col, drop = TRUE], colnames(kin_mat))
  sire_vec <- match(ped[,sire_col, drop = TRUE], colnames(kin_mat))

  inbr_df <- data.frame(id = ped[,id_col, drop = TRUE],
                        f_ped = NA)

  for (i in seq_len(nrow(inbr_df))) {

    if ( any(is.na(c(dam_vec[i], sire_vec[i]))) ) {#check to make sure this works and is correct
      inbr_df[i,'f_ped'] <- 0
    } else {
      inbr_df[i,'f_ped'] <- kin_mat[dam_vec[i], sire_vec[i]]
    }
  }

  return(inbr_df)
}



#' Calculate the partial inbreeding value for each individual derived from each founder
#'
#' @param ped the pedigree from which to calculate the partial kinship matrix
#' @param id_col the name or index of the id column
#' @param sire_col the name or index of the sire column
#' @param dam_col the name or index of the dam column
#' @param founder_val the value in the sire and dam columns that represent founder
#'
#' @return a list that contains a list of partial coancestry matrices and a dataframe of partial inbreeding values
#' @import stats
#' @import dplyr
#' @export
#'
partial_inbreeding <- function(ped,
                               id_col = 1,
                               sire_col = 2,
                               dam_col = 3,
                               founder_val = 0) {

  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  #THIS MIGHT NOT BE NECESSARY ANYMORE: 12/13/23
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)

  founders <- ped[ped[,sire_col] == founder_val & ped[,dam_col] == founder_val, id_col, drop = TRUE]

  partial_coanc_list <- lapply(stats::setNames(nm = founders), function(x) {
    calc_partial_kinmat(ped = ped,
                        id_col = id_col,
                        sire_col = sire_col,
                        dam_col = dam_col,
                        focal_founder = x)
  })

  partial_inbr_list <- lapply(partial_coanc_list, function(COANC_MAT, ped) {
    inbr_from_kinmat(ped = ped,
                     kin_mat = COANC_MAT,
                     id_col = id_col,
                     sire_col = sire_col,
                     dam_col = dam_col)
  }, ped = ped)

  #return list containing 2 elements:
  #(1) list of partial coancestry matrices
  #(2) dataframe of partial inbreeding values
  return(
    list(partial_coanc_list = partial_coanc_list,
         partial_inbr = dplyr::bind_rows(partial_inbr_list, .id = "focal_founder"))
  )

}




################################################
### FUNCTIONS BASED ON GENE DROP SIMULATIONS ###
################################################
#' Calculation of inbreeding from gene dropping
#'
#' @param gdrop_output output from the simple_gene_drop function (specifically, the geno_origin dataframe)
#'
#' @return a dataframe containing the inbreeding values
#' @import dplyr
#' @importFrom rlang .data
#' @export
#'
fped_gdrop <- function(gdrop_output) {
  return(
    gdrop_output %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(f_status = .data$sire_geno_origin == .data$dam_geno_origin) %>%
      dplyr::summarize(f_ped = sum(.data$f_status)/dplyr::n())
  )
}


#' Calculation of partial inbreeding from gene dropping
#'
#' @param gdrop_output output from the simple_gene_drop function (specifically, the geno_origin dataframe)
#'
#' @return a dataframe containing each individual's inbreeding values contributed by each founder
#' @import dplyr
#' @importFrom rlang .data
#' @export
#'
partial_founder_fped_gdrop <- function(gdrop_output) {
  return(
    gdrop_output %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(f_status = .data$sire_geno_origin == .data$dam_geno_origin) %>%
      dplyr::mutate(fped_founder_id = ifelse(.data$f_status, .data$sire_geno_origin, NA)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(.data$id, .data$fped_founder_id) %>%
      dplyr::summarize(#id = first(id),
        #founder_id = first(fped_founder_id),
        founder_count = sum(.data$f_status),
        .groups = 'drop') %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(fped_partial_ancestor = .data$founder_count/sum(.data$founder_count)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(.data$fped_founder_id))
  )
}
