calc_kinmat <- function(ped, id_col = 1, sire_col = 2, dam_col = 3) {
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
  ped_reorder <- reorder_ped(ped)

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
  return(kinmat[-(n_plus_one), -(n_plus_one)])
}



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
    kin_mat <- calc_kinmat(ped = ped,
                           id_col = id_col,
                           sire_col = sire_col,
                           dam_col = dam_col)
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
