#############################
### ALL-PURPOSE FUNCTIONS ###
#############################
harmonic_mean <- function(x) length(x)/sum(1/x)

#######################################
### PEDIGREE MANIPULATON AND CHECKS ###
#######################################
#' Reorder pedigree so that all offspring follow their parents
#'
#' @param ped a pedigree
#' @param id_col name or index of id column
#' @param sire_col name or index of sire column
#' @param dam_col name or index of dam column
#'
#' @return a pedigree ordered so that offspring follow their parents
#' @export
#'
reorder_ped <- function(ped,
                        id_col,
                        sire_col,
                        dam_col) {

  #Zhang et al. 2009
  #An algorithm to sort complex pedigrees chronologically without birthdates

  id_ind <- seq_along(ped[[id_col]])
  id_sire <- match(ped[[sire_col]], ped[[id_col]])
  id_dam <- match(ped[[dam_col]], ped[[id_col]])

  gen_vec <- rep(0, nrow(ped))

  while (TRUE) {
    empty_list <- list()

    for (i in id_ind) {

      if (!is.na(id_sire[i]) && !(id_sire[i] %in% unlist(empty_list))) {
        empty_list[[length(empty_list) + 1]] <- id_sire[i]
        gen_vec[id_sire[i]] <- gen_vec[id_sire[i]] + 1
      }

      if (!is.na(id_dam[i]) && !(id_dam[i] %in% unlist(empty_list))) {
        empty_list[[length(empty_list) + 1]] <- id_dam[i]
        gen_vec[id_dam[i]] <- gen_vec[id_dam[i]] + 1
      }

    }

    if (length(empty_list) == 0) break

    id_ind <- unlist(empty_list)
  }

  return(ped[order(gen_vec, decreasing = TRUE),])
}


#' Evaluate whether the pedigree is ordered. NOTE: this has not been thoroughly vetted.
#'
#' @param ped a ped
#' @param id_col name or index of id column
#' @param sire_col name or index of sire column
#' @param dam_col name or index of dam column
#'
#' @return A logical value indicating whether or not the pedigree is ordered with offspring after parents.
#' @export
#'
check_order_ped <- function(ped,
                            id_col,
                            sire_col,
                            dam_col) {

  fid_index <- match(ped[[sire_col]], ped[[id_col]])
  mid_index <- match(ped[[dam_col]], ped[[id_col]])

  for (i in seq_len(nrow(ped))) {
    if (isTRUE(i < fid_index[i]) | isTRUE(i < mid_index[i])) {
      return(FALSE)
    }
  }
  return(TRUE)
}



#' Create indexed version of the pedigree
#'
#' @param ped pedigree (stored in a dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column3 --> dam id. Founders should have parents coded as 0.
#'
#' @return A dataframe containing the pedigree with ids transformed to the index values of each individual (and founder parents coded as 0s).
#' @export
#'
index_pedigree <- function(ped) {

  return(
    cbind(seq_along(ped[,1,drop=TRUE]),
          match(ped[,2,drop=TRUE], ped[,1,drop=TRUE], nomatch = 0),
          match(ped[,3,drop=TRUE], ped[,1,drop=TRUE], nomatch = 0))
  )
}




########################################
### MISCELLANEOUS PEDIGREE FUNCTIONS ###
########################################

#' Switch character IDs to numeric (for the ID, sire, and dam columns)
#'
#' @param ped pedigree
#' @param id_col name or index of id column
#' @param sire_col name or index of sire column
#' @param dam_col name or index of dam column
#' @param founder_val the value that represents a founder in the sire and dam cols
#'
#' @return a pedigree stored in a dataframe with updated individual ID values
#' @export
#'
character2numeric_id <- function(ped, id_col, sire_col, dam_col, founder_val = 0) {

  id_map <- data.frame(unique(ped[,id_col, drop = TRUE]),
                       seq_along(unique(ped[,id_col, drop = TRUE])))

  for (i in c(id_col, sire_col, dam_col)) {
    colnames(id_map) <- c(i, paste0(i, '_numeric'))
    ped <- dplyr::left_join(ped, id_map, by = i)
    ped[,paste0(i, '_numeric')][is.na(ped[,paste0(i, '_numeric')])] <- founder_val
  }

  return(ped)
}




#######################################
### EXTRACTION OF PEDIGREE ELEMENTS ###
#######################################
#' Extract ancestors for one or more individuals
#'
#' @param ped a pedigree object
#' @param indiv_vec the individuals that for which to extract ancestors
#' @param include_repeats if individuals are duplicate ancestors, should duplicates be removed?
#'
#' @return a vector of ancestors
#' @export
#'
get_ancestors <- function(ped, indiv_vec, include_repeats = TRUE) {

  #if (any(c('tbl_df', "tbl") %in% class(ped))) ped <- as.data.frame(ped)

  #extract all the ancestors for each individual in indiv_vec
  anc_vec <- unlist(ped[match(indiv_vec, ped[,1, drop = TRUE]), c(2, 3), drop = FALSE],
                    use.names = FALSE)
  #anc_vec_removeNA <- anc_vec[!is.na(anc_vec)] #remove NAs (these are founders)
  anc_vec_removeNA <- anc_vec[anc_vec != "0"]
  #if len of anc vec is 0 after removing NAs, no ancestors exist for the indivs
  if (length(anc_vec_removeNA) == 0) return(NULL)

  #output --> include_repeats != TRUE, return vector of ancestors without reps
  if (!isTRUE(include_repeats)) return(unique(anc_vec_removeNA))
  return(anc_vec_removeNA)
}


get_offspring <- function(ped, indiv_vec, include_repeats = TRUE) {

  offspring_vec <- ped[,1,drop = TRUE][ped[,2,drop = TRUE] %in% indiv_vec | ped[,3,drop = TRUE] %in% indiv_vec]
  if (length(offspring_vec) == 0) return(NULL)
  return(offspring_vec)
}


#' Extract the founders from a pedigree
#'
#' @param ped a pedigree object
#' @param founder_vals the founder values
#'
#' @return a list with the founders and partial founders (the individuals with only one known parent)
#' @export
#'
get_founders <- function(ped, founder_vals = '0') {
  founder_list <- list()
  founder_list[['founders']] <- ped[apply(ped[,c(2,3)], 1, function(x) all(x %in% '0')),1,drop=TRUE]
  founder_list[['partial_founders']] <- ped[apply(ped[,c(2,3)], 1, function(x) sum(x %in% '0') == 1),1,drop=TRUE]
  founder_list[sapply(founder_list, function(x) length(x) == 0)] <- 'none'

  return(founder_list)
}



#' Subset pedigree down to a focal set of individuals and their descendants
#'
#' @param ped pedigree (stored in dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column 3 --> dam. Founder parents should be coded as 0s.
#' @param indivs the individuals for whom you want to subset the pedigree
#'
#' @return a dataframe containing pedigree information for the focal individuals and their descendants based on the input pedigree
subset_ped_descendants <- function(ped, indivs) {

  indiv_list <- list(indivs)
  counter <- 1
  while (TRUE) {

    indiv_vec <- get_offspring(ped, indiv_list[[counter]])

    if (is.null(indiv_vec)) break

    indiv_list[[counter + 1]] <- indiv_vec[!indiv_vec %in% unlist(indiv_list)]

    counter <- counter + 1
  }

  return(unlist(indiv_list))

}


#' Subset pedigree down to a focal set of individuals and their ancestors
#'
#' @param ped pedigree (stored in dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column 3 --> dam. Founder parents should be coded as 0s.
#' @param indivs the individuals for whom you want to subset the pedigree
#'
#' @return a dataframe containing pedigree information for the focal individuals and their ancestors based on the input pedigree
#' @export
subset_ped_ancestors <- function(ped, indivs) {

  indiv_list <- list(indivs)
  counter <- 1
  while (TRUE) {

    indiv_vec <- get_ancestors(ped, indiv_list[[counter]],
                               include_repeats = FALSE)

    if (is.null(indiv_vec)) break

    indiv_list[[counter + 1]] <- indiv_vec[!indiv_vec %in% unlist(indiv_list)]

    counter <- counter + 1
  }

  return(unlist(indiv_list))

}






#################################
### CODE CURRENTLY NOT IN USE ###
#################################

# check_order_ped_old <- function(ped) {
#   return_val <- TRUE
#
#   for (i in seq_len(nrow(ped))) {
#
#     if (ped[i,2,drop=TRUE] != '0' && isTRUE(which(ped[i,2,drop=TRUE] == ped[,1,drop=TRUE]) > i)) {
#       return_val <- FALSE
#       break
#     }
#
#     if (ped[i,3,drop=TRUE] != '0' && isTRUE(which(ped[i,3,drop=TRUE] == ped[,1,drop=TRUE]) > i)) {
#       return_val <- FALSE
#       break
#     }
#   }
#
#   return(return_val)
# }


# examp1_pedsim <- sim_ped_repeat_attempt(attempts = 100,
#                                         founders = 5,
#                                         cycles = 10,
#                                         start_year = 1930,
#                                         age_mort_horizontal_shift = -1.5,
#                                         age_mort_slope = 0.5,
#                                         max_hatchling_survive_prob = 0.8,
#                                         prob_pair = 0.7,
#                                         add_immigrants = 'none',
#                                         immigrant_vec = NULL,
#                                         founder_age_lambda = 3,
#                                         offspring_lambda = 2,
#                                         immigrant_age = 5,
#                                         immigrant_prob = 5,
#                                         safeguard = 10000,
#                                         max_pop = 5000,
#                                         report_progress = FALSE)
#
# examp1_extract_ped <- extract_ped(examp1_pedsim,
#                                   remove_nonbreeding_founders = TRUE,
#                                   pedtools_format = TRUE)

