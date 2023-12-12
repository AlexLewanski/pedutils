#############################
### ALL-PURPOSE FUNCTIONS ###
#############################
harmonic_mean <- function(x) length(x)/sum(1/x)

#######################################
### PEDIGREE MANIPULATON AND CHECKS ###
#######################################
reorder_ped <- function(ped) {

  #Zhang et al. 2009
  #An algorithm to sort complex pedigrees chronologically without birthdates

  id_ind <- seq_along(ped$id)
  id_sire <- match(ped$fid, ped$id)
  id_dam <- match(ped$mid, ped$id)

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


check_order_ped <- function(ped) {
  fid_index <- match(ped[,2,drop = TRUE], ped[,1,drop = TRUE])
  mid_index <- match(ped[,3,drop = TRUE], ped[,1,drop = TRUE])

  for (i in seq_len(nrow(ped))) {
    if (isTRUE(i < fid_index[i]) | isTRUE(i < mid_index[i])) {
      return(FALSE)
    }
  }
  return(TRUE)
}



#######################################
### EXTRACTION OF PEDIGREE ELEMENTS ###
#######################################
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

