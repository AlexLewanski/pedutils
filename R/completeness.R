calc_max_gen <- function(ped, indiv = 'all') {

  if (identical('all', indiv)) indiv <- ped[, 1, drop = TRUE]
  max_gen_vec <- stats::setNames(vector(mode = 'integer', length = length(indiv)),
                                 nm = indiv)

  for (i in indiv) {
    gen <- 0
    indiv_vec <- i
    while(TRUE) {
      indiv_vec <- get_ancestors(ped = ped,
                                 indiv_vec = indiv_vec,
                                 include_repeats = FALSE)
      if (is.null(indiv_vec)) break
      gen <- gen + 1
    }

    max_gen_vec[i] <- gen
  }

  return(max_gen_vec)
}

calc_total_gen <- function(complete_info) {
  return(
    stats::aggregate(generation ~ id,
                     data = complete_info,
                     FUN = max,
                     simplify = TRUE)[, 2, drop = TRUE]
  )
}


calc_equiv_gen <- function(complete_info) {
  return(
    stats::aggregate(completeness ~ id,
                     data = complete_info,
                     FUN = function(x) sum(x) - 1,
                     simplify = TRUE)[, 2, drop = TRUE]
  )
}


calc_complete_gen <- function(complete_info) {
  #dplyr::filter(complete_info, abs(1 - completeness) <= .Machine$double.eps),
  return(
    stats::aggregate(generation ~ id,
                     complete_info[abs(1 - complete_info$completeness) <= .Machine$double.eps,],
                     FUN = max)[, 2, drop = TRUE]
  )
}

calc_mean_compl <- function(complete_info, gen = 1) {
  #dplyr::filter(complete_info,generation <= (gen - 1)),
  return(
    stats::aggregate(completeness ~ id,
                     complete_info[complete_info$generation <= (gen - 1),],
                     FUN = function(x) (sum(x))/gen)[, 2, drop = TRUE]
  )
}

calc_compl_index <- function(ped, complete_info, gen = 1) {

  if (!identical(unique(complete_info[, 1, drop = TRUE]), ped[, 1, drop = TRUE])) {
    stop('the order of individuals is different')
  }

  #mean completeness of each individual
  ped$mean_compl <- calc_mean_compl(complete_info, gen = gen)

  #add extra columns that contain the parents' mean completeness values for each indiv
  #col 1 --> id; col 2 --> sire id; col3 --> dam id
  ped$sire_meancompl <- ped$mean_compl[match(ped[,2,drop = TRUE], ped[,1,drop = TRUE])]
  ped$dam_meancompl <- ped$mean_compl[match(ped[,3,drop = TRUE], ped[,1,drop = TRUE])]

  #add 0s for the mean parental completeness vals for indivs w/out known parents
  ped$sire_meancompl[is.na(ped$sire_meancompl)] <- 0
  ped$dam_meancompl[is.na(ped$dam_meancompl)] <- 0

  #calculate completeness index by taking the harmonic mean of the sire and dam
  #mean completeness values
  return(apply(ped[,c('sire_meancompl', 'dam_meancompl')],
               MARGIN = 1,
               harmonic_mean)
  )
}


#' Calculate completeness for individuals in the pedigree
#'
#' @param ped A pedigree dataframe in the format output by the process_ped function.
#' @param id The id of the individual(s) whose completeness you want calculated. The option 'all' can be used if you want the completeness of all individuals to be calculated.
#' @param max_gen For each individual, the maximum number of generations for which to calculate completeness. The option 'all' will result in completeness being calculated for all generations present in the generation for each individual.
#' @return A dataframe containing completeness information for individuals in the pedigree.
#' @import dplyr
#' @export
#'
#' @examples
calc_completeness <- function(ped, id = 'all', max_gen = 'all') {

  ### ARGUMENT CHECKS AND PROCESSING ###
  if (identical(max_gen, 'all')) max_gen <- Inf #set max_gen to inf if no max_gen is given
  if (identical('all', id)) id <- ped[,1,drop=TRUE] #if id = all, calc completeness for all indivs
  #create empty list to store each individuals completeness df
  df_list <- list()

  ### COMPLETENESS CALCULATIONS ###
  for (i in id) {
    gen <- 1 #set gen counter to 1
    anc_vec <- i #assign focal indiv (i) the ancestor vector

    #completeness of gen 0 is always 100%, so set the 0th gen entry to 1
    compl_list <- list(1)
    #compl_list <- list()

    #traverse each gen while gen is <= max_gen
    while (gen <= max_gen) {

      #retrieve parents of individuals in anc_vec
      anc_vec <- get_ancestors(ped = ped,
                               indiv_vec = anc_vec,
                               include_repeats = TRUE)

      #if all indivs have no parents (all are founders), leave loop
      if (is.null(anc_vec)) break


      #completeness = number of ancestors in ped/number of possible ancestors
      compl_list[[gen + 1]] <- length(anc_vec)/(2^gen)
      gen <- gen + 1 #add one to the gen counter
    }

    #add info to dataframe
    #subtract 1 from gen because an extra 1 is added to it in the last iteration

    df_list[[i]] <- data.frame(id = i,
                               generation = 0:(gen - 1),
                               completeness = unlist(compl_list))
  }

  return(dplyr::bind_rows(df_list))
}



#' Calculate completeness and depth measures of a pedigree
#'
#' @param ped A pedigree dataframe in the format output by the process_ped function.
#' @param gen The number of generations to consider for the compl_index calculation
#' @param summary The measures used to summarize the pedigree
#'
#' @return A dataframe with the completeness and/or depth measures of the pedigree
#' @export
#'
#' @examples
ped_summary_stats <- function(ped,
                              gen = 1,
                              summary = c('total_gen', 'equiv_gen', 'complete_gen', 'mean_compl', 'compl_index')) {

  summary <- unique(match.arg(summary,
                              several.ok = TRUE))

  comp_df <- calc_completeness(ped) #calculate completeness for each indiv

  summary_df <- ped[, 1, drop = FALSE] #start output df using the id col of ped

  #iterate through each summary measure and calculate it
  for (i in summary) {
    summary_df[,i] <- switch(i,
                             total_gen = {calc_total_gen(comp_df)},
                             equiv_gen = {calc_equiv_gen(comp_df)},
                             complete_gen = {calc_complete_gen(comp_df)},
                             mean_compl = {calc_mean_compl(comp_df, gen = gen)},
                             compl_index = {calc_compl_index(ped, comp_df, gen = gen)})
  }

  return(summary_df)
}

