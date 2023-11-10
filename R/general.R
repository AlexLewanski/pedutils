harmonic_mean <- function(x) length(x)/sum(1/x)

get_ancestors <- function(ped, indiv_vec, include_repeats = TRUE) {

  #if (any(c('tbl_df', "tbl") %in% class(ped))) ped <- as.data.frame(ped)

  #extract all the ancestors for each individual in indiv_vec
  anc_vec <- unlist(ped[match(indiv_vec, ped[,1]), c(2, 3)], use.names = FALSE)
  #anc_vec_removeNA <- anc_vec[!is.na(anc_vec)] #remove NAs (these are founders)
  anc_vec_removeNA <- anc_vec[anc_vec != "0"]
  #if len of anc vec is 0 after removing NAs, no ancestors exist for the indivs
  if (length(anc_vec_removeNA) == 0) return(NULL)

  #output --> include_repeats != TRUE, return vector of ancestors without reps
  if (!isTRUE(include_repeats)) return(unique(anc_vec_removeNA))
  return(anc_vec_removeNA)
}
