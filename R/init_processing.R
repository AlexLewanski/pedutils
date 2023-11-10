ped_add_missing_indivs <- function(ped, id, fid, mid, sex, founder_par_val = '0') {

  ped_col_names <- colnames(ped)
  parent_vec <- unique(c(ped[, fid, drop = TRUE],
                         ped[, mid, drop = TRUE]))

  parent_vec <- parent_vec[!parent_vec %in% founder_par_val]

  absent_parent <- parent_vec[!parent_vec %in% ped[, id, drop = TRUE]]

  if (length(absent_parent) == 0) return(ped)

  sex_vec <- sapply(absent_parent, function(id, df) {c(1:2)[apply(df, 2, function(x) id %in% x)]},
                    df = ped[,c('fid', 'mid')])

  new_id_df <- stats::setNames(data.frame(absent_parent, founder_par_val, founder_par_val, sex_vec),
                               nm = c(id, fid, mid, sex))

  new_id_df[ped_col_names[!ped_col_names %in% colnames(new_id_df)]] <- NA

  new_id_df_order <- new_id_df[, match(ped_col_names, colnames(new_id_df))]

  return(rbind(ped, new_id_df_order))
}


col_process_internal <- function(ped,
                                 id_col,
                                 sire_col,
                                 dam_col,
                                 sex_col,
                                 keep_extra_cols) {

  colname_vec <- stats::setNames(c('id', 'fid', 'mid', 'sex'),
                                 nm = c(id_col, sire_col, dam_col, sex_col))

  if (any(!names(colname_vec) %in% colnames(ped))) {
    stop('The following column(s) are missing from the pedigree: ',
         paste(names(colname_vec)[!names(colname_vec) %in% colnames(ped)], collapse = ', '))
  }

  #rename columns
  colnames(ped)[match(names(colname_vec), colnames(ped))] <- colname_vec

  if (keep_extra_cols) {
    colname_vec <- c(colname_vec, colnames(ped)[!colnames(ped) %in% colname_vec])
  }

  return(ped[,colname_vec,drop=FALSE])
}

sex_eval_internal <- function(processed_ped) {

  ### dam/sire overlap ###
  sex_overlap <- intersect(processed_ped$fid, processed_ped$mid)

  if (length(sex_overlap[sex_overlap != '0']) > 0) {
    stop('The following individuals are recorded as both a sire and dam: ',
         paste(unique(sex_overlap[sex_overlap != '0']), collapse = ', '))
  }

  ### identifying conflicts between parent status and specified sex ###
  males_unknown <- processed_ped$id[processed_ped$sex %in% c(1, 0)]
  females_unknown <- processed_ped$id[processed_ped$sex %in% c(2, 0)]

  male_unk_dams <- males_unknown[males_unknown %in% processed_ped$mid]
  female_unk_sires <- females_unknown[females_unknown %in% processed_ped$fid]
  combined_sex_mismatch <- c(male_unk_dams, female_unk_sires)

  if (length(combined_sex_mismatch) > 0) {

    sex_correct_vec <- c(rep(2, length(male_unk_dams)),
                         rep(1, length(female_unk_sires)))

    processed_ped$sex[match(combined_sex_mismatch, processed_ped$id)] <- sex_correct_vec
  }

  return(processed_ped)
}


#' Process a pedigree
#'
#' @param ped The unprocessed dataframe. It needs to contain the following columns: id_col (the ID of each individual), sire_col (the ID of the sire), dam_col (the ID of the dam), and sex_col (the sex of the individual).
#' @param id_col The name of the id column.
#' @param sire_col The name of the sire (father) column.
#' @param dam_col The name of the dam (mother) column.
#' @param sex_col The name of the sex column.
#' @param founder_val One or more (stored as a vector) values that represent the value for founder parents.
#' @param sex_vals A list containing the values representing each sex in sex_col. The list should be in the following format: list(male = c('male_val1', male_val2'), female = 'female_val', unknown = c('u1', 'u2', NA)).
#' @param keep_extra_cols Logical value (TRUE or FALSE) indicating whether to keep extra columns in the pedigree beyond the id, sire, dam, and sex columns.
#' @param disable_sex_check Logical value (TRUE or FALSE) indicating whether sex information should be evaluated.
#'
#' @return A processed pedigree stored as a dataframe.
#' @export
#'
#' @examples
process_ped <- function(ped,
                        id_col,
                        sire_col,
                        dam_col,
                        sex_col,
                        founder_val,
                        sex_vals,
                        keep_extra_cols = TRUE,
                        disable_sex_check = FALSE) {

  #===============
  #== OVERVIEW ===
  #===============

  #make the ped object a dataframe (if it isn't already)
  #rename and reorder the id, parent, and sex columns
  #recode the founder parent values to '0'
  #add individuals who are parents but aren't included in the id column
  #process and evaluate sex info:
  #   - recode sex categories to 1 (male), 2 (female), and 0 (unknown)
  #   - make sure single individuals aren't recorded as both sire and dams
  #   - if an individual has as unknown sex but is included as a sire or dam,
  #     update the sex info based on the sex of the parent type that it was
  #     recorded as


  #================================================
  #== CHECKING AND INITIAL PROCESSING OF INPUTS ===
  #================================================
  if (!identical(class(ped), 'data.frame')) ped <- as.data.frame(ped)
  stopifnot(identical(class(ped), 'data.frame'))

  for (i in list(id_col, sire_col, dam_col, sex_col)) {
    stopifnot(length(i) == 1 && class(i) == 'character')
  }

  stopifnot(is.list(sex_vals) && identical(sort(names(sex_vals)), c("female", "male", "unknown") ))

  #this should get all the argument inputs in a list but idk if this a good idea
  #c(as.list(environment()), list(...))
  logic_list <- list(keep_extra_cols, disable_sex_check)
  names(logic_list) <- c('keep_extra_cols', 'disable_sex_check')

  for (i in names(logic_list)) {
    if (! (isFALSE(logic_list[[i]]) | isTRUE(logic_list[[i]])) )
      stop('The following argument is not a single TRUE or FALSE value: ', i)
  }

  if (any(duplicated(ped[,id_col,drop = TRUE]))) {
    stop('The following individuals have >1 entry in the id column: ',
         paste(
           unique(ped[,id_col,drop = TRUE][duplicated(ped[,id_col,drop = TRUE])]),
           collapse = ', ')
    )
  }


  #========================
  #== COLUMN PROCESSING ===
  #========================

  #rename and reorder columns
  ped_processed <- col_process_internal(ped = ped,
                                        id_col = id_col,
                                        sire_col = sire_col,
                                        dam_col = dam_col,
                                        sex_col = sex_col,
                                        keep_extra_cols = keep_extra_cols)

  #if they aren't already, coerce id, fid, and mid to characters
  for (COL in c('id', 'fid', 'mid')) {
    if (!is.character(ped_processed[,COL,drop=TRUE])) {
      ped_processed[,COL] <- as.character(ped_processed[,COL,drop=TRUE])
    }
  }


  #===============================
  #== PEDIGREE INFO PROCESSING ===
  #===============================

  for (i in c('fid', 'mid')) {
    ped_processed[,i][ped_processed[,i,drop = TRUE] %in% founder_val] <- "0"
  }

  if (!all(ped_processed$sex %in% unlist(sex_vals)))
    stop('Sex is encoded with different(or additional) values than those in the sex_vals list')
  ped_processed$sex_old <- ped_processed$sex
  ped_processed$sex <- NA

  for (i in names(sex_vals) ) {
    ped_processed$sex[ped_processed$sex_old %in% sex_vals[[i]]] <- switch(i,
                                                                          male = 1,
                                                                          female = 2,
                                                                          unknown = 0)
  }

  ped_processed <- ped_add_missing_indivs(ped = ped_processed,
                                          id = 'id',
                                          fid = 'fid',
                                          mid = 'mid',
                                          sex = 'sex',
                                          founder_par_val = '0')

  if (!disable_sex_check) {
    ped_processed <- sex_eval_internal(processed_ped = ped_processed)
  }

  return(ped_processed)
}
