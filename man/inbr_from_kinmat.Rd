% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relate_quantities.R
\name{inbr_from_kinmat}
\alias{inbr_from_kinmat}
\title{Calculate inbreeding from a kinship/coancestry matrix}
\usage{
inbr_from_kinmat(ped, kin_mat = NULL, id_col = 1, sire_col = 2, dam_col = 3)
}
\arguments{
\item{ped}{the pedigree from which to calculate inbreeding}

\item{kin_mat}{the kinship matrix from to which to calculate inbreeding. If this is not provided, the matrix is calculated internally.}

\item{id_col}{the name or index of the id column}

\item{sire_col}{the name or index of the sire column}

\item{dam_col}{the name or index of the dam column}
}
\value{
a dataframe with inbreeding values for each individual
}
\description{
Calculate inbreeding from a kinship/coancestry matrix
}
