# Code to create litmz.rda
# Author: Benjamin R. Gordon
#Date: 2019-03-11

# Load literature variables and clean up ---------------------------------------
litmz <- readr::read_csv("./data-raw/metabs_literature.csv", col_names = TRUE)
litmz <- rename(litmz,
                reported_mz = `Reported m/z`,
                reported_ion = `Reported ion`,
                monoiso_mass = `Monoisotopic Mass (Da)`,
                molec_formula = `Molecular Formula`,
                reported_name = `Reported compound name`,
                taxon = `Genus`,
                endnote_ref = `Endnote Reference`,
                ref = `ref`)
save(litmz, file = "./data/litmz.rda", compress = "bzip2")

## Document data ---------------------------------------------------------------

#' litmz
#'
#' A dataset of published and/or known mz ions identified in either coral,
#' Symbiodinium or other closely related species.
#'
#' @format A tibble with 230 rows and 7 variables:
#' \describe{
#'   \item{reported_mz}{the m/z reported}
#'   \item{reported ion}{the ion reported}
#'   \item{monoiso_mass}{the monoisotopic mass}
#'   \item{molec_formula}{the molecular formula}
#'   \item{reported_name}{the reported name}
#'   \item{taxon}{the taxonomic information}
#'   \item{ref}{the reference where the information was reported. This is
#'   currently in custom endnote format and will be updated in the future.}
#' }
#' @source Benjamin R. Gordon
"litmz"
