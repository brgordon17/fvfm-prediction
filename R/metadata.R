# Experiment metadata ------------------------------------------------------------------
day <- c(0:16)
date <- seq(as.Date("2013-05-06"), as.Date("2013-05-22"), by = "days")
avetemp <- c(seq(from = 25, to = 32.3, by = 0.7), 33, 33, 33, 34, 34, 34)
collect <- as.logical(c(1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1))
pam <- as.logical(c(1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
cell_count <- as.logical(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1))
pigment <- as.logical(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1))
exp_metadata <- 
  tibble::tibble(day, date, avetemp, collect, pam, cell_count, pigment)

# save
save(exp_metadata, file = "./data/exp_metadata.rda", compress = "bzip2")

# Export csv
readr::write_csv(exp_metadata, "./inst/extdata/exp_metadata.csv")

#' experiment metadata
#'
#' The metadata for the bleaching experiment performed in May 2013
#'
#' @format A tibble with 17 rows and 7 variables:
#' \describe{
#'   \item{day}{sampling timepoints}
#'   \item{date}{the date}
#'   \item{avetemp}{average temperature}
#'   \item{collect}{were samples collected}
#'   \item{pam}{was PAM fluorometry performed}
#'   \item{cell_count}{were cell counts conducted}
#'   \item{pigment}{was pigment analysis conducted}
#' }
#' @source Benjamin R. Gordon
"exp_metadata"