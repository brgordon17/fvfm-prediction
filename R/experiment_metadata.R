

# Create vectors for day, average set temp, collection (y/n)
day <- c(0:16)
date <- seq(as.Date("2013-05-06"), as.Date("2013-05-22"), by="days")
avetemp <- c(seq(from = 25, to = 32.3, by = 0.7), 33, 33, 33, 34, 34, 34)
collect <- c(1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1)
expmd <- data.frame(day, date, avetemp, collect)

# save
saveRDS(expcond, "./data/processed/expcond.rds")

# Export csv
readr::write_csv(expmd)

## Data documentation ----------------------------------------------------------

#' Preprocessed LCMS data
#'
#' A dataset of positive ESI-LCMS mz variables for coral samples subjected to
#' experimental ocean warming and acidification at Heron Island in Feb 2011.
#'
#' @format A tibble with 101 rows and 1647 variables:
#' \describe{
#'   \item{sample_ids}{ID's for each sample}
#'   \item{class}{Treatment group. One of, control, eT (eleveated temperature),
#'   eCO2 (elevated CO2), eCO2eT (elevated CO2 and temperature), PBQC}
#'   \item{day}{Days of exposure to treatment}
#'   \item{tank}{Tank membership. Either L (left) or R (right)}
#'   \item{rep}{Sample replicate number}
#'   \item{class_day}{interaction variable between class and day variables}
#' }
#' @source Benjamin R. Gordon
"mzdata"