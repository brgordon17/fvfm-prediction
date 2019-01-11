#' Create mzdata.
#'
#' \code{create_mzdata()} pre-processes the LCMS data used for modelling in
#' gordon-C5.
#'
#' Initially, the function takes the raw output from xcms and removes unwanted
#' data (e.g. retention times, isotopes, peak counts etc.). Then, it creates
#' new categorical variables based on the sample information. Finally, it
#' replaces true non-detects with noise, removes poorly resolved mass features
#' and then replaces the small number of remaining missing values using random
#' forest imputaion. The list below details the logic behind the missing values
#' imputation:
#' \itemize{
#' \item \strong{TRUE NON-DETECTS:}
#'  Replace values missing in one class, but not others, with a random number
#'  between zero and the minimum of the matrix (i.e. noise). To be considered a
#'  true non-detect, a class should be missing at least 60 precent of its
#'  values. Achieved with,
#'  \code{metabolomics::MissingValues(group.cutoff = 0.6)}
#'  \item \strong{POORLY RESOLVED MASS FEATURES:}
#'  Remove mass features with more than 90 percent missing values. Achieved
#'  with, \code{metabolomics::MissingValues(column.cutoff = 0.9)}
#'  \item \strong{FALSE NON DETECTS:}
#'  Remaining missing values will be computed using
#'  \code{missForest::missForest()}. Achieved with,
#'  \code{metabolomics::MissingValues(complete.matrix = FALSE)}
#' }
#'
#' @author Benjamin R. Gordon
#'
#' @seealso
#' \code{\link[metabolomics]{MissingValues}}
#' \code{\link[doMC]{registerDoMC}}
#' \code{\link[missForest]{missForest}}
#' 
# Load libraries and data ------------------------------------------------------

library(tidyverse)
library(metabolomics)

mzdata  <-  readr::read_csv("./data-raw/mzdata-raw.csv", na = "0")
  
# remove isotopes --------------------------------------------------------------

mzdata <- mzdata[-grep("[M+1]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+2]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+3]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+4]", mzdata$isotopes, fixed = TRUE),]

# Clean up ---------------------------------------------------------------------

mzdata <- dplyr::select(mzdata, -isotopes, -adduct, -pcgroup)
mzdata <- mzdata[-2:-34]

mz_names <- round(mzdata[, 1], 5)
mz_names <- paste("mz", mz_names$mz, sep = "_")
mz_names <- make.names(mz_names, unique = TRUE)

mzdata <- tibble::as_tibble(t(mzdata), rownames = "sample_ids")
colnames(mzdata)[2:ncol(mzdata)] <- mz_names
mzdata <- mzdata[-1, ]
sample_ids <- mzdata$sample_ids
  
# Create phenodata -------------------------------------------------------------

phenodata <- tibble(sample_id = sample_ids,
                    class = c(rep("2011-T2-C", 6),
                              rep("2011-T2-T", 3),
                              rep("2011-T3-C", 6),
                              rep("2011-T3-T", 6),
                              rep("2011-T4-C", 6),
                              rep("2011-T4-T", 6),
                              rep("2013-PBQC", 18),
                              rep("2013-T0-C", 12),
                              rep("2013-T0-T", 12),
                              rep("2013-T1-C", 12),
                              rep("2013-T1-T", 9),
                              rep("2013-T2-C", 12),
                              rep("2013-T2-T", 12),
                              rep("2013-T3-C", 12),
                              rep("2013-T3-T", 12),
                              rep("2013-T4-C", 12),
                              rep("2013-T4-T", 12),
                              rep("2013-T5-C", 12),
                              rep("2013-T5-T", 12),
                              rep("2014-T1-C", 10),
                              rep("2014-T1-T", 10),
                              rep("2014-T2-C", 10),
                              rep("2014-T2-T", 10),
                              rep("2014-T3-C", 10),
                              rep("2014-T3-T", 8),
                              rep("2014-T4-C", 9),
                              rep("2014-T4-T", 10)
                              ),
                    experiment = c(rep("2011", 33),
                                   rep("2013", 159),
                                   rep("2014", 77)
                                   ),
                    day = c(rep("day 4", 6),
                            rep("day 4", 3),
                            rep("day 6", 6),
                            rep("day 6", 6),
                            rep("day 14", 6),
                            rep("day 14", 6),
                            rep("PBQC", 18),
                            rep("day 1", 12),
                            rep("day 1", 12),
                            rep("day 5", 12),
                            rep("day 5", 9),
                            rep("day 8", 12),
                            rep("day 8", 12),
                            rep("day 10", 12),
                            rep("day 10", 12),
                            rep("day 12", 12),
                            rep("day 12", 12),
                            rep("day 15", 12),
                            rep("day 15", 12),
                            rep("day 3", 10),
                            rep("day 3", 10),
                            rep("day 9", 10),
                            rep("day 9", 10),
                            rep("day 11", 10),
                            rep("day 11", 8),
                            rep("day 13", 9),
                            rep("day 13", 10)
                            ),
                    cont_treat = c(rep("C", 6),
                                   rep("T", 3),
                                   rep("C", 6),
                                   rep("T", 6),
                                   rep("C", 6),
                                   rep("T", 6),
                                   rep("PBQC", 18),
                                   rep("C", 12),
                                   rep("T", 12),
                                   rep("C", 12),
                                   rep("T", 9),
                                   rep("C", 12),
                                   rep("T", 12),
                                   rep("C", 12),
                                   rep("T", 12),
                                   rep("C", 12),
                                   rep("T", 12),
                                   rep("C", 12),
                                   rep("T", 12),
                                   rep("C", 10),
                                   rep("T", 10),
                                   rep("C", 10),
                                   rep("T", 10),
                                   rep("C", 10),
                                   rep("T", 8),
                                   rep("C", 9),
                                   rep("T", 10)
                                   ),
                    FvFm = c(rep(0.613, 6),
                             rep(0.578, 3),
                             rep(0.610, 6),
                             rep(0.523, 6),
                             rep(0.720, 6),
                             rep(0.502, 6),
                             rep(NA, 18),
                             rep(0.652, 12),
                             rep(0.646, 12),
                             rep(0.657, 12),
                             rep(0.641, 9),
                             rep(0.655, 12),
                             rep(0.575, 12),
                             rep(0.678, 12),
                             rep(0.573, 12),
                             rep(0.643, 12),
                             rep(0.304, 12),
                             rep(0.659, 12),
                             rep(0.000, 12),
                             rep(NA, 10),
                             rep(NA, 10),
                             rep(0.666, 10),
                             rep(0.511, 10),
                             rep(NA, 10),
                             rep(NA, 8),
                             rep(0.664, 9),
                             rep(0.600, 10)
                             ))

mzdata <- bind_cols(phenodata, mzdata[-1])


# Below works but is slow when performed on the whole dataset
mzdata[7:9] %>%
  mutate_if(is.character,as.numeric)









# Impute noise and remove unreliable mass features ---------------------------
mzdata_filt <- MissingValues(mzdata[c(-3:-6)],
                             column.cutoff = 0.8,
                             group.cutoff = 0.75,
                             complete.matrix = FALSE,
                             seed = 1978)
  mzdata <- data.frame(sample_ids,
                       class,
                       day,
                       tank,
                       rep,
                       batch,
                       mzdata_filt$output)
  percent_na <- round(sum(is.na(mzdata))/prod(dim(mzdata[-1:-6]))*100, 2)
  
  # Impute remaining missing values --------------------------------------------
  if(parallel) {
    doMC::registerDoMC()
    set.seed(seed)
    mzdata.imp <- missForest::missForest(mzdata[-1], parallelize = "variables")
  }
  
  else {
    set.seed(seed)
    mzdata.imp <- missForest::missForest(mzdata[-1], parallelize = "no")
  }
  
  mzdata <- tibble::as_tibble(cbind(sample_ids, mzdata.imp$ximp))
  mzdata <- dplyr::arrange(mzdata, class, day)
  
  # write data -----------------------------------------------------------------
  save(mzdata, file = "./data/mzdata.rda", compress = "bzip2")
  
  if(saverda) {
    devtools::use_data(mzdata)
  }
  
  if(savecsv) {
    readr::write_csv(mzdata, paste(c("./inst/extdata/", csvname, ".csv"),
                                   collapse = ""))
  }
  
  message("The data contained ", percent_na, "% NAs")########
  message("MissForest NRMSE: ", round(mzdata.imp$OOBerror, 4))
  mzdata


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
