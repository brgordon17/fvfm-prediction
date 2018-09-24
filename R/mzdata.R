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
#' @param parallel Logical indicating if missing values imputation should be run
#' in parallel. If \code{TRUE}, the default number of cores is equal to half the
#' available number of cores
#' @param seed An integer used for setting seeds of random number generation
#' @param savecsv Logical indicating if output should be saved as a \code{.csv}
#' file to the current working directory
#' @param saverda Logical indicating if a .rda file should be saved to /data
#' @param csvname The name of the output .csv file to be saved if TRUE
#'
#' @return Returns a dataframe of class tbl_df
#'
#' @note Using \code{parallel = TRUE} is not reproducible. Future versions of
#' this function may include support for reproducible RNG seeds when using
#' parallel processing. Although this function is exported,
#' \code{create_mzdata()} was not intended to be used outside of this package.
#'
#' @author Benjamin R. Gordon
#'
#' @seealso
#' \code{\link[metabolomics]{MissingValues}}
#' \code{\link[doMC]{registerDoMC}}
#' \code{\link[missForest]{missForest}}
#'
#' @export
#'
create_mzdata <- function(parallel = FALSE,
                          seed = 1978,
                          savecsv = FALSE,
                          saverda = TRUE,
                          csvname = "mzdata") {
  
  mzdata  <-  readr::read_csv("./data-raw/mzdata-raw.csv", na = "0")
  
  # remove isotopes ------------------------------------------------------------
  mzdata <- mzdata[-grep("[M+1]", mzdata$isotopes, fixed = TRUE),]
  mzdata <- mzdata[-grep("[M+2]", mzdata$isotopes, fixed = TRUE),]
  mzdata <- mzdata[-grep("[M+3]", mzdata$isotopes, fixed = TRUE),]
  mzdata <- mzdata[-grep("[M+4]", mzdata$isotopes, fixed = TRUE),]
  
  # Clean up -------------------------------------------------------------------
  mzdata <- dplyr::select(mzdata, -isotopes, -adduct, -pcgroup)
  mzdata <- mzdata[-2:-20]
  mz_names <- round(mzdata[, 1], 4)
  mz_names <- paste("mz", mz_names$mz, sep = "_")
  mz_names <- make.names(mz_names, unique = TRUE)
  mzdata <- tibble::as_tibble(t(mzdata), rownames = "sample_ids")
  colnames(mzdata)[2:ncol(mzdata)] <- mz_names
  mzdata <- mzdata[-1, ]
  sample_ids <- mzdata$sample_ids
  
  # Create categorical variables -----------------------------------------------
  class <- c(rep("PBQC", 18), 
             rep("cont", 12), 
             rep("treat", 12), 
             rep("cont", 12), 
             rep("treat", 9),
             rep("cont", 12), 
             rep("treat", 12),
             rep("cont", 12), 
             rep("treat", 12),
             rep("cont", 12), 
             rep("treat", 12),
             rep("cont", 12), 
             rep("treat", 12))
  class <- factor(class, levels = c("cont", "treat", 
                                    "PBQC"))
  
  day <- c(rep("PBQC", 18),
           rep("day0", 24),
           rep("day5", 21),
           rep("day8", 24),
           rep("day10", 24),
           rep("day12", 24),
           rep("day16", 24))
  day <- factor(day,
                levels = c("day0", "day5", "day8", "day10", "day12", "day16", 
                           "PBQC"))
  
  tank <- c(rep("PBQC", 18),
            c("02",  "02",  "02",  "13",  "13",  "13",  "40",  "40", "40", 
              "48",  "48", "48",  "20",  "20",  "20",  "25", "25",  "25", 
              "54", "54",  "54",  "55", "55", "55"),
            c("02",  "02",  "02",  "13",  "13",  "13",  "40",  "40", "40", 
              "48",  "48", "48",  "25", "25",  "25", "54", "54",  "54", 
              "55", "55", "55"),
            rep(c("02",  "02",  "02",  "13",  "13",  "13",  "40",  "40", "40", 
                  "48",  "48", "48",  "20",  "20",  "20",  "25", "25",  "25", 
                  "54", "54",  "54",  "55", "55", "55"), 4))
  tank <- factor(tank, 
                 levels = c("02", "13", "40", "48", "20", "25", "54", "55", 
                            "PBQC"))
  
  rep <- c(c(1, 1, 10, 11, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9), 
           rep(c(1, 2, 3), 47))
  rep <- factor(rep, levels = (c(1:11)))
  
  batch <- c(1, 2, 2, 2, rep(c(1, 2), 6), 2, 2,
             rep(2, 24), 
             rep(1, 21),
             rep(2, 24),
             rep(2, 24),
             rep(1, 24),
             rep(2, 24))

  class_day <- interaction(class,
                           day,
                           drop = TRUE,
                           sep = ".")
  
  # Impute noise and remove unreliable mass features ---------------------------
  mzdata <- tibble::as_tibble(cbind(class_day, mzdata[-1]))
  mzdata_filt <- metabolomics::MissingValues(mzdata,
                                             column.cutoff = 0.8,
                                             group.cutoff = 0.75,
                                             complete.matrix = FALSE,
                                             seed = seed)
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
}

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
