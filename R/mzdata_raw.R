#' Feature detection using xcms.
#'
#' Method used by Gordon et. al. to extract LCMS features from
#' \code{.mzXML} files using \code{xcms}.
#'
#' This function was used to extract features from \code{.mzXML} data. The
#' function creates an \code{xcmsSet} object, groups peaks, corrects for
#' retention time drift, re-groups the peaks, fills missing peaks, annotates
#' adducts and isotopes, before saving the output to a \code{.csv} file in
#' \code{/data-raw}.
#'
#' @param saveoutput Logical indicating if output should be saved
#' @param outputname The name of the output file to be saved if \code{TRUE}
#' @param ... Other arguments passed on to individual methods
#'
#' @return Returns a dataframe of class \code{tbl_df}
#'
#' @note \code{mzdata_raw()} was not intended to be used outside of this
#' package. To keep the size of this package small, only a small example subset
#' of our LCMS \code{.mzXML} data is made available in
#' \code{./data-raw/example_mzxml_data/}. The full data set is available from
#' the package author.
#'
#' @author Benjamin R. Gordon
#'
#' @seealso
#' \code{\link[xcms]{xcmsSet}}
#' \code{\link[CAMERA]{annotate}}
#'
#' @export
#'
mzdata_raw <- function(saveoutput = FALSE,
                       outputname = "mzdata-raw",
                       ...) {
  
  files <- list.files(path = "./data-raw/mzxml",
                      full.names = TRUE,
                      recursive = TRUE)
  
  xset <- xcms::xcmsSet(files,
                        method = "centWave",
                        ppm = 30,
                        snthr = 10,
                        peakwidth = c(20,60),
                        mzdiff = 0.01,
                        integrate = 2,
                        prefilter = c(3,1100)
  )
  
  # Group peaks
  xsetgr1 <- xcms::group(xset,
                         bw = 5,
                         minfrac = 0.5,
                         minsamp = 1,
                         mzwid = 0.05,
                         max = 100
  )
  
  # Correct RT
  xsetcor <- xcms::retcor(xsetgr1,
                          method = "obiwarp",
                          plottype = "none",
                          profStep = 0.5
  )
  
  # Regroup peaks
  xsetgr2 <- xcms::group(xsetcor,
                         bw = 5,
                         minfrac = 0.5,
                         minsamp = 1,
                         mzwid= 0.05,
                         max = 100
  )
  
  # Fill missing peaks
  xsetmiss <- xcms::fillPeaks(xsetgr2)
  
  # Identify and annotate adducts
  xsetadd <- CAMERA::annotate(xsetmiss,
                              perfwhm = 0.7,
                              cor_eic_th = 0.75,
                              ppm = 10,
                              polarity = "positive"
  )
  
  # Get features
  features <- tibble::as_tibble(CAMERA::getPeaklist(xsetadd))
  
  # Save output
  if (saveoutput) {
    readr::write_csv(features, paste(c("./data-raw/",outputname, ".csv"),
                                     collapse = ""))
  }
  
  features
}
