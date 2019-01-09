# Feature detection and annotation using xcms and CAMERA
# Benjamin Gordon 8-JAN-2019
# See ./dev/mzdata-raw-dev.R for the developmental version of this code.

library(xcms)
library(magrittr)
library(CAMERA)

# Create phenodata (sample group info only) ------------------------------------


# Data import ------------------------------------------------------------------

# path to files
mzfiles <- list.files(path = "./dev/test-data",
                      full.names = TRUE,
                      recursive = TRUE)




















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
