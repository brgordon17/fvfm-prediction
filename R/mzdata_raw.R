# Feature detection and annotation using xcms and CAMERA
# Benjamin Gordon 8-JAN-2019
# See ./dev/mzdata-raw-dev.R for the developmental version of this code.

library(xcms)
library(magrittr)
library(CAMERA)

# Register parallel processing -------------------------------------------------
register(bpstart(MulticoreParam(workers = 4)))

# Data import ------------------------------------------------------------------

# path to files
mzfiles <- list.files(path = "~/Big Data/LCMS Data/chapter5_mzxml",
                      full.names = TRUE,
                      recursive = TRUE)

# create metadata or phenodata
phenod <- data.frame(sample_name = sub(basename(mzfiles), pattern = ".mzXML",
                                       replacement = "", fixed = TRUE),
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
                     stringsAsFactors = FALSE)

# load raw data and save
raw_data <- readMSData(files = mzfiles, pdata = new("NAnnotatedDataFrame", 
                                                    phenod),
                       mode = "onDisk")

saveRDS(raw_data, "./data-raw/temp-saves/OnDiskMSnExp.rds")

# data inspection --------------------------------------------------------------

# ID problematic LCMS runs using boxplots of the TIC
tics <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tics, ylab = "intensity", main = "TIC")

# Peak detection ---------------------------------------------------------------

# define the parameters 
cwp <- CentWaveParam(ppm = 30, 
                     peakwidth = c(10, 60),
                     snthresh = 10,
                     prefilter = c(3, 100),
                     mzCenterFun = "wMean",
                     integrate = 1L,
                     mzdiff = -0.001,
                     fitgauss = FALSE,
                     noise = 1000)

# Do peak detection and save
mzdata <- findChromPeaks(raw_data, param = cwp)
saveRDS(mzdata, "./data-raw/temp-saves/XCMSnExp.rds")

# Alignment --------------------------------------------------------------------

mzdata <- adjustRtime(mzdata, param = ObiwarpParam(binSize = 0.5))

# plot retention time adjustment
plotAdjustedRtime(mzdata)

# Correspondence or Grouping ---------------------------------------------------

pdp <- PeakDensityParam(mzdata$class,
                        bw = 5,
                        minFraction = 0.5,
                        binSize = 0.025)
mzdata <- groupChromPeaks(mzdata, param = pdp)

# fill missing peaks -----------------------------------------------------------
mzdata <- fillChromPeaks(mzdata)
saveRDS(mzdata, "./data-raw/temp-saves/XCMSnExp.rds")
processHistory(mzdata)

# Isotope and adduct annotation ------------------------------------------------

# convert to xcmsset object
mzdata <- as(mzdata, "xcmsSet")

# Create report with with isotope and adduct info
mzdata <- annotate(mzdata, 
                   perfwhm = 0.7, 
                   cor_eic_th = 0.75, 
                   ppm = 10, 
                   polarity = "positive")

# save
saveRDS(mzdata, "./data-raw/temp-saves/xsAnnotate.rds")

# Tidy and save ----------------------------------------------------------------

# create tibble
mzdata_raw <- tibble::as_tibble(CAMERA::getPeaklist(mzdata))

# reattach sample names
names(mzdata_raw)[35:303] <- phenod$sample_name

# save mzdata-raw as csv
readr::write_csv(mzdata_raw, "./data-raw/mzdata-raw.csv")

# save phenodata as csv
readr::write_csv(phenod, "./data-raw/phenodata-raw.csv")

# remove temporary saves
file.remove(list.files("./data-raw/temp-saves", full.names = TRUE))
