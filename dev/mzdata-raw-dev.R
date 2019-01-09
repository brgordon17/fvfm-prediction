# This is a developmental script for creating mzdata using the new xcms
# functions. 16-12-18

library(xcms)
library(magrittr)
library(CAMERA)

# Data import ------------------------------------------------------------------

# path to files
mzfiles <- list.files(path = "./dev/test-data",
                     full.names = TRUE,
                     recursive = TRUE)

# create metadata or phenodata
phenod <- data.frame(sample_name = sub(basename(mzfiles), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 class = c(rep("eCO2eT", 4), rep("recovery", 4)),
                 stringsAsFactors = FALSE)

# load raw data
raw_data <- readMSData(files = mzfiles, pdata = new("NAnnotatedDataFrame", 
                                                    phenod),
                       mode = "onDisk")

# data inspection --------------------------------------------------------------

# plot BPCs
bpcs <- chromatogram(raw_data, aggregationFun = "max")
class_colours <- c("grey80", "red2")
names(class_colours) <- c("eCO2eT", "recovery")
plot(bpcs, col = class_colours[raw_data$class])

# ID problematic LCMS runs using boxplots of the TIC
tics <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tics, col = class_colours[raw_data$class],
        ylab = "intensity", main = "TIC")

# peak inspection --------------------------------------------------------------

# NOTE: THis step is not useful without known compounds or internal standards.
# check peak width by selecting a peak that is common to all samples as IDd in 
# the BPC and mMass.app
# define rt range
# rtr <- c(20, 100)
# mzr <- c(136.05, 136.07)
# 
# # extract chromatogram
# chr_raw <- chromatogram(raw_data, rt = rtr, mz = mzr)
# plot(chr_raw, col = class_colours[raw_data$class])
# 
# # check ppm of above peak
# raw_data %>%
#   filterRt(rt = rtr) %>%
#   filterMz(mz = mzr) %>%
#   plot(type = "XIC")

# Peak detection ---------------------------------------------------------------

# define the parameters used for peak detection (centwave here)
cwp <- CentWaveParam(ppm = 30, 
                     peakwidth = c(10, 60),
                     snthresh = 10,
                     prefilter = c(3, 100),
                     mzCenterFun = "wMean",
                     integrate = 1L,
                     mzdiff = -0.001,
                     fitgauss = FALSE,
                     noise = 1000)

# Do peak detection
mzdata <- findChromPeaks(raw_data, param = cwp)

# # view peaks given an mz and rt range
# chromPeaks(mzdata, mz = c(136, 137), rt = c(0, 200))
# 
# # view peaks given an mz and a ppm error
# chromPeaks(mzdata, mz = 136.1076, ppm = 10 )

# Alignment --------------------------------------------------------------------
# THis method uses the obiwarp method and is much better than the peak groups
# method that was trialled below.
mzdata <- adjustRtime(mzdata, param = ObiwarpParam(binSize = 0.5))

# compare adjusted BPCs with unadjusted
# original
bpcs <- chromatogram(raw_data, aggregationFun = "max")
class_colours <- c("grey80", "red2")
names(class_colours) <- c("eCO2eT", "recovery")
plot(bpcs, col = class_colours[raw_data$class])

# adjusted
bpcs_adj <- chromatogram(mzdata, aggregationFun = "max")
plot(bpcs_adj, col = class_colours[raw_data$class])

# plot retention time adjustment
plotAdjustedRtime(mzdata, col = class_colours[mzdata$class])

# NOTE: THE METHOD BELOW DOES NOT WORK ON MY SAMPLES.
# begin alignment by identifying peak groups (i.e. features or "hooks" that are
# present in most samples) that span the majority of the rt range. This is done
# by running an initial correspondence run to ID groups. Minfrac should be set
# high for this step.

# define grouping parameters and perform grouping
# pdp <- PeakDensityParam(sampleGroups = mzdata$class,
#                         minFraction = 0.8)
# mzdata <- groupChromPeaks(mzdata, param = pdp)
# 
# # define peak groups that will be used for alignment
# pgp <- PeakGroupsParam(minFraction = 0.85)
# 
# # view peak groups that will be used for alignment. Returns a matrix of
# # retention times for each feature (rows) and sample (cols).
# adjustRtimePeakGroups(mzdata, param = pgp)
# 
# # align
# mzdata <- adjustRtime(mzdata, param = pgp)
# 
# # view plots
# class_colours <- c("blue", "red")
# names(class_colours) <- c("eCO2eT", "recovery")
# plotAdjustedRtime(mzdata, 
#                   col = class_colours[mzdata$class],
#                   peakGroupsCol = "grey",
#                   peakGroupsPch = 1)
#

# Correspondence or Grouping ---------------------------------------------------
# define mz range
# mzr <- c(136.0, 138.1)
# 
# # Extract and plot the chromatograms
# chr_mzr <- chromatogram(mzdata, mz = mzr, rt = c(20, 120))
# cols <- class_colours[chr_mzr$class]
# plot(chr_mzr, col = cols)
# 
# # Highlight the detected peaks in that region.
# highlightChromPeaks(mzdata, mz = mzr, col = cols, type = "point", pch = 16)

# Define the parameters for the peak density method
# pdp <- PeakDensityParam(sampleGroups = mzdata$class,
#                         minFraction = 0.5, bw = 30)
# plotChromPeakDensity(mzdata, mz = mzr, col = cols, param = pdp,
#                      pch = 16, xlim = c(20, 100))

# plot the effects of the peakdensityparam
pdp <- PeakDensityParam(sampleGroups = mzdata$class,
                        minFraction = 0.5, 
                        bw = 5,
                        binSize = 0.025)
plotChromPeakDensity(mzdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(20, 100))

# perform grouping
pdp <- PeakDensityParam(mzdata$class,
                        bw = 5,
                        minFraction = 0.5,
                        binSize = 0.025)
mzdata <- groupChromPeaks(mzdata, param = pdp)

# quick look at the results
featureDefinitions(mzdata)

# fill in missing peaks
mzdata <- fillChromPeaks(mzdata)

# examine features
feature_chroms <- featureChromatograms(mzdata, features = 1:10)

# plot individual features
par(mfrow = c(2, 2))
plot(feature_chroms[1, ], col = class_colours[mzdata$class])
plot(feature_chroms[2, ], col = class_colours[mzdata$class])
plot(feature_chroms[3, ], col = class_colours[mzdata$class])
plot(feature_chroms[4, ], col = class_colours[mzdata$class])

# convert to xcmsSet object for CAMERA annotation
mzdata <- as(mzdata, "xcmsSet")

# Create report with with isotope and adduct info
isorep <- annotate(mzdata, 
                   perfwhm = 0.7, 
                   cor_eic_th = 0.75, 
                   ppm = 10, 
                   polarity = "positive")
mzdata <- getPeaklist(isorep)

# Save files
saveRDS(xset.final, "./data/raw/xset.final.v2.rds")
saveRDS(mzdata, "./data/raw/HeronMay13.mzdata.v2.rds")
write.csv(mzdata, "./data/raw/HeronMay13.mzdata.v2.csv")










