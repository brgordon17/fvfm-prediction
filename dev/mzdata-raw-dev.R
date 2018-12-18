# This is a developmental script for creating mzdata using the new xcms
# functions. 16-12-18

library(xcms)
library(magrittr)

# Data import ------------------------------------------------------------------

# path to files
mzfiles <- list.files(path = "./dev/test-data",
                     full.names = TRUE,
                     recursive = TRUE)

# create metadata or phenodata
pd <- data.frame(sample_name = sub(basename(mzfiles), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 class = c(rep("eCO2eT", 4), rep("recovery", 4)),
                 stringsAsFactors = FALSE)

# load raw data
raw_data <- readMSData(files = mzfiles, pdata = new("NAnnotatedDataFrame", pd),
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
rtr <- c(20, 100)
mzr <- c(136.05, 136.07)
# extract chromatogram
chr_raw <- chromatogram(raw_data, rt = rtr, mz = mzr)
plot(chr_raw, col = class_colours[raw_data$class])

# check ppm of above peak
raw_data %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")

# Peak detection ---------------------------------------------------------------

# define the parameters used for peak detection (centwave here)
cwp <- CentWaveParam(ppm = 30, 
                     peakwidth = c(20, 50),
                     snthresh = 10,
                     prefilter = c(3, 100),
                     mzCenterFun = "wMean",
                     integrate = 1L,
                     mzdiff = -0.001,
                     fitgauss = FALSE,
                     noise = 1000)

# Do peak detection
mzdata <- findChromPeaks(raw_data, param = cwp)

# view peaks given an mz and rt range
chromPeaks(mzdata, mz = c(136, 137), rt = rtr)

# view peaks given an mz and a ppm error
chromPeaks(mzdata, mz = 136.1076, ppm = 10 )

# Alignment --------------------------------------------------------------------

mzdata <- adjustRtime(mzdata, param = ObiwarpParam(binSize = 0.5))

# view unadjusted retention times
head(rtime(mzdata, adjusted = FALSE))

# and adjusted
head(rtime(mzdata))

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

# Correspondence or Grouping ---------------------------------------------------

# define mz range
mzr <- c(138.0, 138.1)

# Extract and plot the chromatograms
chr_mzr <- chromatogram(mzdata, mz = mzr, rt = c(20, 300))
cols <- class_colours[chr_mzr$class]
plot(chr_mzr, col = cols)

# Highlight the detected peaks in that region.
highlightChromPeaks(mzdata, mz = mzr, col = cols, type = "point", pch = 16)

# Define the parameters for the peak density method
pdp <- PeakDensityParam(sampleGroups = mzdata$class,
                        minFraction = 0.4, bw = 30)
plotChromPeakDensity(mzdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(20, 100))

# Use a different bw
pdp <- PeakDensityParam(sampleGroups = mzdata$class,
                        minFraction = 0.4, bw = 5)
plotChromPeakDensity(mzdata, mz = mzr, col = cols, param = pdp,
                     pch = 16, xlim = c(20, 100))

# perform grouping
pdp <- PeakDensityParam(mzdata$class,
                        bw = 5,
                        minFraction = 0.5,
                        binSize = 0.05)
mzdata <- groupChromPeaks(mzdata, param = pdp)

# quick look at the results
featureDefinitions(mzdata)

# fill in missing peaks
mzdata <- fillChromPeaks(mzdata)

# get data











