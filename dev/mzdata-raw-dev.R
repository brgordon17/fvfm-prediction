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

# check peak width by selecting a peak that is common to all samples as IDd in 
# the BPC and mMass.app
# define rt range
rtr <- c(20, 100)
mzr <- c(135.1, 136.9)
# extract chromatogram
chr_raw <- chromatogram(raw_data, rt = rtr, mz = mzr)
plot(chr_raw, col = class_colours[raw_data$class])

# check ppm of above peak
raw_data %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")

# Peak detection ---------------------------------------------------------------

# define the parameter used (centwave here)
cwp <- CentWaveParam(ppm = 30, 
                     peakwidth = c(20, 60),
                     noise = 1000)
mzdata <- findChromPeaks(raw_data, param = cwp)

# view peak counts
summary_fun <- function(z) {
  c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
}
peak_summary <- lapply(split.data.frame(chromPeaks(mzdata),
                             f = chromPeaks(mzdata)[, "sample"]),
            FUN = summary_fun)
peak_summary <- do.call(rbind, peak_summary)
rownames(peak_summary) <- basename(fileNames(mzdata))
View(peak_summary)

# ensure previously examined peak was picked correctly
plot(chr_raw, col = class_colours[raw_data$class], lwd = 2)
highlightChromPeaks(mzdata, 
                    border = class_colours[raw_data$class],
                    lty = 3,
                    rt = rtr,
                    mz = mzr)








