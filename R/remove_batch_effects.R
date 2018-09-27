# Script to remove batch effects from mzdata using he Harman package. Conducted
# by B. Gordon on 27SEP18

# Load Libraries
library(tidyverse)
library(Harman)

# Using log transformed data ---------------------------------------------------

# Load data and transform
load("./data/mzdata.rda")
harmdata <- as.data.frame(t(mzdata[-1:-7]))
harmdata <- log(harmdata + 1, 2)
colnames(harmdata) <- mzdata$sample_ids

# Run Harman
expt <- mzdata$class_day
batch <- mzdata$batch
mzharm <- harman(harmdata, expt, batch,
                 limit = 0.95)

# Explore harman results and compare
summary(mzharm)
plot(mzharm)

# Reconstruct the data and compare
mzdata_hm <- reconstructData(mzharm)
par(mfrow = c(1, 2))
prcompPlot(harmdata, 1, 2, colFactor = batch, cex = 1.5, main = "Original")
prcompPlot(mzdata_hm, 1, 2, colFactor = batch, cex = 1.5, main = "Corrected")

# Using raw data ---------------------------------------------------------------

# Load data and transform
load("./data/mzdata.rda")
harmdata <- as.data.frame(t(mzdata[-1:-7]))
colnames(harmdata) <- mzdata$sample_ids

# Run Harman
expt <- mzdata$class_day
batch <- mzdata$batch
mzharm <- harman(harmdata, expt, batch,
                 limit = 0.95)

# Explore harman results and compare
summary(mzharm)
plot(mzharm)

# Reconstruct the data and compare
mzdata_hm <- reconstructData(mzharm)
par(mfrow = c(1, 2))
prcompPlot(harmdata, 1, 2, colFactor = batch, cex = 1.5, main = "Original")
prcompPlot(mzdata_hm, 1, 2, colFactor = batch, cex = 1.5, main = "Corrected")


