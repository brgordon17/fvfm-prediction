# Script to conduct stats anlaysis of VIP features.
# Author: Benjamin R Gordon.

library(tidyverse)
library(caret)

# Load data
load("./data/mzdata.rda")
fvfm_model <- readRDS("./dev/mzrf_model_fvfm.rds")
cvst_model <- readRDS("./dev/mzrf_model_cvst.rds")

# tidy data
fvfm_impvars <- varImp(fvfm_model, scale = TRUE)
fvfm_impvars <- fvfm_impvars$importance
fvfm_impvars <- fvfm_impvars[order(-fvfm_impvars$Overall), ,drop = FALSE]
fvfm_impvars <- bind_cols(mz = rownames(fvfm_impvars), fvfm_impvars)

# top five vips
top5 <- fvfm_impvars[1:5, ][1]
top5.1 <- select(class, )

