# Table that cross references litmz to provide putative IDs
# Author: Benjamin R. Gordon
# Data: 2019-03-11

library(tidyverse)
library(caret)

# Load data --------------------------------------------------------------------
mzrf_fvfm <- readRDS("./dev/mzrf_model_fvfm.rds")
load("./data/litmz.rda")

fvfm_impvars <- varImp(mzrf_fvfm, scale = FALSE)
fvfm_impvars <- as_tibble(fvfm_impvars$importance, rownames = "mz")
fvfm_impvars <-
  fvfm_impvars %>%
  mutate(mz = gsub("mz_", "", mz)) %>%
  rowwise() %>%
  transmute(
    mz = as.numeric(mz),# warnings arise from features with two decimal points
    importance = max(Overall)) %>%
  arrange(desc(importance)) %>%
  slice(1:20)
  
# Add variables for 50 ppm error ranges ----------------------------------------
ppm <- 50
fvfm_impvars <-
  fvfm_impvars %>%
  rowwise() %>%
  mutate(
    model = "regr",
    mz_neutral = mz - 1.007276,
    mz_low = mz_neutral - (mz * ppm/10^6),
    mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  select(model, everything()) %>%
  ungroup()

# Cross reference with litmz ---------------------------------------------------
mz_matches <-
  fvfm_impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high)
