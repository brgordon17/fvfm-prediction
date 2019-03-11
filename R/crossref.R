# Table that cross references litmz to provide putative IDs
# Author: Benjamin R. Gordon
# Data: 2019-03-11

library(tidyverse)
library(caret)

# Load data --------------------------------------------------------------------
mzrf_cvst <- readRDS("./dev/mzrf_model_cvst.rds")
mzrf_fvfm <- readRDS("./dev/mzrf_model_fvfm.rds")
load("./data/litmz.rda")
  
cvst_impvars <- varImp(mzrf_cvst, scale = TRUE)
cvst_impvars <- as_tibble(cvst_impvars$importance, rownames = "mz")
cvst_impvars <- 
  cvst_impvars %>%
  mutate(mz = gsub("mz_", "", mz)) %>%
  rowwise() %>%
  transmute(
    mz = as.numeric(mz), # warnings arise from features with two decimal points
    importance = max(T)) %>%
  arrange(desc(importance)) %>%
  slice(1:20)
  
fvfm_impvars <- varImp(mzrf_fvfm, scale = TRUE)
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
cvst_impvars <-
  cvst_impvars %>%
  rowwise() %>%
  mutate(
    model = "class",
    mz_neutral = mz - 1.007276,
    mz_low = mz_neutral - (mz * ppm/10^6),
    mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  select(model, everything()) %>%
  ungroup()

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

impvars <- bind_rows(cvst_impvars, fvfm_impvars)

# Cross reference with litmz ---------------------------------------------------
mz_matches <-
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high)

# NO HITS
