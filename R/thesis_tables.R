# Script to construct tables for thesis
# Author: Benjamin R. Gordon
# Date: 2019-03-17

library(tidyverse)
library(caret)

# Load and prep data
mzrf_fvfm <- read_rds("./dev/mzrf_model_fvfm.rds")
load("./data/litmz.rda")

impvars <- varImp(mzrf_fvfm, scale = FALSE)
impvars <- as_tibble(impvars$importance, rownames = "mz")
impvars <-
  impvars %>%
  mutate(mz = gsub("mz_", "", mz)) %>%
  rowwise() %>%
  transmute(
    mz = as.numeric(mz),# warnings arise from features with two decimal points
    importance = max(Overall)) %>%
  arrange(desc(importance)) %>%
  slice(1:20)

# Add variables for ppm error ranges -------------------------------------------
ppm <- 50
impvars <-
  impvars %>%
  rowwise() %>%
  mutate(
    mz_neutral = mz - 1.007276,
    mz_low = mz_neutral - (mz * ppm/10^6),
    mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with litmz -------------------------------------------------
matches <-
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high)

# remove generic referencing and replace commas
matches <-
  matches %>%
  mutate(endnote_ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
  select(-ref)


# Save csv ---------------------------------------------------------------------
readr::write_csv(matches, "./tables/important_variable_matches.txt")
