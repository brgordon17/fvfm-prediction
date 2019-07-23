# Code to construct a table of annotated features matching literature (coralmz)
# Author: Benjamin R. Gordon
# Date: 2019-03-17

library(tidyverse)
library(caret)

# Load and prep data -----------------------------------------------------------
mzrf_fvfm <- readRDS("./dev/mzrf_model_fvfm.rds")
coralmz <- coralmz::coralmz
mzdata_raw  <-  readr::read_csv("./data-raw/mzdata-raw.csv", na = "0")
colnames(mzdata_raw)[1] <- "mz_raw"

# cross reference 20 impvars ---------------------------------------------------
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

# Add variables for 50 ppm error ranges 
ppm <- 50
mz_matches <-
  fvfm_impvars %>%
  rowwise() %>%
  mutate(adduct = NA,
         mz_neutral = mz - 1.007276,
         mz_low = mz_neutral - (mz * ppm/10^6),
         mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with coralmz
mz_matches <-
  mz_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(coralmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance)

# cross reference any adducts --------------------------------------------------
# identify adducts from impvars
adduct_matches <- 
  fvfm_impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(mzdata_raw %>% mutate(dummy = TRUE))  %>%
  filter(near(mz, mz_raw, tol = .0001)) %>%
  select(-dummy,
         -5:-182,
         -pcgroup,
         -mz_raw) %>%
  mutate(mz_neutral = as.numeric(str_sub(adduct, str_length(adduct)-6, -1))) %>%
  mutate(adduct = str_sub(adduct, 1, str_length(adduct)-8)) %>%
  filter(!is.na(mz_neutral))

# add variables for 50ppm error
ppm <- 50
adduct_matches <-
  adduct_matches %>%
  rowwise() %>%
  mutate(mz_low = mz_neutral - (mz_neutral * ppm/10^6),
         mz_high = mz_neutral + (mz_neutral * ppm/10^6)) %>%
  ungroup()

# Cross reference with coralmz
adduct_matches <-
  adduct_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(coralmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -isotopes,
         -importance)

# Construct table --------------------------------------------------------------
matches <- bind_rows(mz_matches, adduct_matches)

# remove generic referencing and replace commas
matches <-
  matches %>%
  mutate(endnote_ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
  select(-ref)


# Save csv 
readr::write_csv(matches, "./tables/important_variable_matches.txt")
