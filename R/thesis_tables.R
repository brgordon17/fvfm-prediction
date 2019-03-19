# Script to construct tables for thesis
# Author: Benjamin R. Gordon
# Date: 2019-03-17

library(tidyverse)

# Load and prep data
matches <- read_rds("./dev/impvar_matches.rds")

# remove generic referencing and replace commas
matches <-
  matches %>%
  mutate(endnote_ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
  select(-ref)


# Save csv ---------------------------------------------------------------------
readr::write_csv(matches, "./tables/important_variable_matches.txt")
