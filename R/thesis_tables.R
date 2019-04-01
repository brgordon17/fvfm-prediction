# Script to construct tables for thesis
# Author: Benjamin R. Gordon
# Date: 2019-03-17

# Table of compounds matching literature ---------------------------------------
library(tidyverse)

# Load and prep data
matches <- read_rds("./dev/impvar_matches.rds")

# remove generic referencing and replace commas
matches <-
  matches %>%
  mutate(endnote_ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
  select(-ref)


# Save csv 
readr::write_csv(matches, "./tables/important_variable_matches.txt")

# Table of predictions for unseen samples --------------------------------------
library(tidyverse)
library(caret)

# Load and prep data
load("./data/mzdata.rda")
reg_model <- read_rds("./dev/mzrf_model_fvfm.rds")

# test data
mzdata <- filter(mzdata, cont_treat != "PBQC")
mzdata <- data.frame(droplevels(mzdata))
set.seed(16)
train_index <- createDataPartition(mzdata$FvFm,
                                   p = 0.9,
                                   list = FALSE)
test_data  <- mzdata[-train_index, ]

# predict test data
test_pred <- predict(reg_model, newdata = test_data)

# test data prediction metrics
postResample(pred = test_pred, obs = test_data$FvFm)

# Compare pred vs actual
preds <- bind_cols(day = test_data$day,
                   treatment = test_data$cont_treat,
                   actual = round(test_data$FvFm, 3),
                   predicted = round(test_pred, 3)
                   )
readr::write_csv(preds, "./tables/unseen_predictions.txt")

# remove day samples from predictions
test_data_no8 <- filter(test_data, day != "day 8")
test_pred_no8 <- predict(reg_model, newdata = test_data_no8)
postResample(pred = test_pred_no8, obs = test_data_no8$FvFm)
