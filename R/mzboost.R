# Script for running boosted trees using caret and xgboost
# Benjamin Gordon 18JAN19

library(caret)
library(tidyverse)

# Data -------------------------------------------------------------------------
load("./data/mzdata.rda")

mzdata2013 <-
  mzdata %>%
  filter(experiment == "2013" & cont_treat != "PBQC")
mzdata2013 <- data.frame(droplevels(mzdata2013))

# Partition data into training, test and validation sets
train_index <- createDataPartition(mzdata2013$FvFm,
                                   p = 0.8,
                                   list = FALSE)
train_data <- mzdata2013[train_index, ]
test_data  <- mzdata2013[-train_index, ]
validation_data <- 
  mzdata %>%
  filter(experiment != "2013" & FvFm != "NA")
validation_data <- data.frame(droplevels(validation_data))

# Model ------------------------------------------------------------------------

ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 3,
                            search = "grid")

doMC::registerDoMC()
set.seed(1978)
mzboost <- caret::train(x = train_data[, -1:-7],
                        y = train_data$FvFm,
                        method = "xgbTree",
                        trControl = ctrl,
                        #preProc = c("center", "scale"), # made zero difference
                        allowParallel = TRUE,
                        importance = TRUE,
                        nthread = 1 # stop parallel proc in xgboost
                        )

test_pred <- predict(mzboost, newdata = test_data)

# test data prediction metrics
postResample(pred = test_pred, obs = test_data$FvFm)

# save as temporary file
saveRDS(mzboost, "./dev/mzboost_model.rds")

# mzrf <- read_rds("./dev/mzrf_model.rds")


