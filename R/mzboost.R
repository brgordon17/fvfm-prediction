# Script for running boosted trees using caret and xgboost
# Benjamin Gordon 18JAN19

library(caret)
library(tidyverse)

# Data -------------------------------------------------------------------------
load("./data/mzdata.rda")

mzdata <- filter(mzdata, cont_treat != "PBQC")
mzdata <- data.frame(droplevels(mzdata))

# Partition data into training, test and validation sets
set.seed(1978)
train_index <- createDataPartition(mzdata$FvFm,
                                   p = 0.8,
                                   list = FALSE)
train_data <- mzdata[train_index, ]
test_data  <- mzdata[-train_index, ]

# Model optimisation with auto grid --------------------------------------------

ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 3,
                            search = "grid")

doMC::registerDoMC()
set.seed(1978)
mzboost <- caret::train(x = train_data[, -1:-6],
                        y = train_data$FvFm,
                        method = "xgbTree",
                        trControl = ctrl,
                        #preProc = c("center", "scale"), # made zero difference
                        allowParallel = TRUE,
                        importance = TRUE,
                        nthread = 1 # stop parallel proc in xgboost
                        )

test_pred <- predict(mzboost, newdata = test_data)

# Check pred vs actual
preds <- bind_cols(id = test_data$sample_id,
                  pred = test_pred, 
                  actual = test_data$FvFm)

# test data prediction metrics
postResample(pred = test_pred, obs = test_data$FvFm)

# save as temporary file
saveRDS(mzboost, "./dev/mzboost_model.rds")
