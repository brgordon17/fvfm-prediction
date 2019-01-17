# Script to perform random forests  predicition of FvFm
# Ben Gordon 16JAN19

library(tidyverse)
library(caret)

# data
load("./data/mzdata.rda")

mzdata2013 <-
  mzdata %>%
  filter(experiment == "2013" & cont_treat != "PBQC")
mzdata2013 <- droplevels(mzdata2013)

# Check for near zero variance
nearZeroVar(mzdata2013[-1:-7])

# Check for highly correlated predictors (mz variables)
# Consider removing these later to see if model improves
cormat <- cor(mzdata2013[-1:-7])
highlyCorDescr <-findCorrelation(cormat)
rm(cormat)

# Partition data into training, test and validation sets
train_index <- createDataPartition(mzdata2013$FvFm,
                                   p = 0.8,
                                   list = FALSE)
train_data <- mzdata2013[train_index, ]
test_data  <- mzdata2013[-train_index, ]
validation_data <- 
  mzdata %>%
  filter(experiment != "2013" & FvFm != "NA")

# tuning grid
# Consider changing from sqrt(p) to p/3 (used for regression)
tunegrid <- expand.grid(.mtry = c(25, 50, 75, 100, 200, 300, 400, 500, 1000,
                                  2000, floor(length(train_data)/3)))

# set seeds
set.seed(1978)
seeds <- vector(mode = "list", length = 31)
for(i in 1:30) seeds[[i]] <- sample.int(1000, length(tunegrid[,1]))
seeds[[31]] <- sample.int(1000, 1)

# train control
ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 3,
                            summaryFunction = defaultSummary,
                            seeds = seeds,
                            savePredictions = "all",
                            selectionFunction = "oneSE")

# Create model
doMC::registerDoMC()
set.seed(1978)
mzrf <- caret::train(x = train_data[, -1:-7],
                     y = train_data$FvFm,
                     method = "rf",
                     trControl = ctrl,
                     preProc = c("center", "scale"),
                     allowParallel = TRUE,
                     importance = TRUE,
                     tuneGrid = tunegrid
                     )

# predict test data
test_pred <- predict(mzrf, newdata = test_data)

# test data prediction metrics
postResample(pred = test_pred, obs = test_data$FvFm)

# save as temporary file
saveRDS(mzrf, "./dev/mzrf_model.rds")
