# Script to perform random forests  predicition of FvFm
# Ben Gordon 16JAN19

library(tidyverse)
library(caret)
library(gmailr)

# data
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

# Partition data for multiple k-fold cv's (not used at this stage)
# set.seed(1978)
# train_index <- createDataPartition(mzdata$FvFm,
#                                    p = 0.8,
#                                    times = 3,
#                                    list = FALSE)
# test_index <- apply(train_index, 2, function(x) setdiff(1:141, x))

# Create model for FvFm prediction ---------------------------------------------
# tuning grid p/3 (used for regression)
tunegrid <- expand.grid(.mtry = c(25, 50, 75, 100, 200, 300, 400, 500, 1000,
                                  floor(length(train_data)/3)))

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

# model
doMC::registerDoMC()
set.seed(1978)
mzrf_fvfm <- caret::train(x = train_data[, -1:-6],
                     y = train_data$FvFm,
                     method = "rf",
                     trControl = ctrl,
                     #preProc = c("center", "scale"),
                     allowParallel = TRUE,
                     importance = TRUE,
                     tuneGrid = tunegrid
                     )
# email notification
use_secret_file("~/Documents/R/R_emails/R-emails.json")
send_message(mime(To = "benjamin.gordon@my.jcu.edu.au",
                  From = "brgordon17@gmail.com",
                  Subject = "Analysis complete",
                  body = "Your R process has finished."))

# predict test data
test_pred <- predict(mzrf_fvfm, newdata = test_data)

# test data prediction metrics
postResample(pred = test_pred, obs = test_data$FvFm)

# Compare pred vs actual
preds <- bind_cols(id = test_data$sample_id,
                   pred = test_pred, 
                   actual = test_data$FvFm)

# save as temporary file
saveRDS(mzrf_fvfm, "./dev/mzrf_model_fvfm.rds")

# model for cont_treat classification ------------------------------------------
# tuning grid sqrt(p) (used for classification)
tunegrid <- expand.grid(.mtry = c(25, 50, 75, 100))

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

# model
doMC::registerDoMC()
set.seed(1978)
mzrf_cvst <- caret::train(x = train_data[, -1:-6],
                          y = train_data$cont_treat,
                          method = "rf",
                          trControl = ctrl,
                          #preProc = c("center", "scale"),
                          allowParallel = TRUE,
                          importance = TRUE,
                          tuneGrid = tunegrid
)
# email notification
use_secret_file("~/Documents/R/R_emails/R-emails.json")
send_message(mime(To = "benjamin.gordon@my.jcu.edu.au",
                  From = "brgordon17@gmail.com",
                  Subject = "Analysis complete",
                  body = "Your R process has finished."))

# predict test data
test_pred <- predict(mzrf_cvst, newdata = test_data)

# test data prediction metrics
confusionMatrix(test_pred, test_data$cont_treat)

# Compare pred vs actual
preds <- bind_cols(id = test_data$sample_id,
                   pred = test_pred, 
                   actual = test_data$cont_treat)

# save as temporary file
saveRDS(mzrf_cvst, "./dev/mzrf_model_cvst.rds")
