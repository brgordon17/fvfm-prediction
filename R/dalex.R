# Script to run DALEX analysis
# Author: BEnjamin R. Gordon
# Date: 2019-03-12

library(caret)
library(DALEX)
library(tidyverse)

# Load and prep data -----------------------------------------------------------
load("./data/mzdata.rda")
mzdata <- filter(mzdata, cont_treat != "PBQC")
mzdata <- droplevels(mzdata)
set.seed(1978)
train_index <- createDataPartition(mzdata$FvFm,
                                   p = 0.8,
                                   list = FALSE)
test_data  <- mzdata[-train_index, ]

mzrf_fvfm <- readRDS("./dev/mzrf_model_fvfm.rds")

# Create explainer
explain_fvfm <- explain(mzrf_fvfm, 
                        label = "rf", 
                        data = test_data, 
                        y = test_data$FvFm)

# model performance
mp_fvfm <- model_performance(explain_fvfm)
mp_fvfm
plot(mp_fvfm)

# Variable importance
vi_fvfm <- variable_importance(explain_fvfm, 
                               loss_function = loss_root_mean_square)
plot(vi_fvfm)

# Partial dependance plot
pdp_fvfm  <- variable_response(explain_fvfm, 
                               variable = "mz_170.04765", 
                               type = "pdp")
plot(pdp_fvfm)

ggplot(pdp_fvfm, aes(x = y, y = x)) +
  geom_path()
