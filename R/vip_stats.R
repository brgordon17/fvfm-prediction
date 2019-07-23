# Script to conduct stats anlaysis of annotated features.
# Author: Benjamin R Gordon.

library(tidyverse)

# Load and prep data
load("./data/mzdata.rda")

features <- 
  mzdata %>%
  filter(cont_treat != "PBQC") %>%
  droplevels() %>%
  select(class,
         day,
         cont_treat,
         FvFm,
         mz_249.18855,
         mz_170.04777,
         mz_985.69447)

# test for normality (normal if p > 0.05)
shapiro_data <- apply(features[, 5:ncol(features)], 2, shapiro.test)
shapiro_data <- unlist(lapply(shapiro_data, function(x) x$p.value))
shapiro_data

# all of the variables are non-normal so perform a Kruskal-Wallis test
# for 249
kruskal.test(features$mz_249.18855, features$day)

# for 170
kruskal.test(features$mz_170.04777, features$day)

# for 985
kruskal.test(features$mz_985.69447, features$day)
