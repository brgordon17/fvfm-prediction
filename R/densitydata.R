# Script to create cell density data. Data obtained from 
# ~Documents/Experiment Data/Heron bleaching/physiological data.docx
# Author: Benkamin R. Gordon
# 2019-02-27

library(tidyverse)

# Create data vector
densitydata <- tibble(day = c(1, 1, 5, 5, 8, 8, 10, 10, 12, 12, 15, 15),
                      class = rep(c("control", "treatment"), 6),
                      count = c(1171922, 1073466, 586582, 542247, 1777001, 
                                1730994, 1684769, 1308929, 1781139, 1325810, 
                                1610346, 697616)
                      )

# Save
save(densitydata, file = "./data/densitydata.rda", compress = "bzip2")
