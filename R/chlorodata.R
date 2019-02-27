# Script to create chl a data from 
# ~/Documents/Experimental Data/Heron bleaching/physiological data.docx.
# Author: Benjamin R. Gordon
# Date: 2019-02-27

library(tidyverse)

chlorodata <- tibble(day = c(1, 1, 5, 5, 8, 8, 10, 10, 12, 12, 15, 15),
                     class = rep(c("control", "treatment"), 6),
                     conc = c(0.317626855, 0.53331829, 0.609356481, 0.614901071, 
                              0.317545029, 0.413400328, 0.3679534, 0.551330265, 
                              0.434073463, 0.684514069, 0.456136747, 
                              0.932803169)
                     )

# save
save(chlorodata, file = "./data/chlorodata.rda", compress = "bzip2")
