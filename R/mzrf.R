# Script to perform random forests  predicition of FvFm
# Ben Gordon 16JAN19

# data
load("./data/mzdata.rda")

mzdata2013 <-
  mzdata %>%
  filter(experiment == "2013" & cont_treat != "PBQC")
mzdata2013 <- droplevels(mzdata2013)

