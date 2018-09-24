# Header####
#
# Title: Experimental Conditions and times_HeronMay13
# Date: 13-NOV-2017
# Author: Benjamin R. Gordon
# Owner:Benjamin R. Gordon
# Copyright:Benjamin R. Gordon
# Version: 1.0
# Version of R used: 1.0.143
# Description: Script to create df of conditions and times

# clear global environment
rm(list=ls(all=T))

# Create vectors for day, average set temp, collection (y/n)
day <- c(0:16)
date <- seq(as.Date("2013-05-06"), as.Date("2013-05-22"), by="days")
avetemp <- c(seq(from = 25, to = 32.3, by = 0.7), 33, 33, 33, 34, 34, 34)
collect <- c(1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1)
expcond <- data.frame(day, date, avetemp, collect)

# save
saveRDS(expcond, "./data/processed/expcond.rds")
