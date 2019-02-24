# Script to build a dataframe of the PAM data
# Author: Benjamin R. Gordon
# Date: 22-FEB-

library(tidyverse)

# load all files into a list
fileidx <- 
  list.files(path = "~/Documents/Experiment Data/Heron Bleaching_May2013/PAM Data/",
             pattern = "*.csv",
             full.names = TRUE)
pamdata <- lapply(fileidx, read_delim, delim = ";")
pamdata <- lapply(pamdata, type_convert)
pamdata[1]

# delete all data in each df but keep Y(II)
pamdata <- lapply(pamdata, 
                  function (x) {select(x, Date, Time, starts_with("Y(II)"))})

# keep only the first row of each df
pamdata <- lapply(pamdata, function(x) x[-2:-nrow(x), ])

# convert list to a single df
pamdata <- bind_rows(pamdata)

# add variables for class, day and collect (logical)
class <- factor(c(rep("control", 15),
                rep("treatment", 15))
                )
day <- factor(c(sprintf("day %d", seq(1:15)),
                sprintf("day %d", seq(1:15))
                ),
              levels = sprintf("day %d", seq(1:15))
              )
collect <- day %in% c("day 1", "day 5", "day 8", "day 10", "day 12",
                      "day 15")

pamdata <- bind_cols(class = class, day = day, collect = collect, pamdata)

# merge date and time columns
pamdata$Date <- gsub("\\.", "-", pamdata$Date)
colnames(pamdata)[4] <- "date"
colnames(pamdata)[5] <- "time"
pamdata$date <- as.POSIXct(paste(pamdata$date, pamdata$time), 
                           format = "%d-%m-%y %H:%M:%S")
pamdata$time <- NULL

# remove days where sampling wasnt performed
#pamdata <- filter(pamdata, day %in% c("day 1", "day 5", "day 8", "day 10", 
#                    "day 12", "day 15"))

# gather df to long format
pamdata <- gather(pamdata, key = "rep", value = "yield", "Y(II)1":"Y(II)12")
pamdata$rep <- gsub("Y\\(II\\)", "rep ", pamdata$rep)

# save
save(pamdata, file = "./data/pamdata.rda", compress = "bzip2")
