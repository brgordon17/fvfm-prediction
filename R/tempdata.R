# Script to create tempdata.rda
# Author: Benjamin R. Gordon
# Date: 2019-03-11

library(tidyverse)

# load data
tempdata <- 
  read_csv("~/Documents/Experiment Data/Heron Bleaching_May2013/Temp Data/Temp Data_HeronMay13.csv")

imosdata <- 
  read_csv("~/Documents/Experiment Data/Heron Bleaching_May2013/Temp Data/imosdata_May13.csv")

# The time format changes middway through the df. 
# split df, correct time (POSIXct), and rbind
td1 <- tempdata[1:1686, ]
td2 <- tempdata[1687:3091, ] # dont include last row
td1$time <- as.POSIXct(td1$time, format = "%m/%d/%y %H:%M")
td2$time <- as.POSIXct(td2$time, format = "%m/%d/%y %I:%M:%S %p")
tempdata <- bind_rows(td1, td2)
rm(td1, td2)

# Create columns for average control and treatment temp and calibrate the data. 
# The data loggers read 1 degree higher than the in-tank bulb thermometers.
tempdata$meanC <- rowMeans(tempdata[, 2:5]) - 1
tempdata$meanT <- rowMeans(tempdata[, 6:9]) - 1

# Add imos data after cleaining up ---------------------------------------------
time <- as.POSIXct(paste(imosdata$date, imosdata$time), 
                   format = "%d/%m/%y %H:%M:%S")
imosdata[, 1:2] <- NULL
imosdata <- bind_cols(time = time, imosdata)
imosdata$meanIMOS <- rowMeans(imosdata[, 2:4])

# Remove rows with time occuring after 2013-05-22 18:01:34
imosdata <- imosdata[-6266:-6625, ]

# Remove rows up to (but not including) time of 2013-05-01 07:00:00
imosdata <- imosdata[-1:-84,]

# Remove all rows at five min intervals (i.e. 07:05, 07:15 ...).
r <- 1:nrow(imosdata)
head(r%%2) # labels every second row with a zero (modulo subsetting)
imosdata <- imosdata[r%%2==1 ,]
rm(r)

# bind data
tempdata <- bind_cols(tempdata, meanIMOS = imosdata$meanIMOS)

# Save -------------------------------------------------------------------------
save(tempdata, file = "./data/tempdata.rda", compress = "bzip2")
