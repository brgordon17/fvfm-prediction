# script that returns a table of experiemental conditions
# Benjamin Gordon
# 10 JAN 2019
#
# Data checked and confirmed with field notes, temp data and Sarah's notes

library(tibble)

# Create variables

exp_cond_2013 <- 
  tibble(date = seq(as.Date("2013-05-06"), by = 1, len = 17),
         day = c(0:16),
         do_collect = as.logical(c(1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0,
                                   0, 1)),
         do_pam = as.logical(c(1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)),
         do_dens = as.logical(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1)),
         do_pigment = as.logical(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1)),
         samp_time = c("12:30:00", NA, NA, NA, NA, "12:45:00", NA, NA, 
                       "12:45:00", NA, "12:30:00", NA, "12:45:00", NA, NA, NA, 
                       "12:30:00"),
         # Average FvFm of 12 reps
         cont_FvFm = c(0.652, NA, NA, NA, NA, 0.657, NA, NA, 0.655, NA, 0.678, 
                        NA, 0.643, NA, NA, NA, 0.659),
         # Average FvFm of 12 reps
         treat_FvFm = c(0.646, NA, NA, NA, NA, 0.641, NA, NA, 0.575, NA, 0.573, 
                        NA, 0.304, NA, NA, NA, 0.000),
         cont_dens = c(1171922, NA, NA, NA, NA, 586582, NA, NA, 1777001, NA, 
                       1684769, NA, 1781139, NA, NA, NA, 1610346),
         treat_dens = c(1073466, NA, NA, NA, NA, 542247, NA, NA, 1730994, NA, 
                        1308929, NA, 1325810, NA, NA, NA, 697616),
         # Average control tank temp at time of sampling
         cont_temp = c(23.388, NA, NA, NA, NA, 23.004, NA, NA, 24.859, NA,
                       22.932, NA, 24.545, NA, NA, NA, 23.340),
         # Average treatment tank temp at time of sampling
         treat_temp = c(25.125, NA, NA, NA, NA, 30.0035, NA, NA, 32.911, NA, 
                        32.240, NA, 34.322, NA, NA, NA, 33.066)
         )
