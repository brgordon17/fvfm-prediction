# script that returns a table of experiemental conditions
# Benjamin Gordon
# 11 JAN 2019
#
# Data checked and confirmed with field notes and temp loggers etc

library(tibble)

# Create variables

exp_cond_2014 <- 
  tibble(date = seq(as.Date("2014-03-15"), as.Date("2014-03-28"), by = 1),
         day = c(0:13),
         do_collect = as.logical(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1)),
         do_pam = as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1)),
         samp_time = c(NA, NA, NA, "12:30:00", NA, NA, NA, NA, NA, "12:30:00",
                       NA, "12:30:00", NA, "12:30:00"),
         # Average FvFm of 12 reps
         cont_FvFm = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.714, 0.666, 0.702, NA, 
                       NA, 0.664),
         # Average FvFm of 12 reps
         treat_FvFm = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.520, 0.511, 0.560, NA, 
                        NA, 0.600),
         # Average control tank temp at time of sampling
         cont_temp = (c(0, 0, 0, 27.17, 0, 0, 0, 0, 0, 27.96, 0, 26.59, 0, 
                        26.29)),
         # Average treatment tank temp at time of sampling
         treat_temp = (c(0, 0, 0, 28.95, 0, 0, 0, 0, 0, 33.33, 0, 26.68, 0, 
                         26.12))
  )
