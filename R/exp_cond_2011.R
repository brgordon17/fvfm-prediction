# script that returns a table of experiemental conditions
# Benjamin Gordon
# 11 JAN 2019
#
# Data checked and confirmed with field notes, temp data and Daisie's notes

library(tibble)

# Create variables

exp_cond_2011 <- 
  tibble(date = seq(as.Date("2011-02-09"), by = 1, len = 14),
         day = c(0:13),
         do_collect = as.logical(c(1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1)),
         do_pam = as.logical(c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)),
         do_dens = as.logical(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1)),
         do_pigment = as.logical(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1)),
         # Average FvFm of 12 reps
         cont_FvFm = c(0.693, NA, NA, 0.6125, NA, 0.610, NA, NA, NA, NA, NA, 
                       NA, NA, 0.720),
         # Average FvFm of 12 reps
         treat_FvFm = c(0.635, NA, NA, 0.578, NA, 0.523, NA, NA, NA, NA, NA, 
                        NA, NA, 0.502),
         # Average control tank temp at time of sampling
         cont_temp = c(30.3, NA, NA, 28.5, NA, 30.5, NA, NA, NA, NA, NA, 
                       NA, NA, 29.1),
         # Average treatment tank temp at time of sampling
         treat_temp = c(34.0, NA, NA, 34.5, NA, 32.8, NA, NA, NA, NA, NA, 
                        NA, NA, 30.7)
  )
