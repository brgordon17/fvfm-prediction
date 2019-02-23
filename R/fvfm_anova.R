# Script to run anovas on PSII yield

library(tidyverse)

# load data
load("./data/pamdata.rda")

# students t-test
pam_tests <- list(
  "day_1" = t.test(yield ~ class,
                   data = filter(pamdata, day == "day 1")),
  "day_5" = t.test(yield ~ class,
                   data = filter(pamdata, day == "day 5")),
  "day_8" = t.test(yield ~ class,
                   data = filter(pamdata, day == "day 8")),
  "day_10" = t.test(yield ~ class,
                    data = filter(pamdata, day == "day 10")),
  "day_12" = t.test(yield ~ class,
                    data = filter(pamdata, day == "day 12")),
  "day_15" = t.test(yield ~ class,
                    data = filter(pamdata, day == "day 15"))
  )
str(pam_tests)

# construct df
pam_tests <- data.frame(pval = sapply(pam_tests, "[[", "p.value"),
                        t_stat = sapply(pam_tests, "[[", "statistic"),
                        df = sapply(pam_tests, "[[", "parameter")
                        )
