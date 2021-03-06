# Script to run RM ANOVA and post-hoc analysis on PSII yield

library(tidyverse)
library(car)

# load data
load("./data/pamdata.rda")

# summary stats
summary(pamdata$yield[pamdata$class == "control"])
summary(pamdata$yield[pamdata$class == "treatment"])

# RMANOVA ----------------------------------------------------------------------
pam_aov <- aov(yield ~ class * day, data = pamdata)
summary(pam_aov)
# report as (RMANOVA, F -df_resid, -df_class:day = Fval, p<0.001)

# check the residuals are normally distributed
res <- pam_aov$residuals
hist(res, main = "Histogram of residuals", xlab = "Residuals")

# Tukey's HSD ------------------------------------------------------------------
# check for equality of variance
# Use lower p threshold if equality of variance isnt true (p < 0.05)
leveneTest(yield ~ class, data = pamdata)

# Tukey hsd to see which groups differ
TukeyHSD(pam_aov)

# pull out control vs treatment for each sampling point
pam_tukey <- as_tibble(TukeyHSD(pam_aov)$`class:day`, rownames = "interaction")
pam_tukey <- filter(pam_tukey, interaction == "treatment:day 1-control:day 1" |
                      interaction == "treatment:day 2-control:day 2" |
                      interaction == "treatment:day 3-control:day 3" |
                      interaction == "treatment:day 4-control:day 4" |
                      interaction == "treatment:day 5-control:day 5" |
                      interaction == "treatment:day 6-control:day 6" |
                      interaction == "treatment:day 7-control:day 7" |
                      interaction == "treatment:day 8-control:day 8" |
                      interaction == "treatment:day 9-control:day 9" |
                      interaction == "treatment:day 10-control:day 10" |
                      interaction == "treatment:day 11-control:day 11" |
                      interaction == "treatment:day 12-control:day 12" |
                      interaction == "treatment:day 13-control:day 13" |
                      interaction == "treatment:day 14-control:day 14" |
                      interaction == "treatment:day 15-control:day 15")
