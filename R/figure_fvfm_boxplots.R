# Code to construct the fvfm boxplot
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)

# Load and prep data
reg_model <- read_rds("./dev/mzrf_model_fvfm.rds")

cv_preds <- filter(reg_model$pred, mtry == 2113)
cv_preds$obs <- factor(round(cv_preds$obs, 3))
cv_preds$pred <- round(cv_preds$pred, 3)

ggplot(cv_preds, aes(x = obs, y = pred)) +
  geom_boxplot() +
  xlab(Observed~Mean~Quantum~Yield~(F[v]/F[m])) +
  ylab(Predicted~Quantum~Yield~(F[v]/F[m])) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14)
  )

ggsave("./figs/prediction_boxplot.pdf",
       width = 10,
       height = 6.5,
       units = "in")