# Code to construct the important variables plot
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)
library(caret)

# Load model 
mzrf_fvfm <- readRDS("./dev/mzrf_model_fvfm.rds")

# Identify important variables 
fvfm_impvars <- varImp(mzrf_fvfm, type = 1, scale = FALSE)
fvfm_impvars <- fvfm_impvars$importance
fvfm_impvars <- cbind(vip = apply(fvfm_impvars, 1, max), fvfm_impvars)
fvfm_impvars <- fvfm_impvars[order(-fvfm_impvars$vip), ,drop = FALSE]

# clean up dfs and convert to tibbles
rownames(fvfm_impvars) <- gsub("mz_", "", rownames(fvfm_impvars))
fvfm_impvars <- cbind(mz = factor(rownames(fvfm_impvars),
                                  levels = rev(rownames(fvfm_impvars))
), fvfm_impvars)
fvfm_impvars <- as_tibble(fvfm_impvars[1:20, ], rownames = NULL)

# create plots 
fvfm_plot <- ggplot(fvfm_impvars,
                    aes(x = vip, y = mz)) +
  geom_point(shape = 16,
             colour = phdhelpr::warm_colours[4],
             size = 3) +
  scale_x_continuous(name = "Mean Decrease in Accuracy") +
  scale_y_discrete(name = "m/z") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "italic"),
        axis.title.x = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.position = "none")
fvfm_plot

ggsave("./figs/vip_plot.pdf",
       width = 10,
       height = 6.5,
       units = "in")
