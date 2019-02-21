# Script that contains code for all figures used in chapter 5 (bleaching) of
# thesis.
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)

# PSII FvFm line plot ----------------------------------------------------------

# load data
load("./data/metadata.rda")
metadata <- filter(metadata, cont_treat != "PBQC")
metadata$day <- factor(metadata$day,
                       levels = c("day 1", "day 5", "day 8", "day 10", "day 12",
                                  "day 15"))

# create plot
ggplot(metadata, aes(x = day,
                     y = FvFm,
                     colour = cont_treat,
                     group = cont_treat)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_discrete(name = "Day") +
  scale_y_continuous(name = "Quantum Yield (FvFm)") +
  scale_colour_manual(values = gordon01::qual_colours, 
                      name = NULL,
                      labels = c("control", "treatment")) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.6),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_text(size = 12, colour = "grey65"),
        axis.ticks = element_blank(),
        legend.key = element_blank()
        )
ggsave("./figs/fvfm_plot.pdf",
       width = 10,
       height = 6.5,
       units = "in")









