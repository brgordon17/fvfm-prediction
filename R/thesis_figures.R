# Script that contains code for all figures used in chapter 5 (bleaching) of
# thesis.
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)

# PSII FvFm plot ---------------------------------------------------------------

# load data
load("./data/pamdata.rda")

# summarise data
sum_pamdata <- 
  pamdata %>% 
  group_by(class, day) %>%
  summarise(mean = mean(yield), 
            sd = sd(yield), 
            se = sd(yield)/sqrt(n())
            )

# create plot
ggplot(sum_pamdata, aes(x = day,
                        y = mean,
                        colour = class,
                        shape = class)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = .1,
                show.legend = FALSE) +
  geom_point(size = 3, show.legend = TRUE) +
  scale_x_discrete(name = "Time (days)",
                   labels = c(1:15)) +
  scale_y_continuous(name = expression(Quantum~Yield~PSII~(F[v]/F[m]))) +
  scale_colour_manual(values = gordon01::qual_colours, 
                      name = NULL,
                      labels = c("control", "heated")) +
  scale_shape_manual(values = c(16, 17),
                     name = NULL,
                     labels = c("control", "heated")) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.6),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_text(size = 12, colour = "grey65"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 12)
        )
ggsave("./figs/fvfm_plot.pdf",
       width = 10,
       height = 6.5,
       units = "in")
