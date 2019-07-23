# Code to construct the temperature profile plot

# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

# temperature profile plot -----------------------------------------------------
library(tidyverse)

# load data
load("./data/tempdata.rda")

plot_time <- c(as.POSIXct("2013-05-08", format = "%F"), 
               as.POSIXct("2013-05-22", format = "%F"))

# Plot
ggplot(data = tempdata, aes(x = time)) +
  geom_path(aes(y = meanC, colour = "Control")) +
  geom_path(aes(y = meanT, colour = "Heated")) +
  geom_path(aes(y = meanIMOS, colour = "Reef@0.3m")) +
  labs(title = NULL, x = "Time (days)", y = "Temperature"~(degree~C)) +
  scale_color_manual(values = c("Control" = gordon01::qual_colours[2],
                                "Heated" = gordon01::qual_colours[6],
                                "Reef@0.3m" = "grey80")) +
  scale_x_datetime(date_labels = as.character(c(0:15)), date_breaks = "1 days", 
                   limits = plot_time) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 10,
                                 colour = "grey65"),
        axis.title.y = element_text(size = 14),
        axis.title.x =  element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) 

# save
ggsave("./figs/temp_plot.pdf",
       width = 10,
       height = 6.5,
       units = "in")