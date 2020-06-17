# Code to construct the cell density and chl-a plot
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)

# load data
load("./data/densitydata.rda")
load("./data/chlorodata.rda")

# create cell density plot
dens <- ggplot(filter(densitydata, day != "1" & day != "5"),
               aes(x = day,
                   y = count/10^6,
                   colour = class,
                   shape = class)) +
  geom_point(size = 3, show.legend = TRUE) +
  geom_line(show.legend = FALSE) +
  scale_x_continuous(name = "Time (days)", 
                     breaks = c(8, 10, 12, 15), 
                     labels = c(8, 10, 12, 15)) +
  scale_y_continuous(name = expression(Cell~Density~(10^6~cells~cm^-2)),
                     limits = c(0, 2)) +
  scale_colour_manual(values = phdhelpr::qual_colours[c(2, 6)], 
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
        axis.title = element_text(size = 14),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

# create chl a plot
chl <- ggplot(filter(chlorodata, day != "1" & day != "5"),
              aes(x = day,
                  y = conc,
                  colour = class,
                  shape = class)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_x_continuous(name = "Time (days)", 
                     breaks = c(1, 5, 8, 10, 12, 15), 
                     labels = c(1, 5, 8, 10, 12, 15)) +
  scale_y_continuous(name = expression(Chl~a~(pg~cell^-1)),
                     limits = c(0, 1)) +
  scale_colour_manual(values = gordon01::qual_colours[c(2, 6)], 
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
        axis.title = element_text(size = 14),
        axis.ticks = element_blank()
  )

# Extract legend, create grobs and set titles
g_tab <- ggplot_gtable(ggplot_build(dens))
leg <- which(sapply(g_tab$grobs, function(x) x$name) == "guide-box")
legend <- g_tab$grobs[[leg]]

plot_a <-
  gridExtra::arrangeGrob(dens +
                           theme(legend.position = "none"),
                         top = grid::textGrob("a",
                                              x = grid::unit(0.017, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 16)
                         ))

plot_b <-
  gridExtra::arrangeGrob(
    chl,
    top = grid::textGrob("b",
                         x = grid::unit(0.017, "npc"),
                         y = grid::unit(0.5, "npc"),
                         just = c("left", "top"),
                         gp = grid::gpar(fontsize = 16)
    ))

# print plot
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1),
                        legend,
                        nrow = 2,
                        heights = c(12, 1))

# save plot
grDevices::pdf("./figs/dens_chla_plot.pdf",
               width = 10,
               height = 6.5,
               useDingbats = FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1),
                        legend,
                        nrow = 2,
                        heights = c(12, 1))
grDevices::dev.off()
