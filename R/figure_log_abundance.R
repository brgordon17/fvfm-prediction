# Code to construct the relative log abundance plot
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)

# Load data 
load("./data/mzdata.rda")
box_raw <- readRDS("./dev/box_raw_ggobj.rds")

# RLA of corrected data
box_cor <- ggplot(data = reshape2::melt(filter(mzdata[-6], class != "PBQC")),
                  aes(x = sample_id, 
                      y = log(value, 2) - median(log(value, 2)),
                      fill = batch)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  scale_fill_manual(values = phdhelpr::qual_colours) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Relative log Abundance") +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 11, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")
box_cor

# Extract Legend 
g_tab <- ggplot_gtable(ggplot_build(box_cor))
leg <- which(sapply(g_tab$grobs, function(x) x$name) == "guide-box")
legend <- g_tab$grobs[[leg]]
rm(g_tab, leg)

# create grobs
plot_a <-
  gridExtra::arrangeGrob(box_raw +
                           theme(legend.position = "none",
                                 axis.title.y = element_blank()),
                         top = grid::textGrob("a", 
                                              x = grid::unit(0.017, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 16)
                         ))

plot_b <-
  gridExtra::arrangeGrob(box_cor +
                           theme(legend.position = "none"),
                         top = grid::textGrob("b",
                                              x = grid::unit(0.017, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 16)
                         ))

xgrob <- gridExtra::arrangeGrob(grid::textGrob("Sample",
                                               x = grid::unit(0.6, "npc"),
                                               y = grid::unit(0.5, "npc"),
                                               just = c("centre"),
                                               gp = grid::gpar(fontsize = 14)),
                                right = legend)


ygrob <- grid::textGrob("Relative log Abundance",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        just = c("centre"),
                        rot = 90,
                        gp = grid::gpar(fontsize = 14)
)


# print plot
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 2,
                                               left = ygrob),
                        xgrob,
                        nrow = 2,
                        heights = c(12, 1))

# Save plot
grDevices::pdf("./figs/RLA_plot.pdf",
               width = 10,
               height = 6.5,
               useDingbats = FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 2,
                                               left = ygrob),
                        xgrob,
                        nrow = 2,
                        heights = c(12, 1))
grDevices::dev.off()
