# Script to perform and plot PCA analysis. Conducted on 27SEP18 by B. Gordon
#
# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Load data --------------------------------------------------------------------
load("./data/mzdata.rda")

# PCA --------------------------------------------------------------------------
set.seed(1978)
pca <- stats::prcomp(mzdata[-1:-7], scale = FALSE, center = TRUE)
  
# Define variables for PCA plot ------------------------------------------------
exp_var <- summary(pca)$importance[2 ,]
scores <- data.frame(mzdata[, 2:6], pca$x)
x_lab <- paste("PC1", " (", round(exp_var[1] * 100, 2), "%)", sep =  "")
y_lab <- paste("PC2", " (", round(exp_var[2] * 100, 2), "%)", sep =  "")
custom_colours <- gordon01::qual_colours[c(1, 2, 7)]
  
  
# create class plot ------------------------------------------------------------
classplot <-
  ggplot(data = scores,
          aes(x = PC1,
              y = PC2,
              color = class,
              fill = class,
              shape = class)) +
   geom_point(size = 2.5,
              stroke = 0.7,
              position = position_jitter(width = 0.01 * diff(range(scores$PC1)),
                                         height = 0.01 * diff(range(scores$PC2))
                                         )) +
  labs(x = x_lab, y = y_lab) +
  scale_shape_manual(name = NULL,
                     values = c(21:25)) +
  scale_color_manual(name = NULL,
                     values = grDevices::adjustcolor(custom_colours,
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = NULL,
                    values = grDevices::adjustcolor(custom_colours,
                                                    alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey50"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 10))
classplot
  
# create batch plot --------------------------------------------------------------
custom_colours <- gordon01::qual_colours[c(1, 2)]
batchplot <-
  ggplot(data = scores,
         aes(x = PC1,
             y = PC2,
             color = batch,
             fill = batch,
             shape = batch)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores$PC1)),
                                        height = 0.01 * diff(range(scores$PC2))
                                        )) +
  labs(x = x_lab, y = y_lab) +
  scale_shape_manual(values = c(21, 22)) +
  scale_color_manual(values = grDevices::adjustcolor(custom_colours,
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(custom_colours,
                                                    alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey50"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 10))
batchplot

# Run Harman package to remove batch effects and revisit.
