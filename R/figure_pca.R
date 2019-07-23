# Code to construct the PCA plot
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)
library(GGally)

# Load data 
load("./data/mzdata.rda")
mzdata <- filter(mzdata, class != "PBQC" & cont_treat != "C")

# PCA of batch corrected data (PBQCs excluded)
set.seed(1978)
pca <- stats::prcomp(mzdata[-1:-6], 
                     scale = FALSE, 
                     center = TRUE)

exp_var <- summary(pca)$importance[2 ,]
pc_names <- paste(colnames(pca$x), " (", round(exp_var * 100, 1), "%)", sep =  "")
scores <- data.frame(mzdata[, 1:5], pca$x)
colnames(scores)[6:ncol(scores)] <- pc_names
custom_colours <- gordon01::qual_colours

# create pairs plot of first 5 PCs
pca_pairs <- ggpairs(data = scores,
                     aes(colour = day,
                         fill = day,
                         shape = day),
                     columns = c(6:10),
                     axisLabels = "none",
                     showStrips = FALSE,
                     legend = c(2,1),
                     upper = list(continuous = wrap("smooth", method = "lm",
                                                    se = FALSE,
                                                    size = 1.5)),
                     lower = list(continuous = wrap("points", size = 1.5)),
                     diag = list(continuous = wrap("densityDiag"))
) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11))

# Change colours etc (must be done by extrating each plot)
for(i in 1:pca_pairs$nrow) {
  for(j in 1:pca_pairs$ncol){
    pca_pairs[i,j] <- pca_pairs[i,j] +
      scale_shape_manual(name = NULL,
                         values = c(4, 21:25)) +
      scale_color_manual(name = NULL,
                         values = grDevices::adjustcolor(custom_colours,
                                                         alpha.f = 0.9)) +
      scale_fill_manual(name = NULL,
                        values = grDevices::adjustcolor(custom_colours,
                                                        alpha.f = 0.5)) +
      theme(legend.title = element_blank(),
            legend.key = element_rect(fill = "transparent", colour = NA),
            legend.text = element_text(size = 11))
  }
}
pca_pairs

ggsave("./figs/pca_plot.pdf",
       width = 10,
       height = 6.5,
       units = "in")

# NOT USED
# + stat_ellipse(data = filter(scores, cont_treat == "T" & day == "day 12" |
#                              cont_treat == "T" & day == "day 15"),
#              show.legend = FALSE,
#              type = "norm")
