# Script that contains code for all figures used in chapter 5 (bleaching) of
# thesis.
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

# PSII FvFm plot ---------------------------------------------------------------
library(tidyverse)

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

# Cell density and chlorophyll plot --------------------------------------------
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
        legend.text = element_text(size = 12),
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

# composite figure of relative log abundance -----------------------------------
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
  scale_fill_manual(values = gordon01::qual_colours) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Relative log Abundance") +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "grey65"),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
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
                           theme(legend.position = "none"),
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
                                               gp = grid::gpar(fontsize = 12,
                                                               col = "grey65")),
                                right = legend)

# print plot
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1),
                        xgrob,
                        nrow = 2,
                        heights = c(12, 1))

# Save plot
grDevices::pdf("./figs/RLA_plot.pdf",
               width = 10,
               height = 3.5,
               useDingbats = FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1),
                        xgrob,
                        nrow = 2,
                        heights = c(12, 1))
grDevices::dev.off()

# PCA plot ----------------------------------
library(tidyverse)
library(GGally)

# Load data 
load("./data/mzdata.rda")
mzdata <- filter(mzdata, class != "PBQC")

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
                         shape = cont_treat),
                     columns = c(6:10),
                     axisLabels = "none",
                     showStrips = FALSE,
                     legend = 2,
                     upper = list(continuous = wrap("smooth", method = "lm", 
                                                    se = FALSE)),
                     lower = list(continuous = wrap("points")),
                     diag = list(continuous = wrap("densityDiag"))
                     ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11))

# Change colours etc (must be done by extrating each plot)
for(i in 1:pca_pairs$nrow) {
  for(j in 1:pca_pairs$ncol){
    pca_pairs[i,j] <- pca_pairs[i,j] +
      scale_shape_manual(name = NULL,
                         values = c(21, 24),
                         label = c("control", "heated")) +
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

# Important variables plots ----------------------------------------------------

# Load models 
mzrf_fvfm <- readRDS("./dev/mzrf_model_fvfm.rds")
mzrf_cvst <- readRDS("./dev/mzrf_model_cvst.rds")

# Identify important variables 
fvfm_impvars <- varImp(mzrf_fvfm, scale = TRUE)
fvfm_impvars <- fvfm_impvars$importance
fvfm_impvars <- cbind(vip = apply(fvfm_impvars, 1, max), fvfm_impvars)
fvfm_impvars <- fvfm_impvars[order(-fvfm_impvars$vip), ,drop = FALSE]

cvst_impvars <- varImp(mzrf_cvst, scale = TRUE)
cvst_impvars <- cvst_impvars$importance
cvst_impvars <- cbind(vip = apply(cvst_impvars, 1, max), cvst_impvars)
cvst_impvars <- cvst_impvars[order(-cvst_impvars$vip), ,drop = FALSE]

# clean up dfs and convert to tibbles
rownames(fvfm_impvars) <- gsub("mz_", "", rownames(fvfm_impvars))
fvfm_impvars <- cbind(mz = factor(rownames(fvfm_impvars),
                                   levels = rev(rownames(fvfm_impvars))
                                  ), fvfm_impvars)
fvfm_impvars <- as_tibble(fvfm_impvars[1:20, ], rownames = NULL)

rownames(cvst_impvars) <- gsub("mz_", "", rownames(cvst_impvars))
cvst_impvars <- cbind(mz = factor(rownames(cvst_impvars),
                                  levels = rev(rownames(cvst_impvars))
                                  ), cvst_impvars)
cvst_impvars <- tibble::as_tibble(cvst_impvars[1:20, ], rownames = NULL)

# create plots 
fvfm_plot <- ggplot(fvfm_impvars,
                     aes(x = vip, y = mz)) +
  geom_point(shape = 16,
             colour = gordon01::seq_colours[4],
             size = 3) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10,
                                   colour = "grey70"),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.position = "none")

cvst_plot <- ggplot(cvst_impvars,
                    aes(x = vip, y = mz)) +
  geom_point(shape = 16,
             colour = gordon01::seq_colours[4],
             size = 3) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10,
                                   colour = "grey70"),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.position = "none")

# Create composite plot and axis text grobs
fvfm_plot <-
  gridExtra::arrangeGrob(fvfm_plot,
                         top = grid::textGrob("Quantum Yield",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              gp = grid::gpar(fontsize = 12)))

cvst_plot <-
  gridExtra::arrangeGrob(cvst_plot,
                         top = grid::textGrob("Class",
                         x = grid::unit(0.5, "npc"),
                         y = grid::unit(0.5, "npc"),
                         gp = grid::gpar(fontsize = 12)))

xgrob <- grid::textGrob("Variable Importance",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        just = c("centre"),
                        gp = grid::gpar(fontsize = 12))

ygrob <- grid::textGrob("m/z",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        just = c("centre"),
                        rot = 90,
                        gp = grid::gpar(fontsize = 12))

# save plot 
grDevices::pdf("./figs/vip_plot.pdf",
                width = 10,
                height = 5,
                useDingbats = FALSE)
gridExtra::grid.arrange(fvfm_plot,
                        cvst_plot,
                        nrow = 1,
                        left = ygrob,
                        bottom = xgrob)
grDevices::dev.off()
