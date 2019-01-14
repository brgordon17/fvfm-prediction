# Script to remove batch effects from mzdata using the Harman package. Conducted
# by B. Gordon on 14Jan18

# Load Libraries
library(tidyverse)
library(Harman)
library(reshape2)

# Explore variation ------------------------------------------------------------
# Load data
load("./data/mzdata.rda")
load("./data/phenodata.rda")

# split the data into its individual experiments
mzdata2011 <- 
  mzdata %>% 
  filter(experiment == "2011")

PBQCs <-
  mzdata %>% 
  filter(cont_treat == "PBQC")

mzdata2013 <- 
  mzdata %>% 
  filter(experiment == "2013" & cont_treat != "PBQC")

mzdata2014 <- 
  mzdata %>% 
  filter(experiment == "2014")

# reshape the data
long2013 <- reshape2::melt(mzdata2013)

# Create across group relative log abundance plot
custom_colours <- gordon01::qual_colours
box_orig <- ggplot(data = long2013,
                   aes(x = sample_id, 
                       y = log(value, 2) - median(log(value, 2)),
                       fill = batch)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  labs(x = "Samples", y = "Relative log Abundance") +
  scale_fill_manual(values = custom_colours) +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey70"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10))
box_orig

# Remove batch effects ---------------------------------------------------------
# Using log data

# define variables
phenod <- mzdata2013[1:7] 
harmdata <- as.data.frame(t(mzdata2013[-1:-7]))
harmdata <- log(harmdata, 2)
colnames(harmdata) <- phenod$sample_id

# Run Harman
mzharm <- harman(harmdata, 
                 expt = phenod$cont_treat, 
                 batch = phenod$batch,
                 limit = 0.95)
#summary(mzharm)

# Reconstruct the data
mzdata_cor <- tibble::as_tibble(t(reconstructData(mzharm)))
mzdata_cor <- dplyr::bind_cols(phenod, mzdata_cor)

# Explore variation ------------------------------------------------------------
# reshape the data
long2013_cor <- reshape2::melt(mzdata_cor)

# Create across group relative log abundance plot
custom_colours <- gordon01::qual_colours
box_cor <- ggplot(data = long2013_cor,
                  aes(x = sample_id, 
                      y = value - median(value),
                      fill = batch)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  labs(x = "Samples", y = "Relative log Abundance") +
  scale_fill_manual(values = custom_colours) +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey70"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10))
box_cor

# PCA plot original
set.seed(1978)
pca_orig <- stats::prcomp(mzdata2013[-1:-7], 
                          scale = FALSE, 
                          center = TRUE)

# Define variables
exp_var_orig <- summary(pca_orig)$importance[2 ,]
scores_orig <- data.frame(mzdata2013[, 2:7], pca_orig$x)
x_lab_orig <- paste("PC1", " (", round(exp_var_orig[1] * 100, 2), "%)", sep =  "")
y_lab_orig <- paste("PC2", " (", round(exp_var_orig[2] * 100, 2), "%)", sep =  "")

pca_orig <-
  ggplot(data = scores_orig,
         aes(x = PC1,
             y = PC2,
             color = batch,
             fill = batch,
             shape = batch)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores_orig$PC1)),
                                        height = 0.01 * diff(range(scores_orig$PC2))
             )) +
  labs(x = x_lab_orig, y = y_lab_orig) +
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
pca_orig

# PCA plot corrected
set.seed(1978)
pca_cor <- stats::prcomp(mzdata_cor[-1:-7], scale = FALSE, center = TRUE)

# Define variables
exp_var_cor <- summary(pca_cor)$importance[2 ,]
scores_cor <- data.frame(mzdata2013[, 2:7], pca_cor$x)
x_lab_cor <- paste("PC1", " (", round(exp_var_cor[1] * 100, 2), "%)", sep =  "")
y_lab_cor <- paste("PC2", " (", round(exp_var_cor[2] * 100, 2), "%)", sep =  "")

# plot of control vs treatment
pca_cor_CvsT <-
  ggplot(data = scores_cor,
         aes(x = PC1,
             y = PC2,
             color = cont_treat,
             fill = cont_treat,
             shape = cont_treat)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores_cor$PC1)),
                                        height = 0.01 * diff(range(scores_cor$PC2))
             )) +
  labs(x = x_lab_cor, y = y_lab_cor) +
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
pca_cor_CvsT

# plot of batch
pca_cor_batch <-
  ggplot(data = scores_cor,
         aes(x = PC1,
             y = PC2,
             color = batch,
             fill = batch,
             shape = batch)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores_cor$PC1)),
                                        height = 0.01 * diff(range(scores_cor$PC2))
             )) +
  labs(x = x_lab_cor, y = y_lab_cor) +
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
pca_cor_batch

# rebuild mzdata and save ------------------------------------------------------
mzdata_final <- bind_rows(mzdata2011, PBQCs, mzdata_cor, mzdata2014)
save(mzdata_final, file = "./data/mzdata.rda", compress = "bzip2")
