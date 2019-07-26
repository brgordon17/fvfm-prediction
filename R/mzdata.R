# Code to preprocess mzdata_raw and construct mzdata
# Author: Benjamin R. Gordon
# Date: 2019-03-17

# Load libraries and data ------------------------------------------------------
library(tidyverse)
library(phdhelpr)
library(missForest)
library(Harman)
library(reshape2)

mzdata  <-  readr::read_csv("./data-raw/mzdata-raw.csv", na = "0")
load("./data/pamdata.rda")

# remove isotopes --------------------------------------------------------------
mzdata <- mzdata[-grep("[M+1]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+2]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+3]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+4]", mzdata$isotopes, fixed = TRUE),]

# Clean up ---------------------------------------------------------------------
mzdata <- dplyr::select(mzdata, -isotopes, -adduct, -pcgroup)
mzdata <- mzdata[-2:-20]

mz_names <- round(mzdata[, 1], 5)
mz_names <- paste("mz", mz_names$mz, sep = "_")
mz_names <- make.names(mz_names, unique = TRUE)

mzdata <- tibble::as_tibble(t(mzdata), rownames = "sample_id", 
                            .name_repair = "universal")
colnames(mzdata)[2:ncol(mzdata)] <- mz_names
mzdata <- mzdata[-1, ]
  
# Create metadata --------------------------------------------------------------
newpam <- 
  pamdata %>% 
  group_by(class, day) %>%
  filter(collect == TRUE) %>%
  mutate(mean_yield = mean(yield)) %>%
  arrange(day, class) %>%
  ungroup() %>%
  slice(-46:-48)

metadata <- tibble(sample_id = mzdata$sample_id,
                   class = factor(c(rep("PBQC", 18),
                                    rep("T0-C", 12),
                                    rep("T0-T", 12),
                                    rep("T1-C", 12),
                                    rep("T1-T", 9),
                                    rep("T2-C", 12),
                                    rep("T2-T", 12),
                                    rep("T3-C", 12),
                                    rep("T3-T", 12),
                                    rep("T4-C", 12),
                                    rep("T4-T", 12),
                                    rep("T5-C", 12),
                                    rep("T5-T", 12)
                                    )),
                   batch = factor(c(rep("batch 1", 1),
                                    rep("batch 2", 3),
                                    rep(c("batch 1", "batch 2"), 6),
                                    rep("batch 2", 2),
                                    rep("batch 2", 24),
                                    rep("batch 1", 21),
                                    rep("batch 2", 48),
                                    rep("batch 1", 24),
                                    rep("batch 2", 24)
                                    )),
                   day = factor(c(rep("PBQC", 18),
                                  rep("day 1", 12),
                                  rep("day 1", 12),
                                  rep("day 5", 12),
                                  rep("day 5", 9),
                                  rep("day 8", 12),
                                  rep("day 8", 12),
                                  rep("day 10", 12),
                                  rep("day 10", 12),
                                  rep("day 12", 12),
                                  rep("day 12", 12),
                                  rep("day 15", 12),
                                  rep("day 15", 12)),
                                levels = c("day 1", "day 5", "day 8", "day 10",
                                           "day 12", "day 15", "PBQC")
                                ),
                   cont_treat = factor(c(rep("PBQC", 18),
                                         rep("C", 12),
                                         rep("T", 12),
                                         rep("C", 12),
                                         rep("T", 9),
                                         rep("C", 12),
                                         rep("T", 12),
                                         rep("C", 12),
                                         rep("T", 12),
                                         rep("C", 12),
                                         rep("T", 12),
                                         rep("C", 12),
                                         rep("T", 12)), 
                                       levels = c("C", "T", "PBQC")
                                       ),
                   FvFm = c(rep(NA, 18), c(newpam$mean_yield))
                   )

mzdata <- bind_cols(metadata, mzdata[-1])
mzdata <- type_convert(mzdata)
rm(newpam)

# Impute noise and remove unreliable features ----------------------------------
round(mean(is.na(mzdata))*100, 2)
mzdata_filt <- phdhelpr::missing_values(mzdata[c(-1, -3:-6)],
                                        column.cutoff = 0.8,
                                        group.cutoff = 0.65,
                                        complete.matrix = FALSE,
                                        seed = 1978)
mzdata <- bind_cols(metadata,
                    mzdata_filt$output[-1])
round(mean(is.na(mzdata))*100, 2)

rm(mzdata_filt, mz_names)

# Impute remaining missing values ----------------------------------------------
doMC::registerDoMC()
set.seed(1978)
mzdata.imp <- missForest(as.data.frame(mzdata[-1:-6]),
                         ntree = 100,
                         verbose = TRUE,
                         maxiter = 8,
                         parallelize = "variables")

mzdata <- bind_cols(metadata, mzdata.imp$ximp)
round(mean(is.na(mzdata[-1:-6]))*100, 2)

# remove batch effects ---------------------------------------------------------
# Create across group relative log abundance plot
box_raw <- ggplot(data = reshape2::melt(filter(mzdata[-6], class != "PBQC")),
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
        axis.title.y = element_text(size = 12, colour = "grey65"),
        axis.title.x = element_text(size = 12, colour = "grey65"),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank())
box_raw

# Remove batch effects ---------------------------------------------------------
harmdata <- as.data.frame(t(mzdata[-1:-6]))
harmdata <- log(harmdata, 2)
colnames(harmdata) <- mzdata$sample_id

# Run Harman
harm <- harman(harmdata, 
               expt = mzdata$cont_treat, 
               batch = mzdata$batch,
               limit = 0.99)
summary(harm)

# Reconstruct the data
mzdata_cor <- tibble::as_tibble(t(reconstructData(harm)))
mzdata_cor <- 2^mzdata_cor
mzdata_cor <- dplyr::bind_cols(metadata, mzdata_cor)

# Explore variation ------------------------------------------------------------

# across group relative log abundance plot
box_cor <- ggplot(data = reshape2::melt(mzdata_cor[-6]),
                  aes(x = sample_id, 
                      y = log(value, 2) - median(log(value, 2)),
                      fill = batch)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  labs(x = "Samples", y = "Relative log Abundance") +
  scale_fill_manual(values = gordon01::qual_colours) +
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
pca_orig <- stats::prcomp(mzdata[-1:-6], 
                          scale = FALSE, 
                          center = TRUE)

exp_var_orig <- summary(pca_orig)$importance[2 ,]
scores_orig <- data.frame(metadata, pca_orig$x)
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
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

# PCA plot corrected (PBQCs excluded)
set.seed(1978)
pca_data <- filter(mzdata_cor, cont_treat != "PBQC")
pca_cor <- stats::prcomp(pca_data[-1:-6], scale = FALSE, center = TRUE)

exp_var_cor <- summary(pca_cor)$importance[2 ,]
scores_cor <- data.frame(filter(metadata, cont_treat != "PBQC"), pca_cor$x)
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
  scale_shape_manual(values = c(21:23)) +
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
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

# rebuild mzdata ---------------------------------------------------------------
mzdata <- mzdata_cor

# write data -----------------------------------------------------------------
save(mzdata, file = "./data/mzdata.rda", compress = "bzip2")
save(metadata, file = "./data/metadata.rda", compress = "bzip2")
saveRDS(box_raw, "./dev/box_raw_ggobj.rds")
saveRDS(harm, "./dev/harman_corr.rds")

## Data documentation ----------------------------------------------------------

#' Preprocessed LCMS data
#'
#' A dataset of positive ESI-LCMS mz variables for coral samples subjected to
#' experimental ocean warming and acidification at Heron Island in Feb 2011.
#'
#' @format A tibble with 101 rows and 1647 variables:
#' \describe{
#'   \item{sample_ids}{ID's for each sample}
#'   \item{class}{Treatment group. One of, control, eT (eleveated temperature),
#'   eCO2 (elevated CO2), eCO2eT (elevated CO2 and temperature), PBQC}
#'   \item{day}{Days of exposure to treatment}
#'   \item{tank}{Tank membership. Either L (left) or R (right)}
#'   \item{rep}{Sample replicate number}
#'   \item{class_day}{interaction variable between class and day variables}
#' }
#' @source Benjamin R. Gordon
"mzdata"
