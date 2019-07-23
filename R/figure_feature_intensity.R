# Code to construct the feature intensity plot
# Author: Benjamin R. Gordon
# Date: 21-FEB-2019

library(tidyverse)
library(reshape2)

# Load and prep data
load("./data/mzdata.rda")

sum_mzdata <- 
  mzdata %>% 
  filter(class != "PBQC") %>%
  select(1:5, mz_249.18855, mz_170.04777, mz_985.69447) %>%
  melt(.) %>%
  group_by(variable, day, cont_treat) %>%
  summarise(mean = mean(value/10^4), 
            sd = sd(value/10^4), 
            se = sd(value/10^4)/sqrt(n())
  )
sum_mzdata$variable <- str_replace_all(sum_mzdata$variable, "mz_", "")

names <- c(
  "249.18855" = "dihomomontiporyne H",
  "170.04777" = "alanine betaine",
  "985.69447" = "lyso-PAF C16"
)

# create plot
ggplot(sum_mzdata, aes(x = day,
                       y = mean,
                       colour = cont_treat,
                       shape = cont_treat,
                       group = cont_treat)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = .1,
                show.legend = FALSE) +
  geom_path(size = 0.8, show.legend = FALSE) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.1, linetype = 3) +
  facet_wrap(vars(variable), labeller = as_labeller(names)) +
  scale_x_discrete(name = "Time (days)", 
                   labels = c(1, 5, 8, 10, 12, 15)) +
  scale_y_continuous(name = Mean~Intensity~x~10^4) +
  scale_colour_manual(values = gordon01::qual_colours[c(2,1)], 
                      name = NULL,
                      labels = c("control", "heated")) +
  scale_shape_manual(values = c(16, 17),
                     name = NULL,
                     labels = c("control", "heated")) +
  theme(axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        strip.background = element_blank()
  )

ggsave("./figs/feature_intensities.pdf")
