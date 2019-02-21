# Script to run anovas on PSII yield

# load data
load("./data/metadata.rda")

# Bonferroni post hoc ANOVA
t0aov <- aov(FvFm ~ cont_treat, 
               data = metadata[metadata$day == "day 1", ])
summary(t0aov)

ttest_day1 <- pairwise.t.test(metadata$FvFm[metadata$day == "day 1"],
                              metadata$cont_treat[metadata$day == "day 1"],
                              p.adjust.method = "bonferroni")





pairwise.t.test(celldens$density[celldens$day == "day 5"],
                celldens$class[celldens$day == "day 5"],
                p.adjust.method = "bonferroni")

pairwise.t.test(celldens$density[celldens$day == "day 10"],
                celldens$class[celldens$day == "day 10"],
                p.adjust.method = "bonferroni")

pairwise.t.test(celldens$density[celldens$day == "day 14"],
                celldens$class[celldens$day == "day 14"],
                p.adjust.method = "bonferroni")
