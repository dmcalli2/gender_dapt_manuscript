# 08b_process_complicated_plot
source ("scripts/01_read_data.R")

library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
library(stringr)
library(ggthemes)

## Data and summaries
load (file = "model_summaries/stratified_each_trial.Rdata")
load("model_summaries/jags_samples_main.Rdata")


# Summarise result
ratio_upper_lower <- format(round(inter_q[3,] / inter_q[2,],2),nsmall = 2)
inter_q <- as.data.frame(t(format(round(inter_q,2), nsmall = 2)), stringsAsFactors = FALSE)
myrownames <- c("stratified", "Pooled", "fixed", "Pooled_tx_strat", "Fixed_tx_strat","lau", "lau_sa")
row.names(inter_q) <- myrownames
pooled <- inter_q["Pooled" , ]
names(pooled) <- paste0("q", c(50, 2.5, 97.5))
pooled[] <- lapply(pooled, as.double)


## Plot Figure 2A, gender-specific effects
# Create combined file in inkscape
table_order <- main_wide

table_order$comparison[table_order$indication == "Stroke"] <- "Stroke"
table_order$comparison[table_order$indication == "ACS" & table_order$comparison == "clopidogrel_placebo"] <- 
  "ACS, clopidogrel"
table_order$comparison[table_order$indication == "ACS" & table_order$comparison == "prasugrel_clopidogrel"] <- 
  "ACS, novel P2Y12"

table_order <- table_order %>% 
  mutate(comparison_f = factor(comparison, levels = c("ACS, novel P2Y12", "ACS, clopidogrel",
                                                                  "Stroke"),
                               labels = c("1 ACS, novel P2Y12", "2 ACS, clopidogrel",
                                                                  "3 Stroke"),
                                ordered = TRUE)) %>% 
  arrange(comparison_f, trial) %>% 
  mutate(my_order = 1:nrow(.))
table_order$my_order [table_order$my_order %in% 7:9] <- table_order$my_order [table_order$my_order %in% 7:9]+ 2
table_order$my_order [table_order$my_order %in% 4:6] <- table_order$my_order [table_order$my_order %in% 4:6]+ 1


men_women <- trials %>% 
  mutate(trial = as.character(variable),
         Sex = factor(param, levels = c("men", "women"), labels = c("Men", "Women"))) %>% 
  filter(param %in% c("men", "women")) %>% 
  inner_join(table_order) %>% 
  ggplot(aes(x = my_order, y = q50, ymin = q2.5, ymax = q97.5,
             colour = Sex, shape = Sex,  group = interaction(param, comparison))) +
  geom_point(position = position_dodge(width = 0.75)) + 
  geom_errorbar(position = position_dodge(width = 0.75)) +
  scale_x_continuous(name = NULL, breaks = NULL) +
  scale_y_continuous(name = "Rate ratio") +
  coord_flip() +
  geom_hline(mapping = aes(yintercept = 1), linetype = "dashed") +
  scale_shape(name = NULL) +
    scale_color_colorblind(name = NULL) +
  theme_classic()
win.metafile  (filename = "figures/figure2_colourblindsafe.emf")
men_women
dev.off()

## Table to go with plot
trials_wide <- trials %>% 
  mutate(trial = as.character(variable)) %>% 
  mutate(RR = paste0(q50, " (", q2.5, " to ", q97.5, ")")) %>% 
  dplyr::select(trial, param, RR) %>% 
  spread(key = param, value = RR, sep = "_")
men_women <- trials_wide %>% 
  dplyr::select(trial, param_men, param_women)


men_women_table <- table_order %>% 
  inner_join(men_women) %>% 
  dplyr::select(my_order, trial, comparison,
                rwt, nwt,
                rwc, nwc,
                param_women,
                rmt, nmt,
                rmc, nmc,
                param_men) 

write_csv(men_women_table, path = "figures/table_figure2a.csv")
tapply(table_order$my_order, table_order$comparison, range)

## Summarise interaction effect
inters_table <- trials_wide %>% 
  inner_join(table_order) %>% 
    arrange(-my_order) %>% 
  dplyr::select(my_order, trial, comparison, RR = param_inter)
write_csv(inters_table, path = "figures/table_figure2b.csv")

inters <- trials %>% 
  filter(param == "inter") %>% 
  bind_rows(pooled) %>% 
  mutate(shape_type = if_else(is.na(param), "pooled", "trial"),
         trial = as.character(variable),
         my_size = 1/(log(q97.5/q50)/1.96)^2,
         my_size = if_else(shape_type == "pooled", my_size *2, my_size)) # in R diamond shape sits inside square shape

inters_plot <- inters %>% 
  left_join(table_order) %>% 
  mutate(my_order = if_else(is.na(my_order), 1, my_order + 2)) %>% 
  ggplot(aes(x = my_order, y = q50, ymin = q2.5, ymax = q97.5, shape = shape_type)) +
  geom_point(mapping = aes(size = my_size)) + 
  geom_errorbar() +
  # scale_x_continuous(name = NULL, breaks = NULL) +
  scale_y_continuous(name = "Rate ratio") +
  coord_flip() +
  geom_hline(mapping = aes(yintercept = 1), linetype = "dashed") +
  scale_shape_manual(name = NULL, guide = FALSE, values = c(18, 15)) +
  scale_size(name = NULL, guide = FALSE) +
  theme_classic()
inters_plot

win.metafile  (filename = "figures/figure2b_colourblindsafe.emf")
inters_plot
dev.off()
