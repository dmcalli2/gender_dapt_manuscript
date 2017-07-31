# 08b_process_complicated_plot
source ("scripts/01_read_data.R")

library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
library(stringr)

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

# Summarise treatment effect at trial level in men and women
forest_data <- bind_rows (trials, c(param = "inter", variable = "pooled", pooled))
forest_data <- forest_data %>%
  add_row(param = "men", variable = "pooled") %>%
  add_row(param = "women", variable = "pooled")

# Create table for beside plot
evnt_smry <- merge(main_wide [ , c(1, 6:13)],
                   bleed[, c(1, 6:9)],
                   by = "trial",
                   all.x = TRUE)

# Set order for figure
trial_ordr <- structure(list(trial_order = c("ASPS3", "BCHANCE", "CTRITON-TIMI 38", 
                                             "DTRILOGY ACS", "EPLATO", "FCURE", "GCOMMIT", "HCLARITY-TIMI 28"),
                             trial_label = c("SPS3", "CHANCE", "TRITON-TIMI 38", "TRILOGY ACS", "PLATO", "CURE", 
                                             "COMMIT", "CLARITY-TIMI 28")),
                        class = c("tbl_df", "tbl", "data.frame"),
                        row.names = c(NA, -8L), 
                        .Names = c("trial_order", "trial_label"))

row.names(evnt_smry) <- evnt_smry$trial
evnt_smry$trial <- NULL
evnt_smry <- evnt_smry[rev(trial_ordr$trial_label),]
# Add in totals row
mace <- as.matrix(evnt_smry[,5:8])
n <- as.matrix(evnt_smry[,1:4])
mace_per <- SmryPrcnt(mace, n)
mace_per <- as.data.frame(mace_per)
names(mace_per) <- c("women_i", "women_c", "men_i", "men_c")
mace_per$women <- paste(mace_per$women_i, mace_per$women_c, sep = "/")
mace_per$men <- paste(mace_per$men_i, mace_per$men_c, sep = "/")
mace_per <- mace_per[, c("men", "women")]
row.names(mace_per) <- row.names(evnt_smry)
# Save table so can read in final paper
saveRDS(mace_per, file = "data/Table results for forest plot.Rds")

## Use to draw plot
forest_table <- mace_per

############
# Set ID for ordering
forest_data$trial_id <- group_indices(forest_data, variable)
forest_data$trial_id[forest_data$variable == "pooled"] <- max(forest_data$trial_id) + 1
trial_labels <- forest_data %>% 
  select(trial_id, variable) %>%
  distinct() %>%
  arrange(-trial_id)

# Set into left and right panel
forest_data$panel <- ifelse(forest_data$param != "inter", "Gender-specific", "Interaction")

## Relabel for printing
forest_data$trial_id <- factor(forest_data$trial_id,
                               levels = trial_labels$trial_id,
                               labels = trial_labels$variable)

forest_data$param <- factor(forest_data$param,
                            levels = c("men", "women", "inter"),
                            labels = c("Men", "Women", "Interaction"))

forest_data$pooled <- ifelse(forest_data$variable == "pooled", "pooled", "not")

forest_table <- forest_table[rev(levels(forest_data$trial_id)),]
forest_table$trial_id <- rev(levels(forest_data$trial_id))
forest_table$trial_id <- factor(forest_table$trial_id, levels = levels(forest_data$trial_id))
forest_table <- as_tibble(forest_table) %>%
  select(trial_id, men, women)
forest_table$men <- paste("Men", forest_table$men)
forest_table$women <- paste("Women", forest_table$women)
forest_table$all_text <- paste(forest_table$trial_id, forest_table$men, forest_table$women, sep = "\n")
forest_table[nrow(forest_table), -1] <- rep(NA, 3)
forest_data <- inner_join(forest_data, select(forest_table, trial_id, all_text), by = "trial_id")
forest_data$all_text[nrow(forest_data) - 3] <- paste(c("", "Interaction rate ratio",
                                                       "(shared model)"), collapse = "\n")
# Draw plots
forest_plot <- ggplot(forest_data, aes(x = as.numeric(trial_id), colour = param, shape = pooled,
                                       y = q50, ymin = q2.5, ymax = q97.5)) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(position = position_dodge(width = 0.5)) +
  facet_grid(~ panel) +
  coord_flip() +
  scale_y_continuous("Rate ratio (95% CI)") +
  scale_x_continuous(NULL, breaks = 9:1, labels = unique(na.omit(forest_data$all_text))) +
  scale_color_discrete(NULL, guide = FALSE) +
  scale_shape(NULL, guide = FALSE) +
  geom_hline(mapping = aes(yintercept = 1), linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0),
        axis.text.y = element_text(angle = 0, hjust = 0))
forest_plot

ggsave("figures/forest_plot_with_table.tiff", plot = forest_plot, units = "in",
       height = 6, width = 12, compression = "lzw")

saveRDS(forest_plot, file = "data/Final forest plot stratified and interaction.RDS")


