# 08_process_plots_tables
source ("scripts/01_read_data.R")

library(ggplot2)
library(tidyverse)
library(stringr)
library(rjags)


## Load in model results (inter_q)
load(file = "model_summaries/jags_samples_main.Rdata")

## Prepare tables from Scotland data summaries
scotland_res <- vector(length = 0L) # vector to contain object names for ES results for Scotland

for (bleed_effect in c("0", seq(20, 50, 10))) {
  for (aggrgt in c("", "avd_non_cllps")){
    for (gender_sa in c("", "null_mace")){
      label_file <- paste(bleed_effect, aggrgt, gender_sa, sep = "_")
      print(label_file)
      load(file = paste0("synth_smry", "/", "summary_of_rate_ratio", label_file, ".Rdata" ))
      load(file = paste0("synth_smry", "/", "summary_of_gender_age_spec_arr_", label_file, ".Rdata" ))
      load(file = paste0("synth_smry", "/", "summary_of_arr", label_file, ".Rdata" ))

      arr_smry$res <- paste0(arr_smry$est, " (95% CI ", arr_smry$lci, " to ", arr_smry$uci, ")")
      arr_smry$outcome <- factor (arr_smry$outcome, levels = c("cv", "bleed", "other", "all"),
                                  labels = c("Cardiovascular", "Bleeding", "Other", "Overall"))
      arr_smry$gender <- factor(arr_smry$gender, levels = c("men", "women", "both"),
                                c("Men", "Women", "Both"))
      arr_smry <- reshape2::acast(arr_smry, gender ~ outcome, value.var = "res")
      nnt_smry$nnt <- paste0(nnt_smry$est, " (95% CI ", nnt_smry$lci, " to ", nnt_smry$uci, ")")
      nnt_smry$outcome <- factor (nnt_smry$outcome, levels = c("cv", "bleed", "other", "all"),
                                  labels = c("Cardiovascular", "Bleeding", "Other", "Overall"))
      nnt_smry$gender <- factor(nnt_smry$gender, levels = c("men", "women", "both"),
                                c("Men", "Women", "Both"))
      nnt_smry <- reshape2::acast(nnt_smry, gender ~ outcome, value.var = "nnt")
      
      names(gender_age_spec) <- c("Age", "Sex", "Bleeding", "Cardiovascular", "Other")
   
      for (k in c("rrs", "arr_smry", "nnt_smry", "gender_age_spec")){
        assign(paste0(k, "_", bleed_effect, aggrgt, gender_sa), get(k))
        scotland_res <- c(scotland_res, paste0(k, "_", bleed_effect, aggrgt, gender_sa))
      }
}
}
}
save(list = scotland_res, file = "data/Results_scotland_effects.Rdata")

## Create tables for supplement
arr_nnt <- new.env()
load(file = "data/Results_scotland_effects.Rdata", envir = arr_nnt)
arr_nnt <- as.list(arr_nnt)

arr_nnt <- arr_nnt[grepl("arr", names(arr_nnt))]
arr_nnt <- arr_nnt[grepl("avd_non_cllps", names(arr_nnt))]

men_overall <-  arr_nnt[[1]] ["Men", "Overall"]

arr_nnt <- lapply(arr_nnt, function (x) x["Women", "Overall"])
arr_nnt <- do.call(rbind, arr_nnt)
arr_nnt <- arr_nnt[sort(rownames(arr_nnt)),]
arr_nnt <- as.data.frame(arr_nnt, as.is = TRUE)
arr_nnt$bleed <- substr(rownames(arr_nnt),10,11)
arr_nnt$bleed <- gsub("a", "", arr_nnt$bleed)
arr_nnt$bleed <- format(1+ as.numeric(arr_nnt$bleed) * 0.01, digits = 2, nsmall = 2)
arr_nnt$mace <- grepl("mace", rownames(arr_nnt))
arr_nnt$mace <- as.character(factor(arr_nnt$mace, levels = c(TRUE, FALSE), labels = c("MACE, 1.00", "MACE, 1.09")))
arr_nnt <- reshape2::dcast(arr_nnt, bleed ~ mace, value.var = "arr_nnt")
row.names(arr_nnt) <- paste0("Bleeding, ", arr_nnt$bleed)
arr_nnt$bleed <- NULL

save(arr_nnt, men_overall, file = "data/Summary arr sensitivity analyses interactions.Rdata")

## Plot of baseline risk, absolute risk reduction in women and men for main paper and supplement
## Baseline data estiamted for each stratum
br <- read.csv("data/baseline_risk_actual_data.csv", as.is = TRUE)
br <- br [, names(br) != "surv"]

br <- gather(br, key = otype, value = value, - age, -gender, -n)
PropEst <- function (value, n) 100* qbeta(c(0.025, 0.5, 0.975), 0.5 + value, 0.5 +n - value)
br$key <- formatC(1:nrow(br),flag=0,width=2)
mylist <- by(data = br, br$key, function(x) PropEst(x$value, x$n)) 
mylist <- do.call(rbind, mylist)
colnames(mylist) <- c("lci", "est", "uci")
br <- cbind(br, mylist)
br$age <- factor(br$age, levels = rev(unique(br$age)))
br$wrap <- "Baseline risk"

# Baseline risk data from model
load("synth_smry/baseline risks with data.Rdata")
bline_data <- y_chk 
bline_data$age <- paste0(bline_data$age, "0-", bline_data$age, "9")
bline_data$age <- factor(bline_data$age, levels = unique((bline_data$age)))
bline_data$gender <- factor(bline_data$gender, 0:1, labels = c("Women", "Men"))
bline_data$wrap <- "Baseline risk"
bline_data[ , c("est", "lci", "uci")] <- lapply(bline_data[, c("est", "lci", "uci")], function(x) x/ bline_data$n)

bline <- ggplot(filter(bline_data, otype == "cv"),
                aes(x = age, colour = gender, y = est, ymin = lci, ymax = uci)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
#  coord_flip() +
  scale_y_continuous("Risk (%) (95% CI)", labels = function(x) format(100*x, digits = 1, nsmall = 1)) +
  scale_x_discrete("Age (years)") +
  facet_wrap(~wrap)

## ARR
MakeArrData <- function (analysis_name){
  # Uses analysis name to extract ARR and label with analysis name
  load(file = paste0("synth_smry", "/summary_of_gender_age_spec_arr_", analysis_name, ".Rdata"))
  gender_age_spec <- gather(gender_age_spec, key = "otype", value = "res", -age, -gender)
  gender_age_spec <- inner_join (gender_age_spec, gender_age_spec_est_lci_uci, by = "res")
  names(gender_age_spec) [names(gender_age_spec) %in% c("2.5%", "50%", "97.5%")] <- 
    c("lci", "est", "uci")
  gender_age_spec[, c("lci", "est", "uci")] <- 
    lapply(gender_age_spec[, c("lci", "est", "uci")], as.double)
  gender_age_spec$age <- factor(gender_age_spec$age, levels = (unique(br$age)))
  gender_age_spec$analysis <- analysis_name
  gender_age_spec
}

# Create single long dataset for plots of ARR and baseline risk
# Identify all age, sex and gender-specific results, place in a list and then convert to a single big dataframe
arr_data <- list.files("synth_smry/")
arr_data <- arr_data[str_detect(arr_data, "summary_of_gender_age_spec_arr_")] %>%
  str_replace(fixed("summary_of_gender_age_spec_arr_"), "") %>%
  str_replace(fixed(".rdata"), "")
arr_data_names <- arr_data
arr_data <- as.list(arr_data)
names(arr_data) <- arr_data_names
rm(arr_data_names)

arr_data <- lapply(arr_data, MakeArrData)
arr_data <- do.call(rbind, arr_data)

# Separate the filename into each element 
arr_data$analysis <- str_replace(arr_data$analysis, "avd_non_cllps", "avdnoncllps")
arr_data$analysis <- str_replace(arr_data$analysis, "null_mace", "nullmace")

arr_data <- arr_data %>% separate(analysis, into = c("bleeding_rr", "non_collapse", "mace_rr"), sep = "_")
arr_data$bleeding_rr <- paste0((100+ as.integer(arr_data$bleeding_rr))/100, "-fold")
arr_data <- arr_data %>% filter(non_collapse != "avdnoncllps") %>%
  select(-non_collapse)
# Make string for effect estimate for interaction
pooled_string <-paste0(inter_q["50%", "re"], " (95%CI", inter_q["2.5%", "re"], " to ", inter_q["97.5%", "re"], ")")
arr_data$mace_rr <- ifelse(arr_data$mace_rr == "nullmace", "No interaction for MACE", 
                           paste0("Interaction ", pooled_string, " for MACE") )
arr_plot <- ggplot(filter(arr_data, otype == "MACE"), 
                   aes(x = age, colour = gender, y = est, ymin = lci, ymax = uci)) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid(bleeding_rr ~ mace_rr) +
#  coord_flip() +
  scale_y_continuous("Absolute risk reduction (%) (95% CI)") +
  scale_color_discrete("")

# Total plot
total_plot <- arr_plot %+% filter(arr_data, otype == "Total")

# Flip ARR and draw bleeding and competing risks plots
flip_arr_data <- arr_data
flip_arr_data [ , c("est", "uci", "lci")] <- lapply(arr_data [ , c("est", "lci", "uci")], function (x) -x)
bleed_plot   <- arr_plot %+% filter(flip_arr_data, otype == "Bleeding") + 
  scale_y_continuous("Absolute risk increase, bleeding death (%) (95% CI)")
compete_plot <- arr_plot %+% filter(flip_arr_data, otype == "Competing") + scale_y_continuous("Absolute risk increase, competing deaths (%) (95% CI)")
total_plot <- arr_plot %+% filter(arr_data, otype == "Total") + scale_y_continuous("Absolute risk increase, all deaths (%) (95% CI)")
  
# Create ARR plot comparing MACE interaction and no MACE interaction
mace_no_mace <- arr_plot %+% filter(arr_data, otype == "MACE", bleeding_rr == "1-fold") +
  facet_wrap(~mace_rr)

arr_single_data <-  arr_data %>% 
  filter(otype == "MACE", bleeding_rr == "1-fold", 
         mace_rr != "No interaction for MACE") %>% 
  mutate(treatment_effect = "Treatment effect")

arr_single <- arr_plot %+% 
  arr_single_data +
  facet_wrap(~treatment_effect) +
  scale_x_discrete("Age (years)") 


# Save age and sex specific ARR plots
save(arr_plot, bleed_plot, compete_plot, total_plot, mace_no_mace, bline, arr_single,
     file = "data/Bline, ARR plots.Rdata")

saveRDS(arr_data, file = "data/All_arr_data.Rds")


## Estimate P(arr_women < arr_men)
calc_prob <- readRDS("synth_samples/0___synth_res.Rds")
arr_all_men <- ExtractChains(calc_prob, "arr_all_men")
arr_all_women <- ExtractChains(calc_prob,  "arr_all_women")
arr_men_grt_arr_women <- mean(arr_all_men>arr_all_women) # 0.60
arr_women_ineff <- mean(arr_all_women < 0)
rm(calc_prob, arr_all_men, arr_all_women)
saveRDS(arr_men_grt_arr_women, file = "data/prob_arr_men_grt_arr_women.Rds")
saveRDS(arr_women_ineff, file = "data/arr_women_ineff.Rds")
