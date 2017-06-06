# 05_extract_stratified
# Extract stratifed estimates 
source("scripts/00_functions.r")
load(file = "jags_samples_main/bleed_strat.Rdata")
bleed_trials <- ExtractStratified(LINE.out = LINE.out, mydf = "bleed")
save(bleed_trials, file = "model_summaries/bleeding_stratified.Rdata")

load("jags_samples_main/stratified.Rdata")
trials <- ExtractStratified(LINE.out = LINE.out, mydf = "main_wide")
save (trials, file = paste0("model_summaries/", "stratified_each_trial", ".Rdata"))
