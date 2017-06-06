### Run all scripts

## Run sex-treatment interaction models, note also runs reading in in data and creating JAGS models
diagnostic <- FALSE
source("scripts/03_Run_JAGS_models.R")

diagnostic <- TRUE
source("scripts/03_Run_JAGS_models.R")

## Extract results from sex-treatment interaction models and summarise
source("scripts/04_estimate_gender_specific_tx.R")
source("scripts/05_extract_stratified.R")
source("scripts/06_Summarise_results_using_non_linear_models.R")

## Run baseline risk model and model treatment effects in Scottish population
source("scripts/07a_baseline_risks_and_ES.R")
source("scripts/07b_baseline_risks_and_ES_sens_analysis.R")

## Create summary tables and plots
source("scripts/08a_process_plots_tables.R")
source("scripts/08b_process_complicated_plot.R")

## Create diagnostics plots
source("scripts/09_Run_model_diagnostics.R")

## Create results for paper
## Easiest to knit within Rstudio
"scripts/10_paper.Rmd"
"scripts/11_supplementary_appendix.Rmd"



