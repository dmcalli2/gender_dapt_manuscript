# 03_Run_JAGS_models
# Runs jags models and save samples

source("scripts/02_model_strings.R")

# Switches for diagnostic plots
# diagnostic <- FALSE

## Load packages ----
library(rjags)
library(dplyr)
library(glmmBUGS)

## Load JAGS GLM module to allow block updating ----
load.module("glm")
load.module("dic")

## Create folders for jags models to be saved ---
if(!dir.exists("jags_samples_main")) dir.create ("jags_samples_main")

# List files and make vector
filenames <- c("bleed_re.txt", "bleed_strat.txt", "fixed_effects.txt", "fixed_tx_strat.txt", 
               "lau.txt", "lau_sa.txt", "prasugrel_only_men.txt", "prasugrel_only_men_bleed.txt", 
               "prasugrel_only_women_men.txt", "random_effects.txt", "random_tx_strat.txt", 
               "stratified.txt")

names(filenames) <- filenames

## Sampling, burnin etc for models----
adapt <- 2000
burnin <- 20000
iterations <- 50000
samples_for_r <- 10000

if (diagnostic == TRUE) {
  burnin <- 1
  iterations <- 20000
  samples_for_r <- iterations # ie no thinning
}

################################################### Analyses ----
foldername <- "main"
  
# Read and run data processing and model text based on choice of analysis---
source ("scripts/01_read_data.R")

# Choose which models to run ----
# Selects desired model assigns filename; can loop if want
for (filename in filenames){
  print(filename)
  # for (filename in "random_effects.txt"){
  outname <- gsub(".txt", ".Rdata", filename)
  dicname <- gsub(".txt", "_dic.Rdata", filename)

  if(diagnostic == FALSE) outname <- paste0("jags_samples_", foldername , "/", outname)
  if(diagnostic == TRUE)  outname <- paste0("jags_samples_", foldername , "/", "dgnstc_", outname)
  dicname <- paste0("jags_samples_", foldername , "/", dicname)
  print (paste0("JAGS model running for ", filename, " model will save to ", outname))
  
  ## Models with all trials ----
  if (filename %in% trials_only){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(main_wide[ , c('time', 'nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')],
                                list(n_trials = nrow(main_wide))),
                       n.chains = 2,
                       n.adapt = adapt)
    
    update(jags, burnin)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('d', 'd_delta', 'wd', 'wd_delta', 'wd_sd'),
                             iterations, thin = iterations/samples_for_r)
    # comparison of interest "wd_delta"
  }
  
  ## Models with data from Lau et al in JACC ----
  if (filename %in% lau_choose){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(lau[ , c('time', 'nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')],
                                list(n_trials = nrow(lau))),
                       n.chains = 2,
                       n.adapt = adapt)
    
    update(jags, burnin)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('d', 'd_delta', 'wd', 'wd_delta'),
                             iterations, thin = iterations/samples_for_r)
    # comparison of interest "wd_delta"
  }
  
  ## Models with data from both systematic reviews ----
  if (filename %in% lau_sa_choose){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(lau_sa[ , c('time', 'nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')],
                                list(n_trials = nrow(lau_sa))),
                       n.chains = 2,
                       n.adapt = adapt)
    
    update(jags, burnin)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('d', 'd_delta', 'wd', 'wd_delta'),
                             iterations, thin = iterations/samples_for_r)
    # comparison of interest "wd_delta"
  }
  
  
  ## Models for prasugrel versus clopi in ACS men and women  ----
  if (filename %in% pras_women_men){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(n_trial = nrow(pras), 
                                pras [ , c('time', 'nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')]),
                       n.chains = 2,
                       n.adapt = adapt)
    
    update(jags, burnin)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('mu','women', 'wd_delta', 'd_delta', 'men'),
                             iterations, thin = iterations/samples_for_r)
  }
  
  ## Models for prasugrel versus clopi in ACS men  ----
  if (filename %in% pras_men){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(n_trial = nrow(pras), 
                                pras [ , c('time', 'nmt', 'nmc', 'rmt', 'rmc')]),
                       n.chains = 2,
                       n.adapt = adapt)
    update(jags, burnin)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('men'),
                             iterations, thin = iterations/samples_for_r)
  }
  
  ## BLEEDING Models for prasugrel versus clopi in ACS men  ----
  if (filename %in% pras_men_bleed){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(n_trial = nrow(bleed_pras_men), 
                                bleed_pras_men [ , c('time', 'nmt', 'nmc', 'rmt', 'rmc')]),
                       n.chains = 2,
                       n.adapt = 1000)
    update(jags, 1000)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('men'),
                             iterations, thin = iterations/samples_for_r)
  }
  
  ## BLEEDING Models ----
  if (filename %in% bleed_re_strat){
    jags <- jags.model(paste0("jags_code/", filename),
                       data = c(bleed[ , c('time', 'nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')],
                                list(n_trials = nrow(bleed))),
                       n.chains = 2,
                       n.adapt = adapt)
    
    update(jags, burnin)
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('d', 'd_delta', 'wd', 'wd_delta'),
                             iterations, thin = iterations/samples_for_r)
    # comparison of interest "wd_delta"
  }
  

  ## Nested model treatment and comparisons - there are three of these
  if(filename %in% tx_strat){
    jags <- jags.model (paste0("jags_code/", filename),
                        data = c(list(Ncompar_indic = wide_rag_flat$Ncompar_indic,
                                      Scompar_indic = wide_rag_flat$Scompar_indic),
                                 main_wide_flat [ , c('time', 'nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')]),
                        n.adapt = adapt,
                        n.chains = 2)
    
    update(jags, burnin) # burn-in
    data(LINE)
    LINE$recompile()
    LINE.out <- coda.samples(jags,
                             c('md', 'wd', 'wd_delta', 'wd_delta_mu', 'sd1'),
                             iterations, thin = iterations/samples_for_r)
    # Comparison of interest "wd_delta[4]"
  }
  
  ## Save file for each model in relevant folder----
  # Get DIC for each model, only for models which have converged
  if(diagnostic == FALSE) {
    mydic <- dic.samples (jags, n.iter = 10000, type = "pD")
    save(mydic, file = dicname)
  }
  save (LINE.out, file = outname)
  print ("Expecting warning - Failed to set trace monitor as using same code for multiple models")
  rm(LINE.out)
  
}# End of looping through models
