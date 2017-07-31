# 04_estimate_gender_specific_tx
source ("scripts/00_functions.r")
library(coda)
analysis <- "jags_samples_main"

if(!dir.exists("model_summaries")) dir.create ("model_summaries")
  
  ## 1 Prasugrel/clopidogrel trials in isolation
  load(paste0(analysis, "/prasugrel_only_women_men.Rdata"))
  isln_women <- ExtractChains(LINE.out, "women")
  isln_inter <- ExtractChains(LINE.out, "wd_delta")
  rm(LINE.out)
  
  ## 2 Wider set of trials, 
  # estimate just the effect in men from only 3 trials
  load(paste0(analysis, "/prasugrel_only_men.Rdata"))
  men_only <- ExtractChains(LINE.out)
  rm(LINE.out)
  
  ## 2a add the interaction efficacy from the random effects model
  load(file = paste0(analysis, "/random_effects.Rdata"))
  re_inter <- ExtractChains(LINE.out, "wd_delta")
  rm(LINE.out)
  
  ## 2b add the interaction efficacy from the fixed effects model
  load(file = paste0(analysis, "/fixed_effects.Rdata"))
  fixed_inter <- ExtractChains(LINE.out, "wd")
  rm(LINE.out)
  
  ## 2c add the interaction efficacy from the random effects, treatment stratified model
  load(file = paste0(analysis, "/random_tx_strat.Rdata"))
  # clopidogrel_placebo_ACS clopidogrel_placebo_Stroke  prasugrel_clopidogrel_ACS 
  LINE.out <- LINE.out[ , varnames(LINE.out) %in% c("wd_delta[1]", "wd_delta[2]", "wd_delta[3]")]
  summary(LINE.out)
  re_tx_inter <- ExtractChains(LINE.out, "wd_delta[3]")
  rm(LINE.out)
  
  ## 2d add the interaction efficacy from the fixed effects, treatment stratified model
  load(file = paste0(analysis, "/fixed_tx_strat.Rdata"))
  LINE.out <- LINE.out[ , varnames(LINE.out) %in% c("wd_delta[1]", "wd_delta[2]", "wd_delta[3]")]
  summary(LINE.out)
  fx_tx_inter <- ExtractChains(LINE.out, "wd_delta[3]")
  rm(LINE.out)
  
  ## 2e add the interaction efficacy using data from the Lau paper
  load(file = paste0(analysis, "/lau.Rdata"))
  lau <- ExtractChains(LINE.out, "wd_delta")
  rm(LINE.out)
  
  ## 2f add the trials in both systematic reivews
  load(file = paste0(analysis, "/lau_sa.Rdata"))
  lau_sa <- ExtractChains(LINE.out, "wd_delta")
  rm(LINE.out)
  
  # Arrange interactions as a dataframe
  inters <- cbind(isln_inter, re_inter, fixed_inter, re_tx_inter, fx_tx_inter, lau, lau_sa)
  rm(isln_inter, re_inter, fixed_inter, re_tx_inter, fx_tx_inter, lau, lau_sa)
  
  # Calculate women specific effect by adding male effect to interactions
  women <- apply(inters, 2, function(x) x + men_only)

  # Simplify names of interactions and women specific effects
  colnames(inters) <- gsub("_inter", "", colnames(inters))
  colnames(women) <- gsub("_inter", "", colnames(women))
  # Add in effect in women alone to both interactions and women
  women <- cbind (isln_women, women)
  
  # Summarise results quantiles, and means. Replace median with mean for central estimate
  men_q <- quantile(men_only, probs = c(0.5, 0.025, 0.975))
  men_q[1] <- mean(men_only)
  men_q   <- round(exp(men_q),2)
  
  inter_q <- apply(inters, 2, quantile, probs = c(0.5, 0.025, 0.975))
  inter_q[1,] <- apply(inters, 2, mean)
  inter_q <- round(exp(inter_q),2)
  
  women_q <- apply(women, 2, quantile, probs = c(0.5, 0.025, 0.975))
  women_q[1,] <- apply(women, 2, mean)
  women_q <- round(exp(women_q),2)
  
  # Exponentiate and round both chains (note already taken mean)
  inters <- round(exp(inters),2)
  women <- round(exp(women),2)
  men_only <- round(exp(men_only),2)
  
  # Summarise results probabilities 
  # Probability that treatment less effective in women than men
  inter_p <- apply(inters, 2, function (x) round(mean(x>1),2))
  inter_p
  # probability than treatment harmful in women
  women_p <- apply(women, 2, function(x) round(mean(x>1),2))
  save (men_q, inter_q, women_q, inter_p, women_p, file = paste0("model_summaries/", analysis, ".Rdata"))

## Bleeding ----

## estimate the interaction efficacy for bleeding from the random effects model
load(file = "jags_samples_main/bleed_re.Rdata")
bleed_re_inter <- ExtractChains(LINE.out, "wd_delta")
rm(LINE.out)

## estimate the bleeding effect for prasugrel in ACS only
load(file = "jags_samples_main/prasugrel_only_men_bleed.Rdata")
bleed_pras <- ExtractChains(LINE.out)
rm(LINE.out)

## Estimate the bleeding in women from the RE model
women_bleed <- bleed_re_inter + bleed_pras

## Summarise results for section
inters_bleed <- round(exp(bleed_re_inter),2)
women_bleed <- round(exp(women_bleed),2)
men_bleed <- round(exp(bleed_pras),2)

inter_bleed_q <- quantile(inters_bleed, probs = c(0.5, 0.025, 0.975))
women_bleed_q <- quantile(women_bleed, probs = c(0.5, 0.025, 0.975))
men_bleed_q <-   quantile(men_bleed,   probs = c(0.5, 0.025, 0.975))

save (men_bleed_q, inter_bleed_q, women_bleed_q,
      file = "model_summaries/bleeding.Rdata")
