## 06 Summarise results using non-linear models
source ("scripts/00_functions.r")
library(MASS)
library(rjags)
library(tidyverse)

## Quick summarising function
quickfx <- function(x) exp(quantile(x, probs = c(0.025, 0.5, 0.975))) %>% round(3)

## Effect of tx on MACE ----
# Note already performed this calculation in "04_estimate_gender_specific_tx.r"
# But there I stored it as quantiles, not the full chains

## Sum pooled (RE model) interaction and effect in men to esitmate effect in women
load(file = paste0("jags_samples_main/", "random_effects.Rdata"))
re_inter <- ExtractChains(LINE.out, "wd_delta")
rm(LINE.out)

## Extract treatment effect estimate for men
load(paste0("jags_samples_main/", "prasugrel_only_men.Rdata"))
men_only <- ExtractChains(LINE.out, "men")
rm(LINE.out)

# Calculate effect in women from interaction in RE model
tx_f_mace <- re_inter + men_only
tx_m_mace <- men_only
MakeDens(list(women = tx_f_mace, men = tx_m_mace)) %>% 
  MakeDensPlot()

## Calculate probability that treatment in women (in relative terms) 
## is less effective than the treatment in men
rr_prob_re <- mean(re_inter>0) # higher worse
saveRDS(rr_prob_re, file = "data/rr_prob_re.RDS")

## Effect of tx on Bleeding ----
## Note bleeding interaction uses a strong prior
load(file = "jags_samples_main/bleed_re.Rdata")
bleed_re_inter <- ExtractChains(LINE.out, "wd_delta")
rm(LINE.out)

load(file = "jags_samples_main/prasugrel_only_men_bleed.Rdata")
bleed_men <- ExtractChains(LINE.out, "men")

tx_m_bleed <- bleed_men
tx_f_bleed <- bleed_re_inter + bleed_men

## Generate bleeding sensitivity analyses
bleed_sa_20 <- rnorm(length(bleed_re_inter), log(1.2), log(1.25^(1.96^-1)))
bleed_sa_30 <- rnorm(length(bleed_re_inter), log(1.3), log(1.25^(1.96^-1)))
bleed_sa_40 <- rnorm(length(bleed_re_inter), log(1.4), log(1.25^(1.96^-1)))
bleed_sa_50 <- rnorm(length(bleed_re_inter), log(1.5), log(1.25^(1.96^-1)))

tx_f_bleed_20 <- bleed_sa_20 + bleed_men
tx_f_bleed_30 <- bleed_sa_30 + bleed_men
tx_f_bleed_40 <- bleed_sa_40 + bleed_men
tx_f_bleed_50 <- bleed_sa_50 + bleed_men

rm(bleed_re_inter, bleed_men, bleed_sa_20, bleed_sa_30,
   bleed_sa_40, bleed_sa_50)

## Approximate with nls ----
f_list <- as.list(paste0("tx_f_bleed_", seq(20,50,10)))
mylist <- c(list("tx_f_mace", "tx_m_mace", "tx_f_bleed", "tx_m_bleed"), f_list)
rm(f_list)
names(mylist) <- as.character(mylist)
mylist <- lapply(mylist, Tmod)
mylist <- lapply( mylist, function (x) {
  names(x) <- stringr::str_sub(names(x), 4)
  x}
  )
names(mylist) <- stringr::str_sub(names(mylist), 4)
list2env(mylist, envir = .GlobalEnv)


## Save results
save(f_mace, m_mace, f_bleed, m_bleed, f_bleed_20, f_bleed_30, f_bleed_40, f_bleed_50, 
     file = "data/tx_effects_distr.Rdata")

