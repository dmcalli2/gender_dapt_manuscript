#07b_sensitivity analyses

## Run baseline risk model and evidence synthesis.
library(rjags)
library(stringr)
library(tidyverse)
load.module("glm")

### Create directory
if(!dir.exists("synth_sa_smry")) dir.create ("synth_sa_smry")

## Choose iterations and burnin ----
iterations <- 20000
burnin <- 10000

# Write model ----
modelstring <- "
## Baseline risk 
model{
  for (i in 1:numrow){
    y[i,1:4] ~ dmulti (p[i,1:4], n[i])
    p[i,4] <- 1 - sum(p[i,1:3]) 
    slam[i] <- sum(lamda[i,1:3]) 
    cumslam[i] <- 1 - exp(-slam[i]) 
    for (j in 1:3){
      p[i,j] <- lamda[i,j] * cumslam[i] /slam[i] 
      log(lamda[i,j]) <- alpha0[j] + alpha1[j]*age[i] + alpha2[j]*gender[i] + alpha3[j]*gender[i]*age[i]
    } # outcomes
  } # rows
  # Priors for each outcome
  for (j in 1:3){
    alpha0[j] ~ dnorm(0, 0.0001)
    alpha1[j] ~ dnorm(0, 0.0001)
    alpha2[j] ~ dnorm(0, 0.0001)
    alpha3[j] ~ dnorm(0, 0.0001)
  }
  
  ## Represent treatment effects - parameters for t-distributions supplied as data to model
  # Men
  m_mace ~ dt(m_mace_m, m_mace_prec, m_mace_df)
  m_mace_rr <- exp(m_mace)
  m_bleed ~ dt(m_bleed_m, m_bleed_prec, m_bleed_df)
  m_bleed_rr <- exp(m_bleed)
  # Women
  f_mace ~ dt(f_mace_m, f_mace_prec, f_mace_df)
  f_mace_rr <- exp(f_mace)
  f_bleed ~ dt(f_bleed_m, f_bleed_prec, f_bleed_df)
  f_bleed_rr <- exp(f_bleed)
  
  ## Apply treatment effects
  for (i in 1:numrow){
    # Make new rates for each outcome; gender = 1 for male
    lamda_new [i,1] <- lamda [i, 1] * m_bleed_rr^gender[i] * f_bleed_rr^(1-gender[i]) # Bleeding, 
    lamda_new [i,2] <- lamda [i, 2] * m_mace_rr^gender[i] *  f_mace_rr^(1-gender[i]) # CV risk
    lamda_new [i,3] <- lamda [i, 3] # ie x 1 for non-CV, non bleeding
    # Calculate cumulative incidence and new event numbers based on the intervention for each outcome
    for (j in 1:3){
      p_new[i,j] <- lamda_new[i,j] * cumslam_new[i] /slam_new[i] 
      #arr_scotland[i,j] <- y[i,j] - y_new[i,j] # Use for estimating events in a Scotland-specific prediction
      arr[i,j] <- (p[i,j] - p_new[i,j]) *n[i] #  Difference in proportions, multiplied for weighting
      arr_sex_age[i,j] <- 100 * arr[i,j] / n[i] # age and sex specific ARR, each outcome type
    }# end of loop through outcomes
    # Calculate total proportion of deaths for all outcome types
      arr_sex_age_all[i] <- 100 * sum(arr[i,]) / n[i] # age and sex specific ARR, total
  ## Predict proportions and outcomes based on new model
    y_chk[i,1:4] ~ dmulti (p[i,1:4], n[i]) # Posterior predictive check
    #y_new[i,1:4] ~ dmulti (p_new[i,1:4], n[i]) # Predicted events under treatment
    p_new[i,4] <- 1 - sum(p_new[i,1:3]) # Predicted risk under treatment
    slam_new[i] <- sum(lamda_new[i,1:3]) 
    cumslam_new[i] <- 1 - exp(-slam_new[i])
    
    ## Calculate person time for each row, by subtracting 6 months for any event
    #pt_old[i] <- 1 * n[i] - 0.5 * sum(y_chk[i, 1:3]) # person time is at one-year 
    #pt_new[i] <- 1 * n[i] - 0.5 * sum(y_new[i, 1:3]) # person time is at one-year 
  }# end of loop through rows
  
  ## Summary statistics across rows
  # Sum of age-sex specific absolute risk reductions weighted by the number in each stratum
  # for whole dataset
  arr_bleed <- 100*sum (arr[,1]) / sum(n[]) # note arr is actually a count not a proportion
  arr_cv <-    100*sum (arr[,2]) / sum(n[])
  arr_other <- 100*sum (arr[,3]) / sum(n[])
  arr_all <- (arr_bleed + arr_cv + arr_other) 
  # for men (first half of dataset)
  arr_bleed_men <- 100*sum (arr[1:end_male,1]) / sum(n[1:end_male]) # note arr is actually a count
  arr_cv_men <-    100*sum (arr[1:end_male,2]) / sum(n[1:end_male])
  arr_other_men <- 100*sum (arr[1:end_male,3]) / sum(n[1:end_male])
  arr_all_men <- (arr_bleed_men + arr_cv_men + arr_other_men) 
  # for women (second half of dataset)
  arr_bleed_women <- 100*sum (arr[st_female:numrow,1]) / sum(n[st_female:numrow])
  arr_cv_women <-    100*sum (arr[st_female:numrow,2]) / sum(n[st_female:numrow])
  arr_other_women <- 100*sum (arr[st_female:numrow,3]) / sum(n[st_female:numrow])
  arr_all_women <- (arr_bleed_women + arr_cv_women + arr_other_women) 
  # Summarise rate ratios for men and women
  # (note this is determined by the model inputs, only calculated in JAGS for convenience),
  rr_bleed_men_ss <- lamda_new[1,1] / lamda [1, 1]
  rr_cv_men_ss    <- lamda_new[1,2] / lamda [1, 2]
  rr_other_men_ss <- lamda_new[1,3] / lamda [1, 3]
  rr_bleed_women_ss <- lamda_new[st_female,1] / lamda [st_female, 1]
  rr_cv_women_ss    <- lamda_new[st_female,2] / lamda [st_female, 2]
  rr_other_women_ss <- lamda_new[st_female,3] / lamda [st_female, 3]
}"
writeLines(modelstring,con="jags_code/baseline.txt")

# Read in treatment effects ----
load(file = "data/tx_effects_distr.Rdata")

# Put treatment effects for women into dataframe so can loop through
f_bleed <- rbind(f_bleed, f_bleed_20, f_bleed_30, f_bleed_40, f_bleed_50)
row.names(f_bleed) <- c(0, seq(20, 50, 10))
rm(f_bleed_20, f_bleed_30, f_bleed_40, f_bleed_50)

# Loop through dataset aggregated by age, and left as original strata to examine potential
# effect of non-collapsability of rate ratios on ARR
br <- read.csv("data/baseline_risk_actual_data.csv", as.is = TRUE)
br$age <- as.numeric(substr(br$age,1,1))
br$gender <- ifelse(br$gender == "Male", 1, 0)

br <- br %>% 
  group_by(gender) %>% 
  summarise_all(sum) %>% 
  arrange(desc(gender)) %>% 
  mutate(age = 0)

## Sensitivity analyses
sa <- expand.grid(men = c(1, 5, 10), women = c(1,1.5,2))

br_sa <- lapply(seq_along(sa$men), function (i) {
  sa_select <- sa[i, , drop = FALSE]
  br$bleed_multiply <- ifelse(br$gender ==1, sa_select$men * sa_select$women, sa_select$women)
  br$bleed_new <- round(br$bleed * br$bleed_multiply)
  br$surv <- br$surv - br$bleed_new + br$bleed
  br$bleed <- br$bleed_new
  br$bleed_new <- NULL
  br$bleed_multiply <- NULL
  br$surv <- pmax(0, br$surv)
  br
})

br_sa_res <- lapply(br_sa, function (br_sa_select) {
  y <- as.matrix(br_sa_select[ , c("bleed", "cv", "other", "surv")])
  jags <- jags.model('jags_code/baseline.txt',
                         data = c(list (numrow = length(br_sa_select$age),
                                        end_male = length(br_sa_select$age)/2,
                                        st_female = (length(br_sa_select$age)/2)+1,
                                      y = y,
                                      n =br_sa_select$n,
                                      age = br_sa_select$age,
                                      gender = br_sa_select$gender),
                                  as.list(m_mace),
                                  as.list(f_mace),
                                  as.list(m_bleed),
                                  as.list(f_bleed[1,])
                         ),
                         n.chains = 2,
                         n.adapt = 1000)
      
      update(jags, burnin) # Burn in
      
      data(LINE)
      LINE$recompile()
      LINE.out <- coda.samples(jags, "arr_all_women",
                               n.iter = iterations, n.thin = 20)
      summary(LINE.out)
})
saveRDS(br_sa_res, file = "synth_sa_smry/sensitivity_analyses.Rds")

pt_est <- map_dbl(br_sa_res, function(x) x$statistics["Mean"])
pt_est_rng <- range(pt_est) %>%  round(2) %>%  paste0("%") %>%  paste(collapse = " to ")
saveRDS(pt_est_rng, file = "synth_sa_smry/pt_est_rng.Rds")

      


