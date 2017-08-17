## Run baseline risk model and evidence synthesis.
library(rjags)
library(stringr)
library(tidyverse)
load.module("glm")

### Switches
run_test <- FALSE # Set to TRUE if want to test results are as expected

### Create directory
if(!dir.exists("synth_smry")) dir.create ("synth_smry")
if(!dir.exists("synth_samples")) dir.create ("synth_samples")

## Test datasets ----
# 1 using dataset with identical baseline risks in men and women, increases with age
br_test <- expand.grid(age = seq(30, 90,10), gender = c("Male", "Female"), stringsAsFactors = FALSE)
br_test <- br_test %>%
  mutate(n = 1000,
         death = n * round(plogis(-3 + (age/10) *0.2 ),2),
         bleed = round(0.05*death),
         cv    = round(0.8* death),
         other = death - bleed - cv,
         surv = n - death,
         age = paste0(age, "-", age +9)) %>%
    select(-death)

# Using dataset with much higher baseline risk in women, increases with age
br_test <- expand.grid(age = seq(30, 90,10), gender = c("Male", "Female"), stringsAsFactors = FALSE)
br_test <- br_test %>%
  mutate(n = 1000,
         death = n * round(plogis(-3 + (age/10) *0.2 + 0.69*(gender == "Female")),2),
         bleed = round(0.05*death),
         cv    = round(0.8* death),
         other = death - bleed - cv,
         surv = n - death,
         age = paste0(age, "-", age +9)) %>%
  dplyr::select(-death)


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
for (aggrgt in c("", "avd_non_cllps")){
  count <- 1 # Set to one to capture items only needed once for each dataset
  # Read in baseline dataset or read in fake TEST VERSION of data
  if(run_test == TRUE)  br <- br_test else{
    br <- read.csv("data/baseline_risk_actual_data.csv", as.is = TRUE)
  }
  br$age <- as.numeric(substr(br$age,1,1))

  br$gender <- ifelse(br$gender == "Male", 1, 0)
  ## Aggregates data as per looping
  if (aggrgt == "avd_non_cllps"){
    br <- br %>%
      group_by(gender) %>%
      summarise (n = sum(n), age = mean(age), bleed = sum(bleed), cv = sum(cv), 
                 other = sum(other), surv = sum(surv))
    br <- br[2:1,] # need to re-order so that male is first, otherwise it will get results wrong
  }
  # create matrix of outcomes to loop through
  y <- as.matrix(br [ , c("bleed", "cv", "other", "surv")])
  # Retain original rr for effect of treatment on mace in women
  f_mace_original <- f_mace
  ## Loop through bleeding sensitivity analyses
  for (bleed_effect in row.names(f_bleed)){
    for(gender_sa in c("null_mace", "")){
      if(gender_sa == "null_mace") f_mace[] <- m_mace else f_mace[] <- f_mace_original
    ## Run jags ----
      jags <- jags.model('jags_code/baseline.txt',
                         data = c(list (numrow = length(br$age),
                                        end_male = length(br$age)/2,
                                        st_female = (length(br$age)/2)+1,
                                      y = y,
                                      n =br$n,
                                      age = br$age,
                                      gender = br$gender),
                                  as.list(m_mace),
                                  as.list(f_mace),
                                  as.list(m_bleed),
                                  as.list(f_bleed[bleed_effect,])
                         ),
                         n.chains = 2,
                         n.adapt = 1000)
      
      update(jags, burnin) # Burn in
      
      data(LINE)
      LINE$recompile()
      LINE.out <- coda.samples(jags, c("arr_bleed", "arr_cv", "arr_other","arr_all",
                                       "arr_bleed_men", "arr_cv_men", "arr_other_men","arr_all_men",
                                       "arr_bleed_women", "arr_cv_women", "arr_other_women","arr_all_women",
                                       "arr_sex_age", "arr_sex_age_all",
                                       "y_chk",
                                       "rr_cv_men_ss", "rr_cv_women_ss",
                                       "rr_cv_men_mrg", "rr_cv_women_mrg"),
                               n.iter = iterations, n.thin = 20)
      ## Effect of non-collapsability appears to be small in that applying RR to each stratum 
      ## leads to an RR which is similar to the overall RR 
      ## Summarise arr
      res <- summary(LINE.out)
      
      a <- res$quantile
      # Replace median with mean for all subsequent calculations. 
      # Was previously using median to summarise, now use mean
      a[, "50%"] <- res$statistics[, "Mean"]
      
      # rm(LINE.out)
      gc()
      ## Summarise results of evidence synthesis
      QuickSmrsFx<- function (srch_string = "^rr_cv"){
        my_matrix <- a[ grepl(srch_string, rownames(a)) ,]
        my_matrix <- as.data.frame(my_matrix)
        names(my_matrix) <- gsub(".", "_", names(my_matrix), fixed = TRUE)
        names(my_matrix) <- gsub("%", "", names(my_matrix), fixed = TRUE)
        names(my_matrix) <- paste0("q", names(my_matrix))
        my_matrix[] <- lapply(my_matrix, function(x) format(round(x,2), nsmall = 2))
        my_matrix$res <- paste0(my_matrix$q50, " (", my_matrix$q2_5, "-", my_matrix$q97_5, ")")
        my_matrix
      }
      rrs <- QuickSmrsFx()
      arr_sex_age_all <-  QuickSmrsFx("^arr_sex_age_all") 
      a <- a[!grepl("^arr_sex_age_all|^rr_cv", rownames(a)) , ]
      smry <- as.data.frame(a)
      # Use indexing from summary to create i and j columns
      dtls <- str_split_fixed(rownames(a), "\\[", n = 2)
      dtls2 <- str_split_fixed(dtls[,2], ",", n = 2)
      dtls2[,2] <- gsub ("]", "", dtls2[,2])
      dtls <- cbind (dtls[,1], dtls2)
      rm(dtls2)
      dtls <- as.data.frame(dtls, stringsAsFactors = FALSE)
      names(dtls) <- c("param", "i", "j")
      dtls$i <- as.numeric(dtls$i)
      dtls$j <- as.numeric(dtls$j)
      
      smry <- cbind (smry, dtls)
      
      y_chk <-  smry [ smry$param == "y_chk",] 
      gender_age_spec <- smry [ smry$param == "arr_sex_age", ]
      
      arr_smry <- smry [ !smry$param %in% c("y_chk", "arr_sex_age"),]
      arr_smry$i <- NULL
      arr_smry$j <- NULL
      dtls <- str_split_fixed(arr_smry$param, "_", n = 3)
      dtls[dtls[,3] == "" ,3] <- "both"
      colnames(dtls)[2:3] <- c("outcome", "gender")
      arr_smry <- cbind(arr_smry, dtls[,2:3])
      arr_smry$param <- NULL
      arr_smry <- arr_smry [ , c("gender", "outcome","2.5%", "50%", "97.5%")]
      names(arr_smry) <- c("gender", "outcome", "lci", "est", "uci")
      
      nnt_smry <- arr_smry
      nnt_smry$lci <- 100/arr_smry$uci # note deliberate switch
      nnt_smry$uci <- 100/arr_smry$lci
      nnt_smry$est <- 100/arr_smry$est
      nnt_smry [ , c("est", "lci", "uci")] <- lapply(nnt_smry [ , c("est", "lci", "uci")],
                                                     function (x) format(round(x), trim = TRUE))
      
      arr_smry [ , c("est", "lci", "uci")] <- lapply(arr_smry [ , c("est", "lci", "uci")],
                                                     function (x) format(round(x,2), trim = TRUE))
      
      # Good fit to the data. Same each time so check only once per dataset
      while (count ==1) {
      names(y_chk) <- c("lci", "q25", "est", "q75", "uci", "param", "i", "j")
      otype_lookup <- colnames(y)
      y_chk$otype <- otype_lookup[y_chk$j]
      y_chk <- cbind(br[ , c("age", "gender", "n")], y_chk)
      save(y_chk, file = paste0("synth_smry/","baseline risks with ",  aggrgt, "data.Rdata"))
      count <- count + 1
      }
      # Estimate age and gender-specific, will also work where aggregated by age
      gender_age_spec[ , c("2.5%", "50%", "97.5%")] <- lapply (gender_age_spec[ , c("2.5%", "50%", "97.5%")],
                                                               function(x) format(round((x),2),nsmall = 2, trim = TRUE))
      
      gender_age_spec$res <- paste0(gender_age_spec$`50%`, " (", gender_age_spec$`2.5%`, "-", gender_age_spec$`97.5%`, ")")
      gender_age_spec_est_lci_uci <- gender_age_spec[ , c("res", "2.5%", "50%", "97.5%")]
      gender_age_spec <-  reshape2::dcast (gender_age_spec, i ~ j, value.var = "res")
      gender_age_spec$i <- NULL
      names(gender_age_spec) <- c("Bleeding", "MACE", "Competing")
      gender_age_spec <- cbind (br [ , c("age", "gender")], gender_age_spec )
      gender_age_spec$age <- paste0(10* gender_age_spec$age, "-", 10* gender_age_spec$age+9)
      gender_age_spec$gender <- factor(gender_age_spec$gender, levels = c(0,1), labels = c("Women", "Men"))
    
      # Add in totals column for estimates with 95%CI as string and as separate elements
      gender_age_spec$Total <- arr_sex_age_all$res
      arr_sex_age_all <- arr_sex_age_all[ , c("res", "q2_5", "q50", "q97_5")]
      names(arr_sex_age_all) <- c("res", "2.5%", "50%", "97.5%")
      rownames(arr_sex_age_all) <- gsub(pattern = "]", ",4]", rownames(arr_sex_age_all))
      gender_age_spec_est_lci_uci <- rbind(gender_age_spec_est_lci_uci, arr_sex_age_all)
      # Save age-specific and gender and age specific, save versions for aggregated and non_aggregated data
      label_file <- paste(bleed_effect, aggrgt, gender_sa, sep = "_")
      save(rrs, file = paste0("synth_smry/summary_of_rate_ratio", label_file,".rdata"))
      save(arr_smry, nnt_smry, file = paste0("synth_smry/summary_of_arr", label_file, ".rdata"))
      save(gender_age_spec, gender_age_spec_est_lci_uci,
           file = paste0("synth_smry/summary_of_gender_age_spec_arr_", label_file, ".rdata"))
      saveRDS(LINE.out, file = paste0("synth_samples/", label_file, "_synth_res.Rds"))
      print(paste0("Finished ", bleed_effect, "% worse bleeding in women ", aggrgt, gender_sa))
} # end of loop through MACE sensitivity analysis setting same RR for men and women
} # end of loop through bleeding RR in women sensitivity analyses
} # end of loop through whether aggregated by age or not prior to performing analysis
