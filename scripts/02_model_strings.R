# 02 Write models as strings

## Folder for jags model descriptions to be saved
if(!dir.exists("jags_code")) dir.create ("jags_code")

## Random effects model, ignoring treatment and condition "random_effects.txt" ----
modelstring <- "
model{
for (z in 1:n_trials){
   # likelihood
   rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
   rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
   rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
   rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated 

   # link and linear predictor
   cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
   cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]                # baseline in women
   cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         d[z]         # treatment in men
   cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + d[z] + wd[z] # treatment in women

   # Trial specific priors
   mu[z] ~ dnorm (0, 0.0001)
   wu[z] ~ dnorm (0, 0.0001)
   d[z]  ~ dnorm (0, 0.0001)
   # Shared priors (random effects)
   wd[z] ~ dnorm(wd_delta, wd_tau)
} #end of trials
   # Meta-analysis level priors
 
   wd_delta ~ dnorm(0, 0.0001)
   wd_tau <-  1/wd_var
   wd_var ~ dlnorm(-3.23, 1/1.88^2)
   wd_sd <- wd_var^0.5

   #Residual deviance for all trials
  # totresdev <- sum(resdev[]) #Total Residual Deviance
}# end of model
"
writeLines(modelstring, con = "jags_code/random_effects.txt")

## Fixed effects model, ignoring treatment and condition "fixed_effects.txt" ----
modelstring <- "
model{
  for (z in 1:n_trials){
    # likelihood
    rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
    rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
    rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
    rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated 
    
    # link and linear predictor
    cloglog(pmc[z]) <- log(time[z]) +  mu[z]                   # baseline in men
    cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]           # baseline in women
    cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         d       # treatment in men
    cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + d + wd  # treatment in women
    
    # Trial specific priors
    mu[z] ~ dnorm (0, 0.0001)
    wu[z] ~ dnorm (0, 0.0001)
    
    } #end of trials loop
    # Shared priors (fixed effects)
    d  ~ dnorm (0, 0.0001)
    wd ~ dnorm (0, 0.0001)
}# end of model
"
writeLines(modelstring, con = "jags_code/fixed_effects.txt")

## Stratified model, analysing each trial separately "stratified.txt" ----
  modelstring <- "
  model{
  for (z in 1:n_trials){
  # likelihood
  rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
  rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
  rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
  rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated
  
  # link and linear predictor
  cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
  cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]                # baseline in women
  cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         d[z]         # treatment in men
  cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + d[z] + wd[z] # treatment in women
  
  # Trial specific priors
  mu[z] ~ dnorm (0, 0.0001)
  wu[z] ~ dnorm (0, 0.0001)
  d[z]  ~ dnorm (0, 0.0001)
  wd[z] ~ dnorm(0, 0.0001)
  
  } #end of trials
  }# end of model
  "
writeLines(modelstring, con = "jags_code/stratified.txt")


## Restricted to effect in WOMEN and MEN in prasugrel versus clopidogrel trials "prasugrel_only_women_men.txt" ----
modelstring <- "
model{
for (z in 1:n_trial){
# likelihood
rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated

# link and linear predictor
cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]                # baseline in women
cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         d[z]         # treatment in men
cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + d[z] + wd[z] # treatment in women
# Trial specific priors
mu[z] ~ dnorm (0, 0.0001)
wu[z] ~ dnorm (0, 0.0001)
# d[z]  ~ dnorm (0, 0.0001)
# Shared priors (random effects)
d[z] ~ dnorm (d_delta, d_tau)
wd[z] ~ dnorm(wd_delta, wd_tau)
} #end of trials
# Meta-analysis level priors
d_delta ~ dnorm(0, 0.0001)
d_tau <-  1/d_var
d_var ~ dlnorm(-3.23, 1/1.88^2)

wd_delta ~ dnorm(0, 0.0001)
wd_tau <-  1/wd_var
wd_var ~ dlnorm(-3.23, 1/1.88^2)

# effect in women
women <- d_delta + wd_delta
men <- d_delta
}# end of model
"
writeLines(modelstring, con = "jags_code/prasugrel_only_women_men.txt")

## Restricted to effect in MEN ONLY for prasugrel versus clopidogrel trials  "prasugrel_only_men.txt"  ----
modelstring <- " model{
for (z in 1:n_trial){
  # likelihood
  rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
  rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
  # link and linear predictor
  cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
  cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         d[z]         # treatment in men
  # Trial specific priors
  mu[z] ~ dnorm (0, 0.0001)
  wu[z] ~ dnorm (0, 0.0001)
  # Shared priors (random effects)
  d[z] ~ dnorm (d_delta, d_tau)
  } #end of trials
# Top level priors
d_delta ~ dnorm(0, 0.0001)
d_tau <-  1/d_var
d_var ~ dlnorm(-3.23, 1/1.88^2)
# effect in men
men <- d_delta 
}# end of model
"
writeLines(modelstring, con = "jags_code/prasugrel_only_men.txt")


## AS ABOVE FOR BLEEDING identical to  "prasugrel_only_men.txt",  "prasugrel_only_men_bleed.txt" ----

modelstring <- readLines(con = "jags_code/prasugrel_only_men.txt")
writeLines(modelstring, con = "jags_code/prasugrel_only_men_bleed.txt" )


## Stratified bleeding model - identical to "stratified.txt", called "bleed_strat.txt" ----
modelstring <- readLines(con = "jags_code/stratified.txt")
writeLines(modelstring, con = "jags_code/bleed_strat.txt" )

## Random effects bleeding model - similar to "random_effects.txt", called "bleed_re.txt" ----
# Has informative priors for interaction
modelstring <- "
model{
for (z in 1:n_trials){
   # likelihood
   rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
   rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
   rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
   rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated 

   # link and linear predictor
   cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
   cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]                # baseline in women
   cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         d[z]         # treatment in men
   cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + d[z] + wd[z] # treatment in women

   # Trial specific priors
   mu[z] ~ dnorm (0, 0.0001)
   wu[z] ~ dnorm (0, 0.0001)
   d[z]  ~ dnorm (0, 0.0001)
   # Shared priors (random effects)
   # d[z] ~ dnorm (d_delta, d_tau)
   wd[z] ~ dnorm(wd_delta, wd_tau)
} #end of trials
   # Meta-analysis level priors
   # Meta-analysis level priors
   # d_delta ~ dnorm(0, 0.0001)
   # d_tau <-  1/d_var
   # d_var ~ dlnorm(-3.23, 1/1.88^2)
    
   wd_delta ~ dnorm(0, 1/0.1137755^2) # SD confines effect in women to 1.25 more or less than in men
   wd_tau <-  1/wd_var
   wd_var ~ dlnorm(-3.23, 1/1.88^2)

   #Residual deviance for all trials
  # totresdev <- sum(resdev[]) #Total Residual Deviance
}# end of model
"
writeLines(modelstring, con = "jags_code/bleed_re.txt")

## Random effects model in lau data - identical to "random_effects.txt", called "lau.txt" ----
modelstring <- readLines(con = "jags_code/random_effects.txt")
writeLines(modelstring, con = "jags_code/lau.txt" )

## Random effects model in lau without stable - identical to "random_effects.txt", called "lau_sa.txt" ----
modelstring <- readLines(con = "jags_code/random_effects.txt")
writeLines(modelstring, con = "jags_code/lau_sa.txt" )



## Random effects model stratified by treatment "random_tx_strat.txt" ----
# comparison and indication flattened  - ci (comparison and indication in wide_rag_flat)
# study/trial - z (observations in wide_rag_flat)
# Stratified by treatemnt/comparison combination, NOT pooled across these
modelstring <- "
model{
  for(ci in 1:Ncompar_indic) {
    for(z in Scompar_indic[ci]:(Scompar_indic[ci+1]-1)){
      # likelihood
      rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
      rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
      rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
      rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated
      # link and linear predictor
      cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
      cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]                # baseline in women
      cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         md[z]         # treatment in men
      cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + md[z] + wd[z] # treatment in women
      # Trial specific priors
      mu[z] ~ dnorm (0, 0.0001)
      wu[z] ~ dnorm (0, 0.0001)
      md[z]  ~ dnorm (0, 0.0001)
      wd[z] ~ dnorm(wd_delta[ci],prec1)
    }#trials
    wd_delta[ci] ~ dnorm(0, 0.0001)
  }#comparison and indication
  var1 ~  dlnorm(-3.23, 1/1.88^2) 
  prec1 <- 1 / var1
  sd1 <- var1^0.5
} # model
"
writeLines(modelstring, con = "jags_code/random_tx_strat.txt")

## Fixed effects model stratified by treatment "fixed_tx_strat.txt" ----
# comparison and indication flattened  - ci (comparison and indication in wide_rag_flat)
# study/trial - z (observations in wide_rag_flat)
# Stratified by treatemnt/comparison combination, NOT pooled across these
modelstring <- "
model{
  for(ci in 1:Ncompar_indic) {
    for(z in Scompar_indic[ci]:(Scompar_indic[ci+1]-1)){
      # likelihood
      rmc[z] ~ dbin(pmc[z],nmc[z]) # men controls
      rmt[z] ~ dbin(pmt[z],nmt[z]) # men treated
      rwc[z] ~ dbin(pwc[z],nwc[z]) # women controls
      rwt[z] ~ dbin(pwt[z],nwt[z]) # women treated
      # link and linear predictor
      cloglog(pmc[z]) <- log(time[z]) +  mu[z]                        # baseline in men
      cloglog(pwc[z]) <- log(time[z]) +  mu[z] + wu[z]                # baseline in women
      cloglog(pmt[z]) <- log(time[z]) +  mu[z] +         md[z]         # treatment in men
      cloglog(pwt[z]) <- log(time[z]) +  mu[z] + wu[z] + md[z] + wd_delta[ci] # treatment in women
      # Trial specific priors
      mu[z] ~ dnorm (0, 0.0001)
      wu[z] ~ dnorm (0, 0.0001)
      md[z]  ~ dnorm (0, 0.0001)
    }#trials
    wd_delta[ci] ~ dnorm(0, 0.0001)
  }#comparison and indication
  # ancestor priors
} # model
"
writeLines(modelstring, con = "jags_code/fixed_tx_strat.txt")


## List of files which require different data structures ----
# Models without treatment or indication levels
trials_only <- c("random_effects.txt", "stratified.txt", "fixed_effects.txt")
# Prasugrel alone women
pras_women_men <- "prasugrel_only_women_men.txt"
# Prasugrel alone men
pras_men <- "prasugrel_only_men.txt"
# Prasugrel alone men and bleeding
pras_men_bleed <- "prasugrel_only_men_bleed.txt"
# Bleeding random effects and stratified models
bleed_re_strat <- c("bleed_strat.txt", "bleed_re.txt")
rm(modelstring)
# Lau model
lau_choose <- "lau.txt"
# Lau SA
lau_sa_choose <- "lau_sa.txt"
# Models where treatment and condition is flattened into a single variable
tx_strat <- c("random_tx_strat.txt", "fixed_tx_strat.txt")
