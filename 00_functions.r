# 00_functions.r
### Functions used in analyses

ExtractChains <- function (LINE.out, node_name) {
  # Extract each chain
  chain1 = LINE.out[[1]]
  chain2 = LINE.out[[1]]
  # If only a single column convert to vector, otherwise extract node and then convert to vector
  if (ncol(chain1)>=2) {
    print(colnames(chain1))
    chain1 = as.vector(chain1[ , node_name])
    chain2 = as.vector(chain2[ , node_name])
  } else {
    chain1 = as.vector(chain1)
    chain2 = as.vector(chain2)
  }
  # data.frame(chain1 = chain1, chain2 = chain2)
  c(chain1, chain2)
}


# Function to extract results from models srtatifying by trial ----
ExtractStratified <- function (LINE.out, mydf = "bleed", ntrials = nrow(get(mydf))) {
  # Extract trial-specific estimates
  # Read in raw data
  source ("scripts/01_read_data.R", local = TRUE)
  
  inter_trial <- lapply(LINE.out, function (x) x [ , paste0("wd[", 1:ntrials, "]")])
  inter_trial <- do.call (rbind, inter_trial)
  men_trial <- lapply(LINE.out, function (x) x [ , paste0("d[", 1:ntrials, "]")])
  men_trial <- do.call (rbind, men_trial)
  women_trial <- inter_trial + men_trial
  
  list_trial <- list(inter_trial, men_trial, women_trial)
  list_trial <- lapply(list_trial, function(x) {
    colnames(x) <- get(mydf)$trial
    x}
  )
  rm(inter_trial, men_trial, women_trial)
  list_trial <- lapply(list_trial, function (x) round(exp(x),2) )
  list_trial_q <- lapply(list_trial, function (x) apply(x, 2, quantile, probs = c(0.5, 0.025, 0.975)))
  list_trial_q <- do.call (rbind, list_trial_q)
  list_trial_q <- as.data.frame(list_trial_q)
  row.names(list_trial_q) <- paste (rep(c("inter", "men", "women"), each = 3), rep(c("q50", "q2.5", "q97.5"),3),
                                    sep = "_")
  list_trial_q$measure <- rep(c("q50", "q2.5", "q97.5"),3)
  list_trial_q$param <- rep(c("inter", "men", "women"), each = 3)
  list_trial_q <- reshape2::melt (list_trial_q, id = c("measure", "param"))
  reshape2::dcast (list_trial_q, param + variable ~ measure, value.var = "value")
}

## Function to fit non-linear models to posterior distribution in order to summarise these
Tmod <- function (myvectorchar, df_choose = 3){
  print(myvectorchar)
  myvector <- get(myvectorchar)
  myvector_m <- mean(myvector)
  myvector_s <- sd(myvector)
  a <- fitdistr(myvector, "t", start = list(m=myvector_m,s=myvector_s, df = df_choose),
                lower=c(myvector_m-1*myvector_s, 0.001,1))

  ## Sample from T distribution
  Rrt <- function (m = m, s = s, df = df) {
    rt(1000, df=df)*s + m
  }
  mydt <- Rrt(a$estimate["m"], a$estimate["s"], a$estimate["df"])
  
  plot(density(myvector), main = "", xlab = myvectorchar, lty = "dashed")
  lines(density(mydt), col = "blue", lty = "dashed")

  b <- round(a$estimate,3)
  title(main = paste0(myvectorchar, "\n","t(", "mu = ", b["m"], ", sd = ", b["s"], ", df = ", b["df"],")" ))
  a <- a$estimate
  a["prec"] <- 1/ a["s"]^2
  a <- a [ c("m", "prec", "df")]
  names(a) <- paste(myvectorchar, names(a), sep = "_")
  a
}

# Function to combine two matrices into a single matrix element by element
SmryPrcnt <-  function (x, y, fun, ...){
  # Takes two matrices and combines these
  # Not function applied to eahc element of vecor
 if(any(dim(x) != dim(y))) stop("Error different number of dimensions in matrix")
 store_dim <- dim(x)
 x <- as.vector(x)
 y <- as.vector(y)
 z <- round(100* x/y,1)
 z <- ifelse(!is.na(x),  paste0(x, " (", z, "%)"), "")
 matrix(z, nrow = store_dim[1], ncol = store_dim[2])
 }

## Make density plots
# Make dataframe
MakeDens <- function(mylist){
  mylist <- map(mylist, density)
  mylist <- map(mylist, ~ .x[c("x", "y")] %>%  bind_cols)
  mydf <- bind_rows(mylist, .id = "model_name")
  mydf
}

# Make plot
MakeDensPlot <- function(mydf){
ggplot(mydf, aes(x = x, y = y, colour = model_name, linetype = model_name)) +
  geom_line()
}


