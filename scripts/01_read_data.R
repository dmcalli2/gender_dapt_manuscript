# 01_read_data
## Reads data and reshapes into various arrays for JAGS analysis
## Later in program selects which trials to include for main and sensitivity analyses
## based on whether data_choose is specified

# Packages ----
library(dplyr)
library(reshape2)
library(glmmBUGS)

# Source functions (for all the files with code in this project) ----
source ("scripts/00_functions.r")

# Read in csv files with data ----
main <- read.csv (file = "data/P2Y12meta_combined.csv", as.is = TRUE)
names(main) <- tolower(names(main))
names(main) <- gsub(".", "_", names(main), fixed = TRUE)

bleed <- read.csv(file = "data/P2Y12_bleeding.csv", as.is = TRUE)
names(bleed) <- tolower(names(bleed))
names(bleed) <- gsub(".", "_", names(bleed), fixed = TRUE)

# identify follow-up period ----
main$duration_of_follow_up
main$time <- c("12 months" = 12, "Hospital discharge or 28 days." = 1, "30 days" = 1, "15 months" = 15,  
                "12 months" = 12, "30 months" = 30, "90 days" = 3, "3.4 years" = 3.4*12, "90 days" = 3)/12

# Treat prasugrel and ticagrelor as single treatment  ----
# Collapse ticagrelor and prasugrel into a single comparison
main <- filter(main, indication %in% c("ACS", "Stroke")) %>%
    mutate (drug_intervention  = ifelse (drug_intervention == "ticagrelor", "prasugrel", drug_intervention))

# select study characteristics needed for analyses ----
study <- main [ , c('trial','indication', 'drug_control', 'drug_intervention', 'time')]
study$comparison <- paste0(study$drug_intervention, "_", study$drug_control)

# Collapse all stroke studies into a single comparison as there is only one ticagrelor versus aspirin study
study$comparison[study$indication == "Stroke"] <- "clopidogrel_placebo"

## Create a study ID 
study <- study[, c("trial","indication", "comparison", "time")]
study$study <- 1:nrow(study)

# Relabel events and n for main outcome
events_n <- main[, c(
  'women_intervention_event', 'women_intervention_no_event', 
  'women_control_event', 'women_control_no_event', 
  'men_intervention_event', 'men_intervention_no_event', 
  'men_control_event', 'men_control_no_event')]
events_n$study <- 1:nrow(events_n)
events_n <- melt(events_n, id = "study")
events_n$women <- grepl("women|female", events_n$variable )
events_n$tx <- grepl("intervention", events_n$variable )
events_n$event <- grepl("no_event", events_n$variable)
events_n$event <- factor (events_n$event, c(T,F), labels = c("no_event", "events"))
events_n <- dcast(events_n, study + women + tx ~ event, value.var = "value")
events_n$n <- events_n$events + events_n$no_event
events_n <- events_n[ , c("study","women", "tx", "events", "n")]

# Create a wide format dataset ----
wide <- events_n
wide$women <- factor(wide$women, levels = c(T,F), labels = c("w", "m"))
wide$tx <- factor(wide$tx, levels = c(T,F), labels = c("t", "c"))

events <- reshape2::dcast (wide, study ~ women + tx, value.var = "events" )
names(events) <- gsub("_", "", names(events))
names(events) <- paste0("r", names(events))

n <- reshape2::dcast (wide, study ~ women + tx, value.var = "n" )
names(n) <- gsub("_", "", names(n))
names(n) <- paste0("n", names(n))

main_wide <- cbind(study, n[,-1], events[,-1])

## Set-up ragged array with indication nested within treatment ----
main_wide <- main_wide %>%
  arrange(comparison, indication, trial) 

wide_rag <- winBugsRaggedArray(main_wide, 
                               effects = c("comparison", "indication"), 
                               covariates = list(trial = "trial"),
                               observations = "rwt",
                               returnData = FALSE) 


## Set-up ragged array with trial nested within a single treatment comparison/indication variable ----
# ie flatten the complexity
main_wide_flat <- main_wide %>%
  mutate(compar_indic = paste(comparison, indication, sep = "_")) %>%
  arrange(compar_indic, trial)

wide_rag_flat <- winBugsRaggedArray(main_wide_flat, 
                                      effects = c("compar_indic"), 
                                      covariates = list(trial = "trial"),
                                      observations = "rwt",
                                      returnData = FALSE) # checked ordering it is unchanged

## Select trials with just prasugrel versus clopidogrel ----
pras <- filter(main_wide, comparison == "prasugrel_clopidogrel" &
                 indication == "ACS")

## Set-up bleeding data for treatment comparison within indication ----
bleed <- merge(bleed,
                main_wide_flat [, ! names(main_wide_flat) %in% c('nwt', 'nwc', 'nmt', 'nmc', 'rwt', 'rwc', 'rmt', 'rmc')],
                by = "trial")
bleed$study <- 1:nrow(bleed)
bleed$rwt <- bleed$women_intervention_bleed ; bleed$women_intervention_bleed <- NULL
bleed$rwc <- bleed$women_control_bleed ; bleed$women_control_bleed <- NULL
bleed$rmt <- bleed$men_intervention_bleed ; bleed$men_intervention_bleed<- NULL
bleed$rmc <- bleed$men_control_bleed ; bleed$men_control_bleed <- NULL

bleed$xwt <- bleed$women_intervention_no_bleed ; bleed$women_intervention_no_bleed <- NULL
bleed$xwc <- bleed$women_control_no_bleed ; bleed$women_control_no_bleed <- NULL
bleed$xmt <- bleed$men_intervention_no_bleed ; bleed$men_intervention_no_bleed<- NULL
bleed$xmc <- bleed$men_control_no_bleed ; bleed$men_control_no_bleed <- NULL

bleed$nwt <- bleed$xwt + bleed$rwt
bleed$nwc <- bleed$xwc + bleed$rwc
bleed$nmt <- bleed$xmt + bleed$rmt
bleed$nmc <- bleed$xmc + bleed$rmc
bleed <- bleed [ , ! names(bleed) %in% c("xwt", "xwc", "xmt", "xmc")]
bleed <- bleed %>%
  arrange(indication, comparison, trial) 

bleed_rag <- winBugsRaggedArray(bleed, 
                                      effects = c("indication", "comparison"), 
                                      covariates = list(trial = "trial"),
                                      observations = "rwt",
                                      returnData = FALSE) # checked ordering it is unchanged

## Set-up bleeding data for prasugrel_only comparison ----
bleed_pras_men <- filter (bleed, trial %in% pras$trial)

## Read in data from J Am Coll Cardiol. 2017 Mar 28;69(12):1549-1559. doi: 10.1016/j.jacc.2017.01.028.
lau <- read.csv(file = "data/lau_et_al.csv", as.is = TRUE)
lau_time <- read.csv(file = "data/lau_et_al_time.csv", as.is = TRUE)
lau <- merge(lau_time, lau, by = "trial")
lau$study <- seq_along(lau$trial)
rm(lau_time)

## Create dataset with data from both systematic reviews
lau_sa <- lau[lau$trial %in% setdiff(lau$trial, main_wide$trial),]
lau_sa <- bind_rows(main_wide, lau_sa)
lau_sa$study <- seq_along(lau_sa$trial)

## Save data ----
save.image(file = "data/Data.Rdata")
