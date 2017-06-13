# 09_Run_model_diagnostics
library(rjags)
library(coda)
library(stringr)

# Create directory
if(!dir.exists("figures/diagnostic")) dir.create ("figures/diagnostic")


## First information criteria
# Load all dic results into a single list
all_dic <- list.files("jags_samples_main", patt = "_dic")
all_dic_lst <- as.list(all_dic)
names(all_dic_lst) <- all_dic

for (i in all_dic){
  load(file = paste0("jags_samples_main/", i))
  all_dic_lst [[i]] <- mydic
}
rm(mydic)

# Extract diagnostics result
all_dic_lst2 <- lapply(all_dic_lst, function(x) {
  deviance = sum(x$deviance)
  penalty = sum(x$penalty)
  dic = sum(x$deviance + x$penalty)
  output <- c(deviance, penalty, dic)
  names(output) <- c("Deviance", "pD", "DIC")
  output
}
)
names(all_dic_lst2) <- gsub("_dic.Rdata", "", names(all_dic_lst2))
all_dic <- do.call(rbind, all_dic_lst2)

mace_dic <- all_dic [-grep("age|bleed|prasugrel", rownames(all_dic)), ]
mace_dic <- mace_dic [ c("random_effects",
                         "fixed_effects",
                         "random_tx_strat",
                         "fixed_tx_strat"),]
mace_dic <- as.data.frame(mace_dic)
mace_dic[] <- lapply(mace_dic, round)


bleed_dic <-  all_dic [grep("bleed", rownames(all_dic)), ]
bleed_dic <- round(bleed_dic,1)
bleed_dic <- as.data.frame (bleed_dic)

save(mace_dic, bleed_dic, file = "data/mace_dic.Rdata")

## Diagnostic plots
filenames <- list.files("jags_samples_main/", patt = "dgnstc_")
names(filenames) <- stringr::str_replace(filenames, "\\.Rdata", "")
for (i in seq_along(filenames)){
  load(paste0("jags_samples_main/", filenames[i]))
  print(sapply(list(niter, nchain, nvar), function (x) x(LINE.out) ))
  pdf(paste0("figures/diagnostic/", names(filenames)[i], ".pdf"))
  par(mfrow = c(4,3))
  traceplot(LINE.out)
  par(mfrow = c(1,1))
  autocorr.plot(LINE.out)
  gelman.plot(LINE.out, ylim = c(0.95, 1.2))
  dev.off()
}

