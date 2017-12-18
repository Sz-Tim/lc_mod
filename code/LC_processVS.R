# This script processes the cmdstan output from the cluster for the variable
# selection models. It needs to be run on the cluster since the output from 
# each chain is so large. Looping through each model, it will:
# - Read in each chain
# - Thin each chain using an interval of 20
# - Use loo package to calculate WAIC
# - Calculate HPDIs, mean, median, sd for each parameter
# - Store summarized output as a csv for each model

Packages <- c("rstan", "coda", "tidyverse", "purrr", "loo", "doSNOW")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

chn.dir <- "lc_small/out"
sum.dir <- "lc_small/summaries"


foreach(m=1:7) %dopar% {
  f.m <- list.files("out", pattern=paste0("1pct_", names(d.ls)[1], "_"),
                    full.names=TRUE)
  out <-  f.m %>% read_stan_csv %>% As.mcmc.list
  
}



