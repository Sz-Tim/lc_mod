# This script processes the cmdstan output from the cluster for the variable
# selection models. It needs to be run on the cluster since the output from 
# each chain is so large. Looping through each model, it will:
# - Read in each chain
# - Thin each chain using an interval of 20
# - Use loo package to calculate WAIC
# - Calculate HPDIs, mean, median, sd for each parameter
# - Store summarized output as a csv for each model

Packages <- c("rstan", "coda", "bayesplot", "tidyverse", "loo", "doSNOW")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

n.cores <- 2
n.thin <- 20
chn.dir <- "out/"
sum.dir <- "summaries/"
mods <- c("Cens", "Clim", "ClimCens", "Topo", 
          "TopoCens", "TopoClim", "TopoClimCens")

p.c <- makeCluster(n.cores); registerDoSNOW(p.c)
foreach(m=1:7) %dopar% {
  Packages <- c("rstan", "coda", "bayesplot", "tidyverse", "loo", "doSNOW")
  suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
  
  # read in chains
  f.m <- list.files(chn.dir, 
                    pattern=paste0("_", mods[m], "_"),
                    full.names=TRUE)
  out.stan <-  f.m %>% read_stan_csv
  out <-  out.stan %>% As.mcmc.list
  
  # thin by n.thin
  n.it <- niter(out)
  th.it.chn <- (1:n.it)[(1:n.it %% n.thin) == 0]
  out <- as.mcmc.list(lapply(out, function(x) mcmc(x[th.it.chn,])))
  #out.stan <- apply(rstan::extract(out.stan, permuted=F)[th.it.chn,,], 3, rbind)
  out.all <- as.mcmc(do.call(rbind, out))
  
  # calculate WAIC
  LL <- as.matrix(out.all[,grepl("log_lik", varnames(out.all))])
  waic.m <- waic(LL)
  loo.m <- loo(LL)
  
  # calculate HPD intervals
  probs <- c(0.5, 0.75, 0.9, 0.95, 0.99)
  hpd.ls <- lapply(probs, function(x) HPDinterval(out.all, prob=x))
  for(i in 1:length(probs)) {
    colnames(hpd.ls[[i]]) <- paste(c("HPD_lo", "HPD_hi"), probs[i], sep="_")
  }
  qmn.m <- summary(out.all)
  summary.m <- do.call(cbind, hpd.ls)
  summary.m <- cbind(summary.m, qmn.m$statistics)
  summary.m <- cbind(summary.m, qmn.m$quantiles)
  
  # diagnostics
  out.beta <- out[,grepl("beta", varnames(out))]
  geweke.m <- geweke.diag(out.beta)
  gelman.m <- gelman.diag(out.beta)
  rhats <- rhat(out.stan)
  rhats.beta <- rhat(out.stan, pars=c("beta_d", "beta_p"))
  
  # save output
  write.csv(out.all, paste0(sum.dir, "out_all_", mods[m], ".csv"))
  saveRDS(out, paste0(sum.dir, "out_thin_chains_", mods[m], ".rds"))
  saveRDS(out.beta, paste0(sum.dir, "out_betas_", mods[m], ".rds"))
  write.csv(summary.m, paste0(sum.dir, "summary_", mods[m], ".csv"))
  saveRDS(geweke.m, paste0(sum.dir, "geweke_", mods[m], ".rds"))
  saveRDS(gelman.m, paste0(sum.dir, "gelman_", mods[m], ".rds"))
  write.csv(rhats, paste0(sum.dir, "rhat_all_", mods[m], ".csv"))
  write.csv(rhats.beta, paste0(sum.dir, "rhat_beta_", mods[m], ".csv"))
  saveRDS(waic.m, paste0(sum.dir, "waic_", mods[m], ".rds"))
  saveRDS(loo.m, paste0(sum.dir, "loo_", mods[m], ".rds"))
  
}
stopCluster(p.c)


