# This script processes the cmdstan output from the cluster for the full model. 
# It needs to be run on the cluster since the output from each chain is so 
# large. The output is already thinned. This script will:
# - Read in each chain
# - Calculate HPDIs, mean, median, sd for each parameter
# - Calculate Geweke and Gelman diagnostics
# - Store summarized output as a csv for each model

Packages <- c("coda", "bayesplot", "doSNOW", "data.table")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

chn.dir <- "full_pieces/"
sum.dir <- "summaries/"

write("", "logs/procFull_R.output")
p.c <- makeCluster(2); registerDoSNOW(p.c)
foreach(p=0:15) %dopar% {
  Packages <- c("coda", "data.table")
  suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
  out.dt <- vector("list", 24)
  for(m in 1:24) {
    m.nm <- paste0("full_Clim_", m, "_", p)
    write(paste("Chain", m, "piece", p), "logs/procFull_R.output", append=T)
    out.dt[[m]] <- as.mcmc(fread(paste0(chn.dir, m.nm, ".csv")))
  }
  write("--Finished reading in chains", "logs/procFull_R.output", append=T)
  out <- as.mcmc.list(out.dt)
  rm(out.dt)
  out.all <- as.mcmc(do.call(rbind, out))
  
  # calculate HPD intervals
  write("HPD intervals", "logs/procFull_R.output", append=T)
  probs <- c(0.5, 0.75, 0.9, 0.95, 0.99)
  hpd.ls <- lapply(probs, function(x) HPDinterval(out.all, prob=x))
  for(i in 1:length(probs)) {
    colnames(hpd.ls[[i]]) <- paste(c("HPD_lo", "HPD_hi"), probs[i], sep="_")
  }
  qmn.m <- summary(out.all)
  summary.m <- do.call(cbind, hpd.ls)
  summary.m <- cbind(summary.m, qmn.m$statistics)
  summary.m <- cbind(summary.m, qmn.m$quantiles)
  rm(qmn.m); rm(hpd.ls)
  
  # diagnostics
  write("Running diagnostics", "logs/procFull_R.output", append=T)
  beta.ind <- grepl("beta", varnames(out))
  if(sum(beta.ind) > 0) {
    out.beta <- out[, beta.ind]
    geweke.m <- geweke.diag(out.beta)
    gelman.m <- gelman.diag(out.beta)
    write("Storing beta output", "logs/procFull_R.output", append=T)
    saveRDS(out.beta, paste0(sum.dir, "out_betas_FULL_", p, ".rds"))
    saveRDS(geweke.m, paste0(sum.dir, "geweke_FULL_", p, ".rds"))
    saveRDS(gelman.m, paste0(sum.dir, "gelman_FULL_", p, ".rds"))
  }
  
  # save output
  write("Storing output", "logs/procFull_R.output", append=T)
  write.csv(out.all, paste0(sum.dir, "out_all_FULL_", p, ".csv"))
  write.csv(summary.m, paste0(sum.dir, "summary_FULL_", p, ".csv"))
  write(paste("|||--Finished piece", p), "logs/procFull_R.output", append=T)
}
stopCluster(p.c)
