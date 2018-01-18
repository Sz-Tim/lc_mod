Packages <- c("rstan", "coda", "bayesplot", "tidyverse", "rarhsmm",
              "loo", "doSNOW", "sevcheck", "data.table")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

n.cores <- 2
thin.seq <- seq(4,20,2)
chn.dir <- "out/"
sum.dir <- "thin_summaries/"
mods <- c("Cens", "Clim", "ClimCens", "Topo", 
          "TopoCens", "TopoClim", "TopoClimCens")
write("", "thin_proc.output")

p.c <- makeCluster(n.cores); registerDoSNOW(p.c)
foreach(m=1:7) %dopar% {
  #for(m in 1:7) {
  Packages <- c("rstan", "coda", "bayesplot", "tidyverse", "rarhsmm",
                "loo", "doSNOW", "sevcheck", "data.table")
  suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
  
  # read in chains
  write(paste(m, mods[m], "Reading in chains"), "thin_proc.output", append=T)
  f.m <- list.files(chn.dir, 
                    pattern=paste0("_", mods[m], "_"),
                    full.names=TRUE)
  # read in chains
  out.dt <- vector("list", length(f.m))
  for(i in seq_along(f.m)) {
    write(paste("", m, mods[m], "chain", i), "thin_proc.output", append=T)
    v.nm <- fread(paste("sed '/^#/ d'", f.m[i]), sep=",", nrows=0)
    col.exc <- grep(paste(c("nu", "Y2"), collapse="|"), names(v.nm))
    out.dt[[i]] <- as.mcmc(fread(paste("sed '/^#/ d'", f.m[i]), drop=col.exc))
  }
  out <- as.mcmc.list(out.dt)
  rm(out.dt)
  
  # thin by n.thin
  for(i in seq_along(thin.seq)) {
    n.thin <- thin.seq[i]
    write(paste(m, mods[m], "Thinning", n.thin), "thin_proc.output", append=T)
    n.it <- niter(out)
    th.it.chn <- (1:n.it)[(1:n.it %% n.thin) == 0]
    out.thin <- as.mcmc.list(lapply(out, function(x) mcmc(x[th.it.chn,])))
    out.all <- as.mcmc(do.call(rbind, out.thin))
    
    
    # diagnostics
    write(paste(m, mods[m], "Running diagnostics"), "thin_proc.output", append=T)
    out.beta <- out.thin[, grepl("beta", varnames(out.thin))]
    geweke.m <- geweke.diag(out.beta)
    gelman.m <- gelman.diag(out.beta)
    
    # save output
    write(paste(m, mods[m], "Storing output"), "thin_proc.output", append=T)
    write.csv(out.all, paste0(sum.dir, "out_all_", n.thin, "_", mods[m], ".csv"))
    saveRDS(out.thin, paste0(sum.dir, "out_thin_chains_", n.thin, "_", 
                             mods[m], ".rds"))
    saveRDS(out.beta, paste0(sum.dir, "out_betas_", n.thin, "_", 
                             mods[m], ".rds"))
    saveRDS(geweke.m, paste0(sum.dir, "geweke_", n.thin, "_", 
                             mods[m], ".rds"))
    saveRDS(gelman.m, paste0(sum.dir, "gelman_", 
                             n.thin, "_", mods[m], ".rds"))
    rm(out.all); rm(out.thin); rm(out.beta); rm(geweke.m)rm(gelman.m)
  }
}
stopCluster(p.c)


