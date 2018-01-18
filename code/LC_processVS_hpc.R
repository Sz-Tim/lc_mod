# This script processes the cmdstan output from the cluster for the variable
# selection models. It needs to be run on the cluster since the output from 
# each chain is so large. Looping through each model, it will:
# - Read in each chain
# - Thin each chain using an interval of 20
# - Use loo package to calculate WAIC
# - Calculate HPDIs, mean, median, sd for each parameter
# - Store summarized output as a csv for each model

Packages <- c("rstan", "coda", "bayesplot", "tidyverse", "rarhsmm",
              "loo", "doSNOW", "sevcheck", "data.table")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
tr_gjam_inv <- function(w, a=0.99) {
  eta <- w[-length(w)]
  eta[eta<0] <- 0
  w.p <- sum( (eta > 0 & eta < 1)*eta + (eta > 1) )
  
  if(w.p >= a) {
    while(sum(eta) > 1) {
      D.i <- (w.p^(-1)) * (1 - ((1-a)^(w.p/a)))
      eta <- D.i*eta
    }
  }
  c(eta, 1-sum(eta)) 
}

n.cores <- 2
n.thin <- 20
chn.dir <- "out/"
sum.dir <- "summaries/"
new.dir <- "data/strat_15pct_85oos/"
mods <- c("Cens", "Clim", "ClimCens", "Topo", 
          "TopoCens", "TopoClim", "TopoClimCens")
write("", "LC_proc.output")

p.c <- makeCluster(n.cores); registerDoSNOW(p.c)
foreach(m=1:7) %dopar% {
#for(m in 1:7) {
  Packages <- c("rstan", "coda", "bayesplot", "tidyverse", "rarhsmm",
                "loo", "doSNOW", "sevcheck", "data.table")
  suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
  
  # read in chains
  write(paste(m, mods[m], "Reading in chains"), "LC_proc.output", append=T)
  f.m <- list.files(chn.dir, 
                    pattern=paste0("_", mods[m], "_"),
                    full.names=TRUE)
  # read in chains
  out.dt <- vector("list", length(f.m))
  for(i in seq_along(f.m)) {
    write(paste("", m, mods[m], "chain", i), "LC_proc.output", append=T)
    v.nm <- fread(paste("sed '/^#/ d'", f.m[i]), sep=",", nrows=0)
    col.exc <- grep(paste(c("nu", "Y2"), collapse="|"), names(v.nm))
    out.dt[[i]] <- as.mcmc(fread(paste("sed '/^#/ d'", f.m[i]), drop=col.exc))
  }
  out <- as.mcmc.list(out.dt)
  rm(out.dt)
  
  # thin by n.thin
  write(paste(m, mods[m], "Thinning"), "LC_proc.output", append=T)
  n.it <- niter(out)
  th.it.chn <- (1:n.it)[(1:n.it %% n.thin) == 0]
  out <- as.mcmc.list(lapply(out, function(x) mcmc(x[th.it.chn,])))
  out.all <- as.mcmc(do.call(rbind, out))
  
  # calculate out of sample predictive error
  write(paste(m, mods[m], "OOS validation"), "LC_proc.output", append=T)
  ## set up
  oos.d <- read_rdump(paste0(new.dir, mods[m], ".Rdump"))
  di <- oos.d$di
  nX_d <- di[1]:di[2]
  new.Y1 <- oos.d$Y1
  new.Y1 <- cbind(new.Y1, 1-rowSums(new.Y1))
  beta_d <- out.all[, grepl("beta_d", varnames(out.all))]
  beta_p <- out.all[, grepl("beta_p", varnames(out.all))]
  cov.Y2 <- out.all[, grepl("L_Sigma.2", varnames(out.all), fixed=TRUE)]
  new.pred <- array(dim=c(nrow(oos.d$Y2), ncol(oos.d$Y2)+2, nrow(beta_d)))
  lppd <- array(dim=c(nrow(oos.d$Y2), ncol(oos.d$Y2)+1, nrow(beta_d)))
  for(i in 1:dim(new.pred)[3]) { # iterations
    ## generate predictions
    new.pred[,1,i] <- oos.d$Y2[,1] + oos.d$X[,nX_d] %*% beta_d[i,di[1]:di[2]]
    new.pred[,2,i] <- oos.d$Y2[,2] + oos.d$X[,nX_d] %*% beta_d[i,di[3]:di[4]]
    new.pred[,3,i] <- oos.d$Y2[,3] + oos.d$X[,nX_d] %*% beta_d[i,di[5]:di[6]]
    new.pred[,4,i] <- (oos.d$Y2[,4] + oos.d$X[,nX_d] %*% beta_d[i,di[7]:di[8]])*
      antilogit(oos.d$X %*% beta_p[i,])
    new.pred[,5,i] <- (oos.d$Y2[,4] + oos.d$X[,nX_d] %*% beta_d[i,di[7]:di[8]])*
      (1-antilogit(oos.d$X %*% beta_p[i,]))
    new.pred[,,i] <- t(apply(new.pred[,,i], 1, tr_gjam_inv))
    ## calculate pointwise predictive density
    cov.mx <- matrix(cov.Y2[i,], nrow=5)
    for(j in 1:dim(new.pred)[1]) { # pixels
      lppd[j,,i] <- mvdnorm(rbind(new.pred[j,-6,i]), new.Y1[j,-6], 
                            cov.mx, logd=FALSE)
    }
  }
  ## calculate oos scores
  lpd <- -sum(log(apply(lppd, 1:2, mean)))/prod(dim(new.Y1[,-6]))
  pred.mn <- apply(new.pred, 1:2, mean)
  MSPE <- sum((new.Y1 - pred.mn)^2/prod(dim(new.Y1)))
  rm(new.Y1); rm(beta_d); rm(beta_p); rm(oos.d)
  
  # calculate WAIC
  write(paste(m, mods[m], "WAIC & LOO"), "LC_proc.output", append=T)
  LL <- as.matrix(out.all[,grepl("log_lik", varnames(out.all))])
  n1 <- ncol(LL)/2
  LL <- LL[,1:n1] + LL[,n1+(1:n1)]
  waic.m <- waic(LL)
  loo.m <- loo(LL)
  rm(LL)
  
  # calculate HPD intervals
  write(paste(m, mods[m], "HPD intervals"), "LC_proc.output", append=T)
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
  write(paste(m, mods[m], "Running diagnostics"), "LC_proc.output", append=T)
  out.beta <- out[, grepl("beta", varnames(out))]
  geweke.m <- geweke.diag(out.beta)
  gelman.m <- gelman.diag(out.beta)
  
  # save output
  write(paste(m, mods[m], "Storing output"), "LC_proc.output", append=T)
  write.csv(out.all, paste0(sum.dir, "out_all_", mods[m], ".csv"))
  saveRDS(out, paste0(sum.dir, "out_thin_chains_", mods[m], ".rds"))
  saveRDS(out.beta, paste0(sum.dir, "out_betas_", mods[m], ".rds"))
  write.csv(summary.m, paste0(sum.dir, "summary_", mods[m], ".csv"))
  saveRDS(geweke.m, paste0(sum.dir, "geweke_", mods[m], ".rds"))
  saveRDS(gelman.m, paste0(sum.dir, "gelman_", mods[m], ".rds"))
  saveRDS(waic.m, paste0(sum.dir, "waic_", mods[m], ".rds"))
  saveRDS(loo.m, paste0(sum.dir, "loo_", mods[m], ".rds"))
  saveRDS(MSPE, paste0(sum.dir, "MSPE_", mods[m], ".rds"))
  saveRDS(lpd, paste0(sum.dir, "lpd_", mods[m], ".rds"))
  saveRDS(new.pred, paste0(sum.dir, "new_pred_", mods[m], ".rds"))
  saveRDS(pred.mn, paste0(sum.dir, "pred_mn_", mods[m], ".rds"))
  rm(out.all); rm(out); rm(out.beta); rm(summary.m); rm(geweke.m)
  rm(gelman.m); rm(waic.m); rm(loo.m); rm(MSPE); rm(new.pred); rm(pred.mn)
}
stopCluster(p.c)


