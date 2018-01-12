# This script analyzes the results of LC_processVS_hpc.R

Packages <- c("rstan", "coda", "ggmcmc", "bayesplot", "tidyverse", 
              "loo", "purrr")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

setwd("C:/Users/tms1044/Desktop/lc_mod")

n.cores <- 4
chn.dir <- "out/"
sum.dir <- "summaries/"
mods <- c("Cens", "Clim", "ClimCens", "Topo", 
          "TopoCens", "TopoClim", "TopoClimCens")

# GGMCMC
beta.fls <- list.files(sum.dir, pattern="out_betas_", full.names=TRUE)
names(beta.fls) <- mods
map(beta.fls, readRDS) %>%
  map(., ~(ggs_density(ggs(.x)) + facet_wrap(~Parameter, scales="free")))
map(beta.fls, readRDS) %>%
  map(., ~(ggs_traceplot(ggs(.x)) + facet_wrap(~Parameter, scales="free")))
map(beta.fls, readRDS) %>% map(., ~(ggs_autocorrelation(ggs(.x))))
map(beta.fls, readRDS) %>% map(., ~(ggs_geweke(ggs(.x))))
map(beta.fls, readRDS) %>% map(., ~(ggs_Rhat(ggs(.x))))
map(beta.fls, readRDS) %>% 
  map(., ~(ggs_pairs(ggs(.x), lower=list(continuous="density"))))
map(beta.fls, readRDS) %>%
  map(., ~(ggs_running(ggs(.x)) + facet_wrap(~Parameter, scales="free")))
map(beta.fls, readRDS) %>%
  map(., ~(ggs_compare_partial(ggs(.x)) + facet_wrap(~Parameter, scales="free")))


# effective sample size
par(mfrow=c(3,3))
map(beta.fls, readRDS) %>% 
  walk(., ~hist(effectiveSize(.), main="", xlab="Effective Sample Size"))
par(mfrow=c(3,3))
map(beta.fls, readRDS) %>% 
  walk(., ~hist(effectiveSize(.)/(250*6), main="",
                xlab="n_eff/n_iter", xlim=c(0,1.5)))

# WAIC
waic.ls <- map(list.files(sum.dir, pattern="waic", full.names=TRUE), readRDS) 
names(waic.ls) <- mods
map_dbl(waic.ls, ~(.$waic)) %>% sort
compare(waic.ls[[1]], waic.ls[[2]], waic.ls[[3]], waic.ls[[4]], 
        waic.ls[[5]], waic.ls[[6]], waic.ls[[7]])


# LOO
loo.ls <- map(list.files(sum.dir, pattern="loo", full.names=TRUE), readRDS) 
names(loo.ls) <- mods
map_dbl(loo.ls, ~(.$looic)) %>% sort
compare(loo.ls[[1]], loo.ls[[2]], loo.ls[[3]], loo.ls[[4]], 
        loo.ls[[5]], loo.ls[[6]], loo.ls[[7]])


# OOS prediction
mspe.ls <- map(list.files(sum.dir, pattern="MSPE", full.names=TRUE), readRDS)
names(mspe.ls) <- mods
map_dbl(mspe.ls, sqrt) %>% sort


# Gelman
par(mfrow=c(3,3))
map(list.files(sum.dir, pattern="gelman", full.names=TRUE), readRDS) %>% 
  walk2(.x=., .y=mods, ~({hist(x=.$psrf[,2], main=.y, xlim=c(0.95, 1.2)); 
    abline(v=1.1, lty=3, col="red")}))
map(list.files(sum.dir, pattern="out_betas", full.names=TRUE), readRDS) %>% 
  walk2(.x=., .y=mods,  ~(gelman.plot(x=.x, ylab=paste(.y, "shrink factor"),
                                      ylim=c(0.95, 1.2))))


# Geweke
gew.ls <- map(list.files(sum.dir, pattern="geweke", full.names=TRUE), readRDS) 
for(m in 1:length(mods)) {
  par(mfrow=c(2,3))
  for(i in 1:length(gew.ls[[m]])) {
    plot(density(gew.ls[[m]][[i]]$z), main=mods[m], xlim=c(-5,5), ylim=c(0,0.5))
    curve(dnorm(x, 0, 1), from=-5, to=5, add=TRUE, col="dodgerblue")
    abline(v=c(-1.96, 1.96), lty=3)
  }
}
map(list.files(sum.dir, pattern="out_betas", full.names=TRUE),
    readRDS) %>% walk(geweke.plot)










