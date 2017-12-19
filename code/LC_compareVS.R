# This script analyses the results of LC_processVS_hpc.R

Packages <- c("rstan", "coda", "ggmcmc", "bayesplot", "tidyverse", 
              "loo", "purrr")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

n.cores <- 4
chn.dir <- "lc_small/out/"
sum.dir <- "lc_small/summaries/"
mods <- c("Cens", "Clim", "ClimCens", "Topo", 
          "TopoCens", "TopoClim", "TopoClimCens")
out.type <- c("gelman", "geweke", "loo", "out_all", "out_betas",
              "out_thin_chains", "rhats_all", "rhats_beta", "summary", "waic")



# Rhats
rhat.beta.ls <- map(list.files(sum.dir, pattern="rhat_beta", full.names=TRUE),
                   ~read.csv(., col.names=c("Parameter", "Rhat")))
walk(rhat.beta.ls, ~hist(.$Rhat))
rhat.all.ls <- map(list.files(sum.dir, pattern="rhat_all", full.names=TRUE),
                    ~read.csv(., col.names=c("Parameter", "Rhat")))
walk(rhat.all.ls, ~hist(.$Rhat))


# WAIC
waic.ls <- map(list.files(sum.dir, pattern="waic", full.names=TRUE),
               readRDS) 
names(waic.ls) <- mods
map_dfr(waic.ls, ~(.$waic))


# LOO
loo.ls <- map(list.files(sum.dir, pattern="loo", full.names=TRUE),
               readRDS) 
names(loo.ls) <- mods
map_dfr(loo.ls, ~(.$looic))


# Gelman
par(mfrow=c(3,3))
map(list.files(sum.dir, pattern="gelman", full.names=TRUE), readRDS) %>% 
  walk2(.x=., .y=mods, ~(hist(x=.$psrf[,2], main=.y, xlim=c(0.95, 2))))
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

