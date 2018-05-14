# This script runs the variable selection version of the nonspatial and spatial
# models. Chains run in parallel if `options(mc.cores=parallel::detectCores())`
# is run first, with one chain on each available core. If you have excess cores
# (e.g., on a cluster), set n.core > 1.

########
## setup
########
n.core <- 2  # for running *models* in parallel in addition to chains
n.chain <- 2
n.iter <- 10
n.thin <- 10
n.refr <- 50
exclude <- c("nu", "Z_", "Z_new_")
opt.non <- c("ClimTopo")
opt.car <- c("CensClim")


########
## run models
########
library(rstan)
options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
out.non <- stan(file="R/nonspatial.stan",
                data=read_rdump(paste0("data_stan/", opt.non, ".Rdump")),
                iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
                pars=exclude, include=F,
                sample_file=paste0("out_stan/", opt.non, "_nonspatial.csv"))
out.car <- stan(file="R/spatial.stan",
                data=read_rdump(paste0("data_stan/", opt.car, ".Rdump")),
                iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
                pars=exclude, include=F,
                sample_file=paste0("out_stan/", opt.car, "_spatial.csv"))



########
## munge output
########
# load full model for optimal models
L <- read.csv("data/L_full.csv")
opt <- paste0(c(opt.non, opt.car), c("_nonspatial", "_spatial"))
full.non <- read_stan_csv(list.files("out_stan", paste0("^", opt[1]), full.names=T))
full.car <- read_stan_csv(list.files("out_stan", paste0("^", opt[2]), full.names=T))
full.df <- L %>% gather(LC, grnt, 4:9) %>%
  mutate(nlcd=with(L, c(nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Mxd,
                        nlcd_Evg*(1-pWP/100), nlcd_Evg*pWP/100))) %>%
  select(-(4:8))
full.df$LC <- factor(full.df$LC, 
                     levels=paste0("grnt_", 
                                   c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd")),
                     labels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
full.df <- full.df %>% arrange(desc(Fit), cellID, LC) %>% 
  mutate(non=c(get_posterior_mean(full.non, pars="eta")[,4]),
         car=c(get_posterior_mean(full.car, pars="eta")[,4])) %>%
  gather(Space, mn, 17:18) %>%
  cbind(., rbind(summary(full.non, pars="eta")$summary[,c(4:8)],
                 summary(full.car, pars="eta")$summary[,c(4:8)])) %>%
  rename(q025="2.5%", q25="25%", q50="50%", q75="75%", q975="97.5%")
write.csv(full.df, "out/3km_out.csv", row.names=F)










