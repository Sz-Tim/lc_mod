



########
## Run variable selection models
########
# set parameters
res <- c("200ha", "2000ha")[2]
# n.core <- 7  # for running *models* in parallel in addition to chains
# n.chain <- 2
# n.iter <- 1000
# n.thin <- 1
# n.refr <- 50
# exclude <- c("PSI_", "Z_", "Z_new_")

# # load workspace
# library(tidyverse)
# v <- str_remove(list.files("data_stan", paste0(res, "_varSel")), ".Rdump")

# run models
library(doSNOW)
# p.c <- makeCluster(n.core); registerDoSNOW(p.c)
# foreach(i=seq_along(v), .packages="rstan") %dopar% {
#   options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
#   # nonspatial model
#   mod <- paste0(v[i], "_nonspatial_varSel")
#   out <- stan(file="ms/nonspatial_varSel.stan",
#               data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
#               iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
#               pars=exclude, include=F,
#               sample_file=paste0("out_stan/", mod, ".csv"))
#   
#   # spatial model
#   mod <- paste0(v[i], "_spatial_varSel")
#   out <- stan(file="ms/nonspatial_varSel.stan",
#               data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
#               iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
#               pars=exclude, include=F,
#               sample_file=paste0("out_stan/", mod, ".csv"))
# }
# stopCluster(p.c)






########
## Calculate lppd
########
# setup workspace
library(tidyverse); library(rstan); library(loo)
v <- str_remove(list.files("data_stan", paste0(res, "_varSel")), ".Rdump")
grnt.df <- read.csv(paste0("data/L_varSel_", res, ".csv")) %>%
  select(lon, lat, Train, CellID, grnt_Opn, grnt_Oth, grnt_Dec,
         grnt_WP, grnt_Evg, grnt_Mxd) %>%
  gather(LC, grnt, 5:10) %>%
  arrange(CellID)

lpd.non <- PSI.sum.non <- waic.non <- setNames(vector("list", length(v)), 
                                               paste0("non", v))
lpd.car <- PSI.sum.car <- waic.car <- setNames(vector("list", length(v)), 
                                               paste0("car", v))

for(i in seq_along(v)) {
  stan_d <- read_rdump(paste0("data_stan/", v[i], ".Rdump"))
  non.f <- paste0(v[i], "_nonspatial_varSel")
  car.f <- paste0(v[i], "_spatial_varSel")
  
  # load output
  out.non <- read_stan_csv(list.files("out_stan", paste0("^", non.f), full.names=T))
  out.car <- read_stan_csv(list.files("out_stan", paste0("^", car.f), full.names=T))
  
  # extract pointwise predictive density
  lppd.non <- extract_log_lik(out.non)
  lpd.non[[i]] <- -sum(colMeans(exp(lppd.non)))/(ncol(lppd.non))
  waic.non[[i]] <- waic(lppd.non[[i]])
  lppd.car <- extract_log_lik(out.car)
  lpd.car[[i]] <- -sum(colMeans(exp(lppd.car)))/(ncol(lppd.car))
  waic.car[[i]] <- waic(lppd.car[[i]])
  
  # extract eta summaries
  PSI.sum.non[[i]] <- summary(out.non, pars="PSI")$summary[,c(1,4,8)]
  PSI.sum.car[[i]] <- summary(out.car, pars="PSI")$summary[,c(1,4,8)]
  colnames(PSI.sum.non[[i]]) <- colnames(PSI.sum.car[[i]]) <- paste0("PSI_", 
                                                        c("mn", "025", "975"))
  cat("\n\nFinished", i, "of", length(v), "\n")
  waic.non[[i]]
  lpd.non[[i]]
  waic.car[[i]]
  lpd.car[[i]]
}

saveRDS(waic.non, paste0("out/waic_non_", res, ".rds"))
saveRDS(waic.car, paste0("out/waic_car_", res, ".rds"))
saveRDS(lppd.non, paste0("out/lppd_non_", res, ".rds"))
saveRDS(lppd.car, paste0("out/lppd_car_", res, ".rds"))
saveRDS(lpd.non, paste0("out/lpd_non_", res, ".rds"))
saveRDS(lpd.car, paste0("out/lpd_car_", res, ".rds"))

sort(unlist(lpd.non))  # ClimTopo, but models are very similar
sort(unlist(lpd.car))  # ClimTopo, but models are very similar

opt.non <- names(sort(unlist(lpd.non)))[1]
opt.car <- names(sort(unlist(lpd.car)))[1]

opt.non
opt.car




########
## Visualize
########
psi.df <- tibble(CellID=rep(grnt.df$CellID, 14),
                 lon=rep(grnt.df$lon, 14),
                 lat=rep(grnt.df$lat, 14),
                 Train=rep(grnt.df$Train, 14),
                 LC=rep(c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"), 
                        nrow(grnt.df)/6*14),
                 model=rep(str_split_fixed(v, "_", 2)[,1], each=nrow(grnt.df)),
                 grnt=rep(grnt.df$grnt, 14),
                 PSI.non=unlist(map(PSI.sum.non, ~.[,1])),
                 PSI.non.lo=unlist(map(PSI.sum.non, ~.[,2])),
                 PSI.non.hi=unlist(map(PSI.sum.non, ~.[,3])),
                 PSI.non.CI=PSI.non.hi-PSI.non.lo,
                 PSI.car=unlist(map(PSI.sum.car, ~.[,1])),
                 PSI.car.lo=unlist(map(PSI.sum.car, ~.[,2])),
                 PSI.car.hi=unlist(map(PSI.sum.car, ~.[,3])),
                 PSI.car.CI=PSI.car.hi-PSI.car.lo) %>%
  group_by(model, CellID) %>%
  mutate(PSI.non=PSI.non/sum(PSI.non),
         PSI.car=PSI.car/sum(PSI.car))
write.csv(psi.df, paste0("out/psi_df_", res, ".csv"))


# library(viridis); theme_set(theme_bw())
# psi.df <- read.csv(paste0("out/psi_df_", res, ".csv"))
# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=PSI.car)) + facet_grid(LC~model) +
#   scale_fill_viridis(limits=c(0,1))
# ggplot(filter(psi.df, !Train)) +
#   geom_tile(aes(lon, lat, fill=PSI.car.CI)) + facet_grid(LC~model) +
#   scale_fill_viridis(limits=c(0,1))
# ggplot(filter(psi.df, !Train)) +
#   geom_tile(aes(lon, lat, fill=PSI.non-grnt)) + facet_grid(LC~model) +
#   scale_fill_gradient2(limits=c(-1,1))
# ggplot(filter(psi.df, !Train)) + geom_density(aes(x=PSI.car-grnt, colour=model)) +
#   facet_wrap(~LC) + scale_colour_viridis(discrete=TRUE)
# ggplot(filter(psi.df, !Train)) +
#   geom_density(aes(x=PSI.car.CI), colour="black") +
#   geom_density(aes(x=PSI.non.CI), colour="red") +
#   facet_grid(LC~model, scales="free")
# 
# ggplot(filter(psi.df, !Train)) + xlim(0,1) + ylim(-1,1) +
#   geom_point(aes(x=grnt, y=PSI.non-grnt), colour="red", alpha=0.2, shape=1) +
#   geom_point(aes(x=grnt, y=PSI.car-grnt), colour="black", alpha=0.2, shape=1) +
#   geom_hline(yintercept=0, linetype=3) +
#   facet_grid(LC~model) + labs(x="GRANIT", y="PSI-GRANIT")
# ggplot(filter(psi.df, !Train)) + xlim(0,1) + ylim(0,1) +
#   geom_point(aes(x=grnt, y=PSI.non), colour="red", alpha=0.2, shape=1) +
#   geom_point(aes(x=grnt, y=PSI.car), colour="black", alpha=0.2, shape=1) +
#   geom_abline(linetype=3) + 
#   facet_grid(LC~model) + labs(x="GRANIT", y="PSI")
# 
# 
# psi.df %>% filter(!Train) %>% group_by(model) %>%
#   summarise(mnDiff.non=mean(abs(PSI.non-grnt)),
#             mnDiff.car=mean(abs(PSI.car-grnt)),
#             mdDiff.non=median(abs(PSI.non-grnt)),
#             mdDiff.car=median(abs(PSI.car-grnt)),
#             sumDiff.non=sum(abs(PSI.non-grnt)),
#             sumDiff.car=sum(abs(PSI.car-grnt))) %>%
#   arrange(sumDiff.car)

# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=PSI.car)) + facet_wrap(~LC) +
#   scale_fill_viridis(limits=c(0,1))
# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=PSI.CI)) + facet_wrap(~LC) +
#   scale_fill_viridis()
# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=PSI.car-grnt)) + facet_wrap(~LC) +
#   scale_fill_gradient2(limits=c(-1,1))
#
#
#
#
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(rstan); library(tidyverse); library(viridis)
# out <- stan(file="ms/nonspatial.stan",
#             data=read_rdump(paste0("data_stan/Cens_2000ha_varSel.Rdump")),
#             iter=1000, chains=1, refresh=10)
# 
# psi.df <- read.csv("data/L_varSel_2000ha.csv") %>% 
#   select(lon, lat, Train, CellID, grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg, grnt_Mxd) %>% 
#   gather(LC, grnt, 5:10) %>% 
#   arrange(CellID) %>% 
#   mutate(psi=c(get_posterior_mean(out, pars="PSI")[,1]),
#          psi_025=summary(out, pars="PSI")$summary[,4],
#          psi_975=summary(out, pars="PSI")$summary[,8],
#          psi_CI=psi_975-psi_025)
# 
# psiError.df <- read.csv("data/L_varSel_2000ha.csv") %>% 
#   select(lon, lat, Train, CellID, grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg, grnt_Mxd) %>% 
#   filter(!Train) %>% 
#   gather(LC, grnt, 5:10) %>% 
#   arrange(CellID) %>% 
#   mutate(error=c(get_posterior_mean(out, pars="PSI_Y_error")[,3]))
# 
# phi.df <- read.csv("data/L_varSel_2000ha.csv") %>% 
#   select(lon, lat, Train, CellID, grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg) %>% 
#   # filter(Train) %>% 
#   gather(LC, grnt, 5:9) %>% 
#   arrange(CellID) %>% 
#   mutate(phi=c(get_posterior_mean(out, pars="phi")[,1]),
#          phi_025=summary(out, pars="phi")$summary[,4],
#          phi_975=summary(out, pars="phi")$summary[,8],
#          phi_CI=phi_975-phi_025)
# 
# phiNew.df <- read.csv("data/L_varSel_2000ha.csv") %>% 
#   select(lon, lat, Train, CellID, grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg) %>% 
#   # filter(!Train) %>% 
#   gather(LC, grnt, 5:9) %>% 
#   arrange(CellID) %>% 
#   mutate(phi=c(get_posterior_mean(out, pars="new_phi")[,1]),
#          phi_025=summary(out, pars="new_phi")$summary[,4],
#          phi_975=summary(out, pars="new_phi")$summary[,8],
#          phi_CI=phi_975-phi_025)
# phi.df <- rbind(phi.df, phiNew.df)
# 
# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=psi)) + facet_wrap(~LC) + 
#   scale_fill_viridis(limits=c(0,1))
# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=psi_CI)) + facet_wrap(~LC) + 
#   scale_fill_viridis(limits=c(0,1))
# 
# ggplot(phi.df) + geom_tile(aes(lon, lat, fill=phi)) + facet_wrap(~LC)
# ggplot(phi.df) + geom_tile(aes(lon, lat, fill=phi_CI)) + facet_wrap(~LC)
# 
# ggplot(phiNew.df) + geom_tile(aes(lon, lat, fill=phi)) + facet_wrap(~LC)
# ggplot(phiNew.df) + geom_tile(aes(lon, lat, fill=phi_CI)) + facet_wrap(~LC)
# 
# 
# ggplot(psi.df) + geom_tile(aes(lon, lat, fill=error)) + 
#   scale_fill_gradient2(limits=c(-1,1)) + facet_wrap(~LC)
