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
# load workspace
library(stringr)
v <- str_remove(list.files("data_stan", "_VS"), ".Rdump")


########
## run models
########
if(n.core > 1) {  # run models in parallel
  library(doSNOW); p.c <- makeCluster(n.core); registerDoSNOW(p.c)
  foreach(i=seq_along(v)) %dopar% {
    library(rstan)
    options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
    mod <- paste0(v[i], "_nonspatial_vs")
    out <- stan(file="R/nonspatial_vs.stan",
                data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
                iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
                pars=exclude, include=F,
                sample_file=paste0("out_stan/", mod, ".csv"))
  }
  foreach(i=seq_along(v)) %dopar% {
    library(rstan)
    options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
    mod <- paste0(v[i], "_spatial_vs")
    out <- stan(file="R/spatial_vs.stan",
                data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
                iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
                pars=exclude, include=F,
                sample_file=paste0("out_stan/", mod, ".csv"))
  }
  stopCluster(p.c)
} else {  # run models in serial
  library(rstan)
  options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
  for(i in seq_along(v)) {
    mod <- paste0(v[i], "_nonspatial_vs")
    out <- stan(file="R/nonspatial_vs.stan",
                data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
                iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
                pars=exclude, include=F,
                sample_file=paste0("out_stan/", mod, ".csv"))
  }
  for(i in seq_along(v)) {
    mod <- paste0(v[i], "_spatial_vs")
    out <- stan(file="R/spatial_vs.stan",
                data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
                iter=n.iter, chains=n.chain, refresh=n.refr, thin=n.thin,
                pars=exclude, include=F,
                sample_file=paste0("out_stan/", mod, ".csv"))
  }
}


########
## compare lppd
########
library(rarhsmm); library(tidyverse); library(rstan)

# storage structures
lpd.non <- lpd.car <- eta.sum.non <- eta.sum.car <- phi <- vector("list", length(v))
names(lpd.non) <- names(eta.sum.non) <- paste0("non_", v)
names(lpd.car) <- names(eta.sum.car) <- names(phi) <- paste0("car_", v)
for(i in seq_along(v)) {
  non <- paste0(v[i], "_nonspatial_vs")
  car <- paste0(v[i], "_spatial_vs")
  stan_d <- read_rdump(paste0("data_stan/", v[i], ".Rdump"))
  Y_new <- stan_d$Y_new
  
  # load output
  out.non <- read_stan_csv(list.files("out_stan", paste0("^", non), full.names=T))
  out.car <- read_stan_csv(list.files("out_stan", paste0("^", car), full.names=T))
  
  # extract predictions, covariate matrices
  eta.non <- rstan::extract(out.non, pars="eta")$eta[,stan_d$n2:stan_d$n3,-stan_d$D]
  eta.car <- rstan::extract(out.car, pars="eta")$eta[,stan_d$n2:stan_d$n3,-stan_d$D]
  Sigma.non <- rstan::extract(out.non, pars="Sigma")$Sigma[,2,,]
  Sigma.car <- rstan::extract(out.car, pars="Sigma")$Sigma[,2,,]
  n.j <- dim(eta.non)[1]
  n.cell <- dim(eta.non)[2]
  
  # calculate pointwise predective density
  lppd.non <- lppd.car <- matrix(0, nrow=n.j, ncol=n.cell)
  for(j in 1:n.j) {
    lppd.non[j,] <- sapply(1:n.cell, function(k) 
      mvdnorm(rbind(eta.non[j,k,]), Y_new[k,], Sigma.non[j,,], logd=F))
    lppd.car[j,] <- sapply(1:n.cell, function(k) 
      mvdnorm(rbind(eta.car[j,k,]), Y_new[k,], Sigma.car[j,,], logd=F))
  }
  
  # calculate log predictive density
  lpd.non[[i]] <- -sum(log(apply(lppd.non, 2, mean)))/prod(dim(Y_new[,-6]))
  lpd.car[[i]] <- -sum(log(apply(lppd.car, 2, mean)))/prod(dim(Y_new[,-6]))

  # extract eta summaries
  eta.sum.non[[i]] <- summary(out.non, pars="eta")$summary[,c(1,4,8)]
  colnames(eta.sum.non[[i]]) <- paste("non", 
                                      str_split_fixed(v[i], "_", 2)[,1],
                                      c("mn", "025", "975"), sep="_")
  eta.sum.car[[i]] <- summary(out.car, pars="eta")$summary[,c(1,4,8)]
  colnames(eta.sum.car[[i]]) <- paste("car", 
                                      str_split_fixed(v[i], "_", 2)[,1],
                                      c("mn", "025", "975"), sep="_")
  
  # extract spatial random effects
  phi[[i]] <- c(get_posterior_mean(out.car, pars="phi")[,4])
}

sort(unlist(lpd.non))  # ClimTopo, but models are very similar
sort(unlist(lpd.car))  # CensClim, but models are very similar

opt.non <- names(sort(unlist(lpd.non)))[1]
opt.car <- names(sort(unlist(lpd.car)))[1]

for(i in 1:length(v)) {
  colnames(eta.sum.non[[i]]) <- paste("non", 
                                      str_split_fixed(v[i], "_", 2)[,1],
                                      c("mn", "025", "975"), sep="_")
  colnames(eta.sum.car[[i]]) <- paste("car", 
                                      str_split_fixed(v[i], "_", 2)[,1],
                                      c("mn", "025", "975"), sep="_")
}


# munge output
L.YZ <- read.csv("data/L_YZ.csv")
comp.df <- L.YZ %>% gather(LC, grnt, 4:9) %>%
  mutate(nlcd=with(L.YZ, c(nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Mxd,
                           nlcd_Evg*(1-pWP/100), nlcd_Evg*pWP/100))) %>%
  select(-(4:8))
comp.df$LC <- factor(comp.df$LC, 
                     levels=paste0("grnt_", 
                                   c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd")),
                     labels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
comp.df <- comp.df %>% arrange(desc(Fit), cellID, LC) %>% 
  cbind(., do.call("cbind", eta.sum.non), do.call("cbind", eta.sum.car))
comp.long <- comp.df %>% select(-contains("025")) %>% select(-contains("975")) %>%
  gather(., model, mn, 17:30)  %>%
  cbind(., rbind(do.call("rbind", eta.sum.non)[,2:3],
                 do.call("rbind", eta.sum.car)[,2:3])) %>%
  mutate(Space=str_split_fixed(model, "_", 3)[,1],
         Cov=str_split_fixed(model, "_", 3)[,2]) %>%
  rename(q025=non_Cens_025, q975=non_Cens_975)
comp.space <- comp.long  %>% select(-model) %>% spread(Space, mn)
opt.long <- comp.long %>% filter(model %in% c("non_ClimTopo_mn", "car_CensClim_mn"))
phi.df <- comp.long %>% filter(LC != "Mxd" & Space=="car") %>%
  mutate(phi=unlist(phi))


# load full model for optimal models
L <- read.csv("data/L_full.csv")
opt <- paste0(str_split_fixed(c(opt.non, opt.car), "_", 3)[,2], 
              c("_nonspatial", "_spatial"))
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



# maps of LC proportions
prop.p <- ggplot(comp.df, aes(x=lon, y=lat)) + facet_wrap(~LC) + 
  scale_fill_gradient("", low="white", high="red", limits=c(0,1))
prop.p + geom_tile(aes(fill=grnt)) + ggtitle("GRANIT")
prop.p + geom_tile(aes(fill=nlcd)) + ggtitle("NLCD")
prop.p + geom_tile(aes(fill=Fit)) + scale_fill_manual(values=c("gray", "blue"))
prop.p + geom_tile(aes(fill=non_ClimTopo_mn)) + 
  ggtitle("Nonspatial optimal model mean")
prop.p + geom_tile(aes(fill=non_ClimTopo_975-non_ClimTopo_025)) + 
  ggtitle("Nonspatial optimal model 95% range")
prop.p + geom_tile(aes(fill=car_CensClim_mn)) +
  ggtitle("CAR optimal model mean")
prop.p + geom_tile(aes(fill=car_CensClim_975-car_CensClim_025)) +
  ggtitle("CAR optimal model 95% range")


# maps of residuals
comp.p <- ggplot(comp.df, aes(x=lon, y=lat)) + facet_wrap(~LC) + 
  scale_fill_gradient2("Residual\nerror", low="blue", mid="white",
                       high="red", limits=c(-1,1))
comp.p + geom_tile(aes(fill=non_ClimTopo_VS-grnt)) + ggtitle("No spatial random effects")
comp.p + geom_tile(aes(fill=car_CensClim_VS-grnt)) + ggtitle("CAR")


# table of RMSE
opt.long %>%  
# comp.long %>% 
  filter(!Fit) %>%
  group_by(Space, Cov, LC) %>% 
  summarise(RMSE_nlcd=sqrt(mean( (nlcd - grnt)^2 )),
            RMSE_mod=sqrt(mean( (mn - grnt)^2 ))) %>%
  mutate(pct_diff=(RMSE_mod - RMSE_nlcd)/RMSE_nlcd*100)
comp.long %>% filter(!Fit) %>% group_by(Space, Cov, LC) %>% 
  summarise(mn_95int=mean(q975-q025)) %>%
  ggplot(aes(x=Cov, y=mn_95int, colour=Space)) + geom_point() + 
  facet_wrap(~LC) + ylim(0,1) + 
  theme(axis.text.x=element_text(angle=300, hjust=0, vjust=1))


# density of residuals
ggplot(filter(comp.long, !Fit), aes(x=mn-grnt, colour=Space)) + 
  geom_vline(xintercept=0, linetype=3) + geom_density() + facet_grid(Cov~LC) +
  xlim(-.5,.5) + ggtitle("OOS prediction error")
ggplot(filter(comp.long, Space=="car" & !Fit), aes(x=mn-grnt, colour=Cov)) + 
  geom_vline(xintercept=0, linetype=3) + geom_density() + facet_wrap(~LC) +
  xlim(-.5,.5) + ggtitle("CAR: OOS prediction error")
ggplot(filter(comp.long, Space=="non" & !Fit), aes(x=mn-grnt, colour=Cov)) + 
  geom_vline(xintercept=0, linetype=3) + geom_density() + facet_wrap(~LC) +
  xlim(-.5,.5) + ggtitle("No spatial random effects: OOS prediction error")


# biplot of grnt, mean
ggplot(filter(comp.long, !Fit), aes(y=mn, x=grnt, colour=Space)) + 
  geom_abline() + geom_point(alpha=.25) + geom_rug(alpha=0.2) +
  facet_grid(LC~Cov) + xlim(0,1) + ylim(0,1) +
  labs(x="GRANIT", y="Posterior mean")
ggplot(filter(comp.long, !Fit), aes(y=mn, x=grnt, colour=Space)) + 
  geom_abline() + geom_density2d() + geom_rug(alpha=0.2) +
  stat_smooth(se=F, method="loess") + 
  facet_grid(LC~Cov) + xlim(0,1) + ylim(0,1) +
  labs(x="GRANIT", y="Posterior mean")


# map of spatial random effects
ggplot(phi.df, aes(x=lon, y=lat, fill=phi)) +
  geom_tile() + facet_grid(LC~Cov) + 
  scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-1,1))










# old

########
## set up workspace
########
library(tidyverse); library(rstan)
options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
data_f <- "data/landscape_3km.csv"
v <- str_remove(list.files("data_stan"), ".Rdump")
m <- c("nonspatial", "spatial")

L <- read.csv(data_f) %>% arrange(desc(Fit), cellID)
comp.df <- L %>% gather(LC, grnt, 4:9) %>%
  mutate(nlcd=with(L, c(nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Mxd,
                        nlcd_Evg*(1-pWP/100), nlcd_Evg*pWP/100))) %>%
  select(-(4:8))
comp.df$LC <- factor(comp.df$LC, 
                     levels=paste0("grnt_", 
                                   c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd")),
                     labels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
comp.df <- arrange(comp.df, cellID, LC)
phi.df <- cbind(comp.df, 
                matrix(ncol=length(v)*length(m), nrow=nrow(L), 
                       dimnames=list(NULL, 
                                     paste0(rep(v, each=2), "_", m, "_phi"))))
phi.df <- filter(phi.df, LC != "WP")
comp.df <- cbind(comp.df, 
                 matrix(ncol=length(v)*length(m), nrow=nrow(L), 
                        dimnames=list(NULL, paste0(rep(v, each=2), "_", m))))


########
## run models
########
n.chain <- 4
n.iter <- 500
for(j in 1:2) {
  for(i in seq_along(v)) {
    mod <- paste0(v[i], "_", m[j])
    phi <- paste0(mod, "_phi")
    # out <- stan(file=paste0("ms/", m[j], ".stan"),
    #             data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
    #             pars=c("nu", "Z_", "Z_new_"), include=F,
    #             iter=n.iter, chains=n.chain,
    #             sample_file=paste0("out_stan/", mod, ".csv"))
    out <- read_stan_csv(list.files("out_stan", paste0("^", mod), full.names=T))
    comp.df[[mod]] <- get_posterior_mean(out, pars="eta")[,n.chain+1]
    if(m[j]=="spatial") {
      phi.df[[phi]] <- get_posterior_mean(out, pars="phi")[,n.chain+1]
    }
  }
}

comp.long <- gather(comp.df, model, mn, 17:30) %>%
  mutate(Cov=str_split_fixed(comp.long$model, "_", 2)[,1],
         Space=str_split_fixed(comp.long$model, "_", 2)[,2])
comp.sp <- select(comp.long, -model) %>%
  spread(Space, mn)
  
ggplot(comp.sp, aes(x=lon, y=lat, fill=spatial)) + geom_tile() + 
  facet_grid(LC~Cov) + ggtitle("CAR random effects") + 
  scale_fill_gradient("Post. Mn", low="white", high="red", limits=c(0,1))
ggplot(comp.sp, aes(x=lon, y=lat, fill=nonspatial)) + geom_tile() + 
  facet_grid(LC~Cov) + ggtitle("No spatial random effects") + 
  scale_fill_gradient("Post. Mn", low="white", high="red", limits=c(0,1))
ggplot(comp.long, aes(x=lon, y=lat, fill=mn-nlcd)) + geom_tile() + 
  facet_grid(LC~model) +
  scale_fill_gradient2(low="blue", high="red", mid="white")
ggplot(comp.long, aes(x=lon, y=lat, fill=mn-grnt)) + geom_tile() + 
  facet_grid(LC~model) +
  scale_fill_gradient2(low="blue", high="red", mid="white")
ggplot(comp.sp, aes(x=lon, y=lat, fill=nonspatial-spatial)) + geom_tile() + 
  facet_grid(LC~Cov) +
  scale_fill_gradient2(low="blue", high="red", mid="white")

comp.long %>% group_by(Cov, Space, LC) %>%
  summarise(mn_grnt=mean(mn-grnt, na.rm=T),
            mn_nlcd=mean(mn-nlcd))
comp.sp %>% group_by(Cov, LC) %>%
  summarise(car_grnt=mean(spatial-grnt, na.rm=T),
            non_grnt=mean(nonspatial-grnt, na.rm=T),
            car_nlcd=mean(spatial-nlcd, na.rm=T),
            non_nlcd=mean(nonspatial-nlcd, na.rm=T)) %>%
  mutate(car_non_grnt=(abs(car_grnt)-abs(non_grnt))/abs(non_grnt),
         car_non_nlcd=(abs(car_nlcd)-abs(non_nlcd))/abs(non_grnt)) %>%
  mutate_if(is.numeric, round, 3) %>%
  print.AsIs

