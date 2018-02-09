# This script creates the datasets to be used by the land cover model. It uses
# the output from GIS (a row for each cell with left, top, CellID, pWP, climate,
# GRANIT, NLCD, topo, census, & Set=(fit vs predict)) and creates datasets in 
# the correct format for the stan model.

library(tidyverse); library(magrittr)
lc_5km <- read_csv("data/landcover/landcover_5km.csv")


# clean GIS output
for(i in 1:nrow(lc_5km)) {
  nlcd_sum <- sum(lc_5km[i,4:8])
  if(nlcd_sum==0 | is.na(nlcd_sum)) {
    lc_5km[i,4:8] <- NA
  } else {
    lc_5km[i,4:8] <- lc_5km[i,4:8]/nlcd_sum
  }
  nhlc_sum <- sum(lc_5km[i,9:14])
  if(nhlc_sum==0 | is.na(nhlc_sum)) {
    lc_5km[i,9:14] <- NA
  } else {
    lc_5km[i,9:14] <- lc_5km[i,9:14]/nhlc_sum
  }
}
lc_5km %<>% filter(!is.na(nlcd1_mean) & !is.na(rdLen) & 
                     !is.na(pop_sum) & !is.na(bio1_mean))
lc_5km$Set <- !is.na(lc_5km$nhlc1_mean)  # T: fit, F: predict
lc_5km %<>% arrange(!Set, left, top)


# make base stan_data
nFit <- sum(lc_5km$Set)
d.b <- list(n1=nFit, n2=nFit+1, n3=nrow(lc_5km), D=6, 
            Y1=as.matrix(lc_5km[1:nFit, c(9:11,13:14)]), # OpI,OpU,Dec,Evg,WP
            Y2=as.matrix(lc_5km[, c(4:6,8)]))  # OpI,OpU,Dec,Evg

# arrange covariate sets
X.all <- scale(with(lc_5km, cbind(el_mean, rugg_mean,
                                   bio1_mean, bio7_mean, bio12_mean,
                                   pop_sum, hms_sum, rdLen,
                                   pWP_mean)))
colnames(X.all) <- c("el", "rugg",
                     "temp", "seas", "precip",
                     "pop", "homes", "roads",
                     "pWP")
X.id <- tibble(var=colnames(X.all),
               cat=rep(c("topo", "clim", "census", "pWP"), c(2,3,3,1)))
X.ls <- list(
  Topo=rep(list(which(X.id$cat %in% "topo")), 5),
  Clim=rep(list(which(X.id$cat %in% "clim")), 5),
  Cens=rep(list(which(X.id$cat %in% "census")), 5),
  TopoClim=rep(list(which(X.id$cat %in% c("topo", "clim"))), 5),
  TopoCens=rep(list(which(X.id$cat %in% c("topo", "census"))), 5),
  ClimCens=rep(list(which(X.id$cat %in% c("clim", "census"))), 5),
  TopoClimCens=rep(list(which(X.id$cat %in% c("topo", "clim", "census"))), 5)
)
for(i in 1:length(X.ls)) {
  X.ls[[i]][[5]] <- c(9, X.ls[[i]][[5]])  # add pWP
}
nX.ls <- map(X.ls, ~map_int(., length))

# store in list of lists and dump
make_d <- function(X.i, mod.nm, X.all, d.b) {
  nX <- map_int(X.i, length)
  d_i <- rep(NA, 8)
  d_i[c(1,3,5,7)] <- seq(1, nX[1]*3+1, by=nX[1])
  d_i[c(2,4,6,8)] <- seq(nX[1], nX[1]*4, by=nX[1])
  d.b$nB=nX[1]
  d.b$nTh=nX[5]
  d.b$ri=d_i
  d.b$Z=X.all[,X.i[[5]]]
  d.b$mod=mod.nm
  return(d.b)
}
d.ls <- map2(X.ls, names(X.ls), make_d, X.all, d.b)
for(i in 1:7) {
  rstan::stan_rdump(ls(d.ls[[i]][1:10]),
                    file=paste0("data/landcover/", names(d.ls)[i], ".Rdump"),
                    envir=list2env(d.ls[[i]]))
}

# run models
library(rstan)
options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
out <- stan(file="code/stan/lc_full.stan", data=d.ls[[1]][1:10], iter=1000)






