
library(tidyverse)
load("data/lc_base.Rdata")


# For making different sized datasets
data_full <- data_df
nFit <- sum(!data_full$Set)

# For different sizes
sizes <- c("50k"=50000, "75k"=75000, "100k"=100000, "125k"=125000)
id <- purrr::map(sizes, ~sample(1:nFit, ., replace=FALSE))

for(s in 1:length(sizes)) {
  
  data_df <- data_full[id[[s]], ]
  # set base data that is identical across all runs
  
  d.b <- list(n1=sizes[s], L=6, #n2=nFit+1, n3=nrow(data_df), L=6,
              Y1=as.matrix(data_df[, c(8:10,12:13)]), # OpI,OpU,Dec,Evg,WP
              Y2=as.matrix(data_df[, c(14:16,18)]))  # OpI,OpU,Dec,Evg
  
  # arrange covariate sets
  X.all <- scale(with(data_df, cbind(el_mean, rugg_mean,
                                     bio1_mean, bio7_mean, bio12_mean,
                                     pop, hms, rdLen_,
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
  
  make_d <- function(X.i, mod.nm, X.all, d.b) {
    nX <- map_int(X.i, length)
    d_i <- rep(NA, 8)
    d_i[c(1,3,5,7)] <- seq(1, nX[1]*3+1, by=nX[1])
    d_i[c(2,4,6,8)] <- seq(nX[1], nX[1]*4, by=nX[1])
    d.b$nB_d=nX[1]; d.b$nB_p=nX[5]; d.b$di=d_i
    d.b$X=X.all[,X.i[[5]]]
    d.b$mod=mod.nm
    return(d.b)
  }
  
  d.ls <- map2(X.ls, names(X.ls), make_d, X.all, d.b)
  
  for(i in c(1,3,5,6,7)) {
    rstan::stan_rdump(ls(d.ls[[i]][1:8]),
                      file=paste0("~/Desktop/lc/data/comp/", names(d.ls)[i], 
                                  "_", names(sizes)[s], ".Rdump"),
                      envir=list2env(d.ls[[i]]))
  }
}


