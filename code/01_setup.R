# This script reads in a csv with a row for each cell in the landscape, 
# reshaping and generating the data objects required to run the stan models. 
# This requires the following input:
#   - data_f: csv file with col(landcover, covariates) and row(cells)
#   - var_f: csv file with variable categories for model comparison
# An Rdump object is created in the folder data_stan/ for each covariate set.

########
## set up workspace
########
D <- 6
d <- D-1
# identify csv's
data_f <- "data/landscape_3km.csv"
var_f <- "data/X_cat.csv"
# load workspace
library(tidyverse); library(rstan)
L <- read.csv(data_f) %>% arrange(desc(Fit), cellID)
X_i <- read.csv(var_f)


########
## identify covariate sets
########
# X.all: covariate matrix (scaled)
# X.ls: list(covariateSets=list(rho[1:(d-1)], p)) with X.all indexes
X.all <- scale(as.matrix(L[,15:23]))
X.ls <- list(
  Cens=rep(list(which(X_i$cat %in% "cens")), d),
  Clim=rep(list(which(X_i$cat %in% "clim")), d),
  Topo=rep(list(which(X_i$cat %in% "topo")), d),
  CensClim=rep(list(which(X_i$cat %in% c("cens", "clim"))), d),
  CensTopo=rep(list(which(X_i$cat %in% c("cens", "topo"))), d),
  ClimTopo=rep(list(which(X_i$cat %in% c("clim", "topo"))), d),
  CensClimTopo=rep(list(which(X_i$cat %in% c("cens", "clim", "topo"))), d))
# add pWP as first covariate
for(i in seq_along(X.ls)) X.ls[[i]][[5]] <- c(9, X.ls[[i]][[5]])  


########
## construct adjacency matrix
########
# following the Stan exact sparse CAR case study
# <http://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html>
# W: adjacency matrix with W[i,i]=0 and W[i,j]=ifelse(neighbors, 1, 0) 
# W_n: number of adjacent pairs (count only W[i,j] or W[j,i] -- not both)
L$x <- as.numeric(as.factor(L$lon)) # grid col index
L$y <- as.numeric(as.factor(L$lat)) # grid row index
dist_mx <- as.matrix(dist(data.frame(L$x, L$y)))
W <- dist_mx <= sqrt(2) # max 8 neighbors per cell
diag(W) <- 0
W_n <- sum(W)/2 


########
## build stan datasets
########
n1 <- sum(L$Fit)  # number of cells for fitting
# base data list
base <- list(
  n1=n1, 
  n2=n1+1, 
  n3=nrow(L), 
  D=D,
  Y=as.matrix(select(L[1:n1,], grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Z=as.matrix(select(L, nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Evg)),
  W=W,
  W_n=W_n)
# function for assembling covariate sets
make_d <- function(X.j, X.all, base) {
  nX <- map_int(X.j, length) # number of covariates for rho, p in cov set
  ri <- rep(NA, 8)  
  ri[c(1,3,5,7)] <- seq(1, nX[1]*3+1, by=nX[1])
  ri[c(2,4,6,8)] <- seq(nX[1], nX[1]*4, by=nX[1])
  base$n_b=nX[1] # number of betas (same for all rho's in cov set)
  base$n_t=nX[5] # number of thetas
  base$ri=ri
  base$X=X.all[1:base$n1,X.j[[5]]] # covariates in cov set
  base$X_new=X.all[base$n2:base$n3,X.j[[5]]]
  return(base)
}

data.ls <- map(X.ls, make_d, X.all, base)


########
## store stan datasets
########
for(j in seq_along(data.ls)) {
  rstan::stan_rdump(ls(data.ls[[j]]),
                    file=paste0("data_stan/", names(data.ls)[j], ".Rdump"),
                    envir=list2env(data.ls[[j]]))
}


