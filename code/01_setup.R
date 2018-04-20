# This script reads in a csv with a row for each cell in the landscape, 
# reshaping and generating the data objects required to run the stan models. 
# This requires the following input:
#   - data_f: csv file with col(landcover, covariates) and row(cells)
#   - var_f: csv file with variable categories for model comparison
# A stratified random sample is generated for model comparison where the sample
# is representative of the ntiles within tol %. 
# An Rdump object is created in the folder data_stan/ for each covariate set,
# where datasets for variable selection are appended "_VS".

########
## set up
########
pSamp <- 0.85  # proportion of cells in overlap region to fit model with 
tol <- 0.05  # % error tolerance for ntile representation in sample
q <- 5  # number of bins for assessing representativeness of sample
D <- 6  # number of land cover categories
data_f <- "data/landscape_3km.csv"  # landscape csv
var_f <- "data/X_cat.csv"  # variable set identification
samp.f <- paste0("data/stratSamp_", pSamp*100, ".txt") # for VS sample indexes

# load workspace
library(tidyverse); library(rstan)
d <- D-1
L <- read.csv(data_f) %>% arrange(desc(Fit), cellID)
X_i <- read.csv(var_f)
L.YZ <- L %>% filter(Fit==1) %>% arrange(lon, lat) # df of cells with Y & Z
X <- L.YZ %>% select(pWP, temp, seas, precip, el, rugg, roads, pop, homes)


########
## identify indexes for representative sample for variable selection
########
X.q <- X %>% mutate_all(ntile, n=q)  # identify ntiles for each covariate
n <- floor(nrow(X)*pSamp)  # n.cells in sample
samp.id <- sample(1:nrow(X), n, replace=FALSE)  # sample row indexes
x.sum <- sapply(X.q[samp.id,], table)  # quartile representation in sample

i <- 0  # brute force: resample until acceptibly representative ntile counts
while( sum(abs(x.sum - n/q) > tol*n/q) != 0 ) {
  samp.id <- sample(1:nrow(X), n, replace=FALSE)
  x.sum <- sapply(X.q[samp.id,], table)
  i <- i+1
  if((i %% 1000) == 0) cat("Still trying...", i, "attempts\n")
}
x.sum
cat(samp.id, file=samp.f)  # store row IDs in sample


########
## partition landscape for variable selection
########
samp.id <- scan(samp.f)
L.YZ$Fit <- 1:nrow(L.YZ) %in% samp.id
L.YZ <- L.YZ %>% arrange(desc(Fit), cellID)


########
## identify covariate sets
########
# X.all: covariate matrix (scaled)
# X.ls: list(covariateSets=list(rho[1:(d-1)], p)) with X.all indexes
X.all <- scale(as.matrix(L[,15:23]))  # full landscape
X.YZ <- scale(as.matrix(L.YZ[,15:23]))  # overlapping cells only
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
L.YZ$x <- as.numeric(as.factor(L.YZ$lon)) # grid col index
L.YZ$y <- as.numeric(as.factor(L.YZ$lat)) # grid row index
dist_mx <- as.matrix(dist(data.frame(L$x, L$y)))
dist_mx.YZ <- as.matrix(dist(data.frame(L.YZ$x, L.YZ$y)))
W <- dist_mx <= sqrt(2) # max 8 neighbors per cell
W.YZ <- dist_mx.YZ <= sqrt(2) # max 8 neighbors per cell
diag(W) <- diag(W.YZ) <- 0
W_n <- sum(W)/2 
W_n.YZ <- sum(W.YZ)/2 


########
## build stan datasets
########
n1 <- sum(L$Fit)  # number of cells for fitting
n1.YZ <- sum(L.YZ$Fit)  # number of cells for fitting
b.ls <- list(  # base data list: all cells
  n1=n1, 
  n2=n1+1, 
  n3=nrow(L), 
  D=D,
  Y=as.matrix(select(L[1:n1,], grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Z=as.matrix(select(L, nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Evg)),
  W=W,
  W_n=W_n)
b.YZ <- list(  # base data list: variable selection
  n1=n1.YZ, 
  n2=n1.YZ+1, 
  n3=nrow(L.YZ), 
  D=D,
  Y=as.matrix(select(L.YZ[1:n1.YZ,], 
                     grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Y_new=as.matrix(select(L.YZ[(n1.YZ+1):nrow(L.YZ),], 
                         grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Z=as.matrix(select(L.YZ, nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Evg)),
  W=W.YZ,
  W_n=W_n.YZ)

# function for adding covariate sets 
make_d <- function(X.j, # covariate set indexes
                   X.all, # all covariates
                   l # base data list
                   ) {
  nX <- map_int(X.j, length) # number of covariates for rho, p in cov set
  ri <- rep(NA, 8)  
  ri[c(1,3,5,7)] <- seq(1, nX[1]*3+1, by=nX[1])
  ri[c(2,4,6,8)] <- seq(nX[1], nX[1]*4, by=nX[1])
  l$n_b=nX[1] # number of betas (same for all rho's in cov set)
  l$n_t=nX[5] # number of thetas
  l$ri=ri
  l$X=X.all[1:l$n1,X.j[[5]]] # covariates in cov set
  l$X_new=X.all[l$n2:l$n3,X.j[[5]]]
  return(l)
}

all.ls <- map(X.ls, make_d, X.all, b.ls)
YZ.ls <- map(X.ls, make_d, X.YZ, b.YZ)


########
## store stan datasets
########
for(j in seq_along(all.ls)) {
  rstan::stan_rdump(ls(all.ls[[j]]),
                    file=paste0("data_stan/", names(all.ls)[j], ".Rdump"),
                    envir=list2env(all.ls[[j]]))
  rstan::stan_rdump(ls(YZ.ls[[j]]),
                    file=paste0("data_stan/", names(YZ.ls)[j], "_VS.Rdump"),
                    envir=list2env(YZ.ls[[j]]))
}


