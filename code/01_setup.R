# This script reads in a csv with a row for each cell in the landscape, 
# reshaping and generating the data objects required to run the stan models. 
# This requires the following input:
#   - data_f: csv file with col(landcover, covariates) and row(cells)
#   - var_f: csv file with variable categories for model comparison
# Two sets of data are created: one for variable selection and one for the full
# run. The variable selection dataset is limited to cells with data for both Y
# and Z, and these cells are then split into a training (80%) and testing subset
# (20%), with the testing subset in eastern New Hampshire nearest to the
# extrapolation region. In the dataset for the full run, the 'Train' column
# indicates whether the region is used for fitting or predicting (i.e.,
# extrapolating). An Rdump object is created in the folder data_stan/ for each
# covariate set, where datasets for variable selection are appended "_varSel".



########
## set up
########
pTrain <- 0.8  # proportion of cells in overlap region to fit model with 
K <- 6  # number of land cover categories
res <- c("200ha", "2000ha")[2]
data_f <- paste0("data/landscape_", res, ".csv")  # landscape csv
var_f <- "data/X_cat.csv"  # variable set identification

# load workspace
library(tidyverse); library(rstan)
k <- K-1
L <- read.csv(data_f) %>% arrange(desc(Train), CellID)
X_i <- read.csv(var_f)



########
## partition landscape for variable selection
########
L_varSel <- L %>% filter(Train==1) %>%  # df of cells with Y
  mutate(Train=CellID < quantile(CellID, probs=pTrain)) %>%
  arrange(desc(Train), CellID)
write_csv(L_varSel, paste0("data/L_varSel_", res, ".csv"))
write_csv(L, paste0("data/L_", res, ".csv"))



########
## identify covariate sets
########
# V: covariate matrix (scaled)
# V_varSel: list(covariateSets=list(delta[1:(k-1)], rho)) with V indexes
# X = V, except V includes pWP as an extra possible covariate. For simplicity,
# we define only V in the data for the model, then index V as needed to isolate
# the columns for X
V <- scale(as.matrix(L[,15:23]))  # full landscape
V_varSel <- scale(as.matrix(L_varSel[,15:23]))  # overlapping cells only
V.ls <- list(
  Cens=rep(list(which(X_i$cat %in% "cens")), k),
  Clim=rep(list(which(X_i$cat %in% "clim")), k),
  Topo=rep(list(which(X_i$cat %in% "topo")), k),
  CensClim=rep(list(which(X_i$cat %in% c("cens", "clim"))), k),
  CensTopo=rep(list(which(X_i$cat %in% c("cens", "topo"))), k),
  ClimTopo=rep(list(which(X_i$cat %in% c("clim", "topo"))), k),
  CensClimTopo=rep(list(which(X_i$cat %in% c("cens", "clim", "topo"))), k))
V.ls <- append(setNames(V.ls, paste0(names(V.ls), "WP")), V.ls)
# add pWP as first covariate in first rep of datasets
for(i in 1:(length(V.ls)/2)) V.ls[[i]][[5]] <- c(V.ls[[i]][[5]], 9)
  


########
## construct adjacency matrix
########
# following the Stan exact sparse CAR case study
# <http://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html>
# W: adjacency matrix with W[i,i]=0 and W[i,j]=ifelse(neighbors, 1, 0) 
# W_n: number of adjacent pairs (count only W[i,j] or W[j,i] -- not both)
n1 <- sum(L$Train)  # number of cells for fitting
n1_varSel <- sum(L_varSel$Train)  # number of cells for fitting
L$x <- as.numeric(as.factor(L$lon)) # grid col index
L$y <- as.numeric(as.factor(L$lat)) # grid row index
L_varSel$x <- as.numeric(as.factor(L_varSel$lon)) # grid col index
L_varSel$y <- as.numeric(as.factor(L_varSel$lat)) # grid row index
dist_mx <- as.matrix(dist(data.frame(L$x, L$y)))
dist_mx_varSel <- as.matrix(dist(data.frame(L_varSel$x, L_varSel$y)))
W <- dist_mx <= sqrt(2) # max 8 neighbors per cell
W_varSel <- dist_mx_varSel <= sqrt(2) # max 8 neighbors per cell
diag(W) <- 0 
diag(W_varSel) <- 0
W_n <- sum(W)/2 
W_n_varSel <- sum(W_varSel)/2 


########
## build stan datasets
########
b.ls <- list(  # base data list: all cells
  n1=n1, 
  n2=n1+1, 
  n3=nrow(L), 
  K=K,
  Y=as.matrix(select(L[1:n1,], grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Z=as.matrix(select(L, nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Evg)),
  W=W,
  W_n=W_n)
b_varSel <- list(  # base data list: variable selection
  n1=n1_varSel, 
  n2=n1_varSel+1, 
  n3=nrow(L_varSel), 
  K=K,
  Y=as.matrix(select(L_varSel[1:n1_varSel,], 
                     grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  new_Y=as.matrix(select(L_varSel[(n1_varSel+1):nrow(L_varSel),], 
                         grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg, grnt_Mxd)),
  Z=as.matrix(select(L_varSel, nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Evg)),
  W=W_varSel,
  W_n=W_n_varSel)

# function for adding covariate sets 
make_d <- function(V.j, # covariate set indexes
                   V, # all covariates
                   l # base data list
                   ) {
  nXV <- map_int(V.j, length) # number of covariates for delta, rho in cov set
  di <- rep(NA, 8)  
  di[c(1,3,5,7)] <- seq(1, nXV[1]*3+1, by=nXV[1])
  di[c(2,4,6,8)] <- seq(nXV[1], nXV[1]*4, by=nXV[1])
  l$nX=nXV[1] # number of betas (same for all delta's in cov set)
  l$nV=nXV[5] # number of thetas
  l$di=di
  l$V=V[1:l$n1,V.j[[5]]] # covariates in cov set
  l$new_V=V[l$n2:l$n3,V.j[[5]]]
  return(l)
}

Full.ls <- map(V.ls, make_d, V, b.ls)
varSel.ls <- map(V.ls, make_d, V_varSel, b_varSel)


########
## store stan datasets
########
for(j in seq_along(V.ls)) {
  rstan::stan_rdump(ls(Full.ls[[j]]),
                    file=paste0("data_stan/", names(Full.ls)[j], 
                                "_", res, ".Rdump"),
                    envir=list2env(Full.ls[[j]]))
  rstan::stan_rdump(ls(varSel.ls[[j]]),
                    file=paste0("data_stan/", names(varSel.ls)[j], 
                                "_", res, "_varSel.Rdump"),
                    envir=list2env(varSel.ls[[j]]))
}


