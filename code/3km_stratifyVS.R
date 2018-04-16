# This script performs stratified random sampling to create the variable 
# selection grid. The full grid takes too long to run, so variable selection is
# done on a subset of the cells, with each variable represented proportionally
# in the subset as it is found on the landscape. As the variable selection is 
# only possible for cells with both GRANIT & NLCD, it is restricted to cells 
# with data_df$Set == FALSE.

library(tidyverse)
D <- 6
d <- D-1
# identify csv's
data_f <- "data/landscape_3km.csv"
var_f <- "data/X_cat.csv"
L <- read.csv(data_f) %>% arrange(desc(Fit), cellID)
X_i <- read.csv(var_f)

# extract covariates
NH <- L %>% filter(Fit==1) %>%
  arrange(lon, lat) %>% 
  mutate(id=row_number(), 
         x=as.numeric(as.factor(str_pad(lon, 3, "left", "0"))),
         y=as.numeric(as.factor(str_pad(lat, 3, "left", "0"))))
X <- NH %>% 
  select(pWP, temp, seas, precip,
         el, rugg, roads, pop, homes)

# calculate quartiles
X.q <- X %>% mutate_all(ntile, n=5)

n <- floor(nrow(X)*0.15)

samp.id <- sample(1:nrow(X), n, replace=FALSE)
x.sum <- sapply(X.q[samp.id,], table)
x.sum

n.q <- colSums(x.sum)[1]/n_distinct(X.q$el)

i <- 0
while( sum(abs(x.sum - n.q) > 0.125*n.q) != 0 ) {
  samp.id <- sample(1:nrow(X), n, replace=FALSE)
  x.sum <- sapply(X.q[samp.id,], table)
  i <- i+1
  if((i %% 100) == 0) cat("Still trying...", i, "attempts\n")
}
x.sum

X.q.samp <- X.q[samp.id,]
X.q.inv <- X.q[-samp.id,]
names(X.q.samp) <- names(X.q.inv) <- paste0(names(X.q.samp), "_q")
X_df <- cbind(NH[samp.id,], X.q.samp)
X_inv <- cbind(NH[-samp.id,], X.q.inv)

write_csv(data.frame(id=samp.id), "data/stratified_sample_15pct_3km_rowID.csv")

theme_set(theme_bw())
ggplot(X_df, aes(x=x, y=y, fill=factor(pWP_q))) + 
  geom_tile() + scale_fill_brewer(palette=8)
ggplot(X_inv, aes(x=x, y=y, fill=factor(pWP_q))) + 
  geom_tile() + scale_fill_brewer(palette=8)


########
## set up L
########
NH$Fit <- ifelse(1:nrow(NH) %in% samp.id, 0, 1)
NH <- NH %>% arrange(desc(Fit), cellID)


########
## identify covariate sets
########
# X.all: covariate matrix (scaled)
# X.ls: list(covariateSets=list(rho[1:(d-1)], p)) with X.all indexes
X.all <- scale(as.matrix(NH[,15:23]))
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
NH$x <- as.numeric(as.factor(NH$lon)) # grid col index
NH$y <- as.numeric(as.factor(NH$lat)) # grid row index
dist_mx <- as.matrix(dist(data.frame(NH$x, NH$y)))
W <- dist_mx <= sqrt(2) # max 8 neighbors per cell
diag(W) <- 0
W_n <- sum(W)/2 


########
## build stan datasets
########
n1 <- sum(NH$Fit)  # number of cells for fitting
# base data list
base <- list(
  n1=n1, 
  n2=n1+1, 
  n3=nrow(NH), 
  D=D,
  Y=as.matrix(select(NH[1:n1,], 
                     grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Y_new=as.matrix(select(NH[(n1+1):nrow(NH),], 
                         grnt_Opn, grnt_Oth, grnt_Dec, grnt_WP, grnt_Evg)),
  Z=as.matrix(select(NH, nlcd_Opn, nlcd_Oth, nlcd_Dec, nlcd_Evg)),
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
                    file=paste0("data_stan/", names(data.ls)[j], "_VS.Rdump"),
                    envir=list2env(data.ls[[j]]))
}




