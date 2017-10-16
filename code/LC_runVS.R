# This script generates a list of covariate sets, then runs the latent model
# without spatial random effects


##########
## set up
##########

if(!require("sevcheck")) devtools::install_github("Sz-Tim/sevcheck")
p_load("tidyverse", "magrittr", "forcats", "rstan", "ggmcmc", "stringr")
rstan_options(auto_write=TRUE); options(mc.cores=parallel::detectCores())
source("code/LC_fun.R"); source("code/stan_utilities.R"); theme_set(theme_bw())

grdSz <- "01_1a"
blockSize <- 10  # block = (blockSize x blockSize) grid cells
maxGridRow <- 40  # number of blocks per row; 1-40
maxGridCol <- 40  # number of blocks per column; 1-54

# cell-block reference tibble
cb.i <- read_csv(paste0("data/roads_", grdSz, ".csv")) %>% 
  mutate(CellRow=1:n_distinct(top) %>% rep(n_distinct(left)),
         CellCol=1:n_distinct(left) %>% rep(each=n_distinct(top))) %>%
  filter((CellRow <= max((CellRow %/% blockSize) * blockSize)) &
           (CellCol <= max((CellCol %/% blockSize) * blockSize))) %>%
  mutate(BlockRow=((CellRow-1)%/%blockSize)+1, 
         BlockCol=((CellCol-1)%/%blockSize)+1,
         BlockID=paste(str_pad(BlockCol, 7, "left", "0"), 
                       str_pad(BlockRow, 7, "left", "0"), sep=".") %>% 
           as.numeric %>% factor %>% as.numeric) %>%
  select(c(CellID, CellRow, CellCol, BlockID, BlockRow, BlockCol, left, top))
Block.inc <- cb.i$BlockID[cb.i$BlockCol <= maxGridCol &
                            cb.i$BlockRow <= maxGridRow] %>% unique

# covariates summarized to blocks
pop00 <- read_csv(paste0("data/pop00_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(popTot=log(sum(sum)+0.001)) %>%
  filter(BlockID %in% Block.inc)
hous00 <- read_csv(paste0("data/housing00_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(secHome=log(sum(sum)+0.001)) %>%
  filter(BlockID %in% Block.inc)
rdLen <- read_csv(paste0("data/roads_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(rdLen=log(sum(roadLen)+0.001))  %>%
  filter(BlockID %in% Block.inc)
clim <- read_csv(paste0("data/clim_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% 
  summarise(b1=mean(bio1_mean), b7=mean(bio7_mean), b12=mean(bio12_mean)) %>%
  filter(BlockID %in% Block.inc)
topo <- read_csv(paste0("data/topo_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% 
  summarise(el=mean(el_mean), rugg=mean(rugg_mean)) %>%
  filter(BlockID %in% Block.inc)
pWP <- read_csv(paste0("data/pWP_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>%
  add_blocks(cb.i=cb.i) %>% summarise(mnWP=mean(mean)/100) %>%
  filter(BlockID %in% Block.inc)

# land cover summarized to blocks
grnt <- read_csv(paste0("data/out_", grdSz, "_grnt.csv")) %>% 
  mutate(CellID=1:nrow(.)) %>% add_blocks(cb.i=cb.i) %>% 
  summarise(Dev=sum(V1)/n(), Oth=sum(V2)/n(), Hwd=sum(V3)/n(), 
            WP=sum(V4)/n(), Evg=sum(V5)/n(), Mxd=sum(V6)/n()) %>%
  filter(BlockID %in% Block.inc) %>%
  select(-BlockID) %>% as.matrix
nlcd <- read_csv(paste0("data/out_",grdSz,"_nlcd.csv"))  %>% 
  mutate(CellID=1:nrow(.)) %>% add_blocks(cb.i=cb.i) %>% 
  summarise(Dev=sum(V1)/n(), Oth=sum(V2)/n(), Hwd=sum(V3)/n(), 
            Evg=sum(V4)/n(), Mxd=sum(V5)/n()) %>%
  filter(BlockID %in% Block.inc) %>%
  select(-BlockID) %>% as.matrix

set.seed(22222)
nFit <- 960
nNew <- 640
n <- sampleCells(nFit, nNew, nrow(grnt), partition=TRUE)

# Y1 & Y2
Y1.fit <- grnt[n$fit,]
Y1.new <- grnt[n$new,]
Y2 <- nlcd[n$all,]



##########
## generate covariate sets
##########

X.all <- scale(cbind(topo$el[n$all], topo$rugg[n$all],
               clim$b1[n$all], clim$b12[n$all], clim$b7[n$all],
               pop00$popTot[n$all], hous00$secHome[n$all], rdLen$rdLen[n$all],
               pWP$mnWP[n$all]))
colnames(X.all) <- c("el", "rug", 
                     "temp", "precip", "seas",
                     "pop", "homes", "rdLen",
                     "pWP")
X.ls <- list(
  m1=list(c(1,2,5,6), c(2,4,5,7), c(1,3,4,5), c(5,6,7,8), c(1,3,5,9)),
  m2=list(c(1,2,5,6), c(2,4,5,7), c(1,3,4), c(5,6,7,8), c(1,3,5,9)),
  m3=list(c(1,2,5,6), c(2,4,5,7), c(1,3), c(5,6,7,8), c(1,3,5,9))
)
nX.ls <- map(X.ls, ~map_int(., length))

make_d <- function(X.i, X.all, nFit, n, L, Y1.fit, Y2) {
  # X.i = X.ls[[i]] = covariate column indices for each LC
  # X.all = covariates
  nX <- map_int(X.i, length)
  d <- list(n1=nFit, n2=nFit+1, n3=n$tot, L=6, Y1=Y1.fit, Y2=Y2,
            nB_d=nX[1:4], nB_p=nX[5],
            X_d1=X.all[,X.i[[1]]], X_d2=X.all[,X.i[[2]]], X_d3=X.all[,X.i[[3]]],
            X_d4=X.all[,X.i[[4]]], X_p=X.all[,X.i[[5]]])
  return(d)
}

d.ls <- d.ls <- map(X.ls, make_d, X.all, nFit, n, 6, Y1.fit[,-6], Y2[,-5])





##########
## run stan in batch
##########

run_stan <- function(d, mod, nChain=1, iter=2, warmup=1) {
  out <- stan(file=paste0("code/", mod), 
              data=d, 
              iter=iter, warmup=warmup, chains=nChain, seed=4337, init=0,
              include=FALSE, pars=c("Y2_", "Y2new_", "nu"))
  return(out)
}

out.ls <- map(d.ls, run_stan, "sc_theta/direct.stan",
              nChain=2, iter=2000, warmup=1000)

beta.order <- map(X.ls, `[`, 1:4) %>% map(unlist)
beta.index <- map(nX.ls, ~rep(1:4, times=.[1:4]))
beta.labels <- map(beta.order, ~paste0("theta_d_z[", 1:length(.), "]"))
beta.vars <- map(beta.order, ~colnames(X.all)[.])

add_beta_id <- function(x, y) {
  x$var <- y[match(x$Parameter, y[,1]), 3]
  x$bias <- y[match(x$Parameter, y[,1]), 2]
  return(x)
}
gg.beta.ls <- map2(map(out.ls, ~ggs(., "theta_d_z")), 
                   pmap(list(beta.labels, beta.index, beta.vars), cbind), 
                   add_beta_id)
map(gg.beta.ls, ~ggplot(., aes(x=value)) + geom_density() + xlim(-5, 5) + 
      facet_grid(bias~var) + geom_vline(xintercept=0, linetype=3))


