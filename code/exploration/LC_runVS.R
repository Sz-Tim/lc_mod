# This script generates a list of covariate sets, then runs the latent model
# without spatial random effects


##########
## set up
##########

if(!require("sevcheck")) devtools::install_github("Sz-Tim/sevcheck")
p_load("tidyverse", "magrittr", "forcats", "rstan", "ggmcmc", 
       "stringr", "parallel", "loo")
rstan_options(auto_write=TRUE); options(mc.cores=parallel::detectCores())
source("code/LC_fun.R"); source("code/stan_utilities.R"); theme_set(theme_bw())

grdSz <- "01_1a"
blockSize <- 2  # block = (blockSize x blockSize) grid cells
maxGridRow <- 204  # number of blocks per row; 1-40
maxGridCol <- 271  # number of blocks per column; 1-54


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
nFit <- 47000
nNew <- 8284
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
X.id <- tibble(var=colnames(X.all),
              cat=rep(c("topo", "clim", "census", "pWP"), c(2,3,3,1)))
X.ls <- list(
  tp=rep(list(which(X.id$cat %in% "topo")), 5), 
  cl=rep(list(which(X.id$cat %in% "clim")), 5), 
  cn=rep(list(which(X.id$cat %in% "census")), 5), 
  tp_cl=rep(list(which(X.id$cat %in% c("topo", "clim"))), 5), 
  tp_cn=rep(list(which(X.id$cat %in% c("topo", "census"))), 5),
  cl_cn=rep(list(which(X.id$cat %in% c("clim", "census"))), 5),
  tp_cl_cn=rep(list(which(X.id$cat %in% c("topo", "clim", "census"))), 5)
)
for(i in 1:length(X.ls)) {
  X.ls[[i]][[5]] <- c(9, X.ls[[i]][[5]])
}
nX.ls <- map(X.ls, ~map_int(., length))

make_d <- function(X.i, mod.nm, X.all, nFit, n, L, Y1.fit, Y2) {
  # X.i = X.ls[[i]] = covariate column indices for each LC
  # X.all = covariates
  nX <- map_int(X.i, length)
  d_i <- rep(NA, 8)
  d_i[c(1,3,5,7)] <- seq(1, nX[1]*3+1, by=nX[1])
  d_i[c(2,4,6,8)] <- seq(nX[1], nX[1]*4, by=nX[1])
  d <- list(n1=nFit, n2=nFit+1, n3=n$tot, L=6, Y1=Y1.fit, Y2=Y2,
            nB_d=nX[1], nB_p=nX[5], di=d_i,
            X=X.all[,X.i[[5]]], 
            mod=mod.nm)
  return(d)
}

d.ls <- map2(X.ls, names(X.ls), make_d, X.all, nFit, n, 6, Y1.fit[,-6], Y2[,-5])
for(i in 1:7) {
  stan_rdump(ls(d.ls[[i]]),
             file=paste0("../null_dir/lc/data/", names(d.ls)[i],".Rdump"),
             envir=list2env(d.ls[[i]]))
}





##########
## run stan in batch
##########


run_stan <- function(d, nChain=1, iter=2, warmup=1) {
  out <- stan(file="code/sc_theta/latent_bp.stan", 
              data=d, 
              iter=iter, warmup=warmup, chains=nChain, seed=4337, init=0,
              include=FALSE, pars=c("Y2_", "Y2new_", "nu", 
                                    "bias", "bias_new", "pWP", "pWP_new"))
  return(out)
}

mod.ls <- paste0("sc_theta/latent_", c("bp.stan", "vs.stan"))

out.ls <- map(d.ls[1:2], run_stan, nChain=3, iter=2000, warmup=1000)

out.ls <- mclapply(d.ls[1:3], FUN=run_stan, nChain=2, iter=2000, warmup=1000, 
                   mc.preschedule=FALSE)

map(out.ls, ~waic(extract_log_lik(.)))
map_df(out.ls, ~rowSums(get_elapsed_time(.))) %>% t %>% cbind(rowMeans(.))
map_dbl(out.ls, ~(sum(get_elapsed_time(.))/summary(., pars="lp__")$summary[9]))


beta.order <- map(X.ls, `[`, 1:5) %>% map(unlist)
beta.index <- map(nX.ls, ~c(rep(1:5, times=.[1:5])))
beta.labels <- map(beta.order, ~c(paste0("theta_d_z[", 
                                         1:(which(.==9)-1), 
                                         "]"),
                                  paste0("theta_p[", 
                                         1:(length(.)-which(.==9)+1),
                                         "]")))
beta.vars <- map(beta.order, ~colnames(X.all)[.])

add_beta_id <- function(x, y) {
  x$var <- y[match(x$Parameter, y[,1]), 3]
  x$bias <- y[match(x$Parameter, y[,1]), 2]
  return(x)
}
gg.ls <- map(out.ls, ~(ggs(., "theta") %>% 
                         filter(Parameter %in% beta.labels$tp_cl_cn)))
gg.beta.ls <- map2(gg.ls, 
                   pmap(list(beta.labels, beta.index, beta.vars), cbind), 
                   add_beta_id)
gg.beta <- bind_rows(gg.beta.ls, .id="model")
ggplot(gg.beta, aes(x=value, colour=model)) + geom_density() + 
      facet_grid(bias~var) + geom_vline(xintercept=0, linetype=3)
ggplot(gg.beta, aes(x=bias, y=value, colour=model)) + geom_violin() +
      facet_wrap(~var) + geom_hline(yintercept=0, linetype=3)

