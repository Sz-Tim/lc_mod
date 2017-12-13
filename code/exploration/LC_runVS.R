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
load("data/lc_base.Rdata")
strat.id <- read_csv("data/stratified_sample_05pct_20a_rowID.csv")

# For extracting a stratified random sample within NH (LC_stratifyVS.R)
NH_df <- data_df %>% filter(!Set) %>%
  arrange(left, top)
strat_df <- NH_df[strat.id$id,]

d.b <- list(n1=nrow(strat_df), L=6, #n2=nFit+1, n3=nrow(data_df), L=6,
            Y1=as.matrix(strat_df[, c(8:10,12:13)]), # OpI,OpU,Dec,Evg,WP
            Y2=as.matrix(strat_df[, c(14:16,18)]))  # OpI,OpU,Dec,Evg

# arrange covariate sets
X.all <- scale(with(strat_df, cbind(el_mean, rugg_mean,
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





##########
## run stan in batch
##########


for(i in 1:length(d.ls)) {
  invisible(stan(file="code/lc_VS.stan",
              data=d.ls[[i]][1:8], 
              model_name=d.ls[[i]]$mod,
              sample_file=paste0("out/1pct_", d.ls[[i]]$mod, ".csv"),
              iter=5000, warmup=2500, chains=6, 
              seed=4337, init=0))
  cat("Finished", i, "\n")
}

out.ls <- map(names(d.ls), ~list.files("out", pattern=paste0("1pct_", ., "_"), 
                                       full.names=TRUE)) %>%
  map(read_stan_csv)

out <- list.files("out", pattern=paste0("1pct_", names(d.ls)[1], "_"),
                  full.names=TRUE) %>% read_stan_csv %>% As.mcmc.list

map(out.ls, ~waic(extract_log_lik(.)))
map_df(out.ls, ~rowSums(get_elapsed_time(.))) %>% t %>% cbind(rowMeans(.))
map_dbl(out.ls, ~(sum(get_elapsed_time(.))/summary(., pars="lp__")$summary[9]))


beta.order <- map(X.ls, `[`, 1:5) %>% map(unlist)
beta.index <- map(nX.ls, ~c(rep(1:5, times=.[1:5])))
beta.labels <- map(beta.order, ~c(paste0("theta_d_z.", 
                                         1:(which(.==9)-1)),
                                  paste0("theta_p.", 
                                         1:(length(.)-which(.==9)+1))))
beta.vars <- map(beta.order, ~colnames(X.all)[.])

add_beta_id <- function(x, y) {
  x$var <- y[match(x$Parameter, y[,1]), 3]
  x$bias <- y[match(x$Parameter, y[,1]), 2]
  return(x)
}
gg.ls <- map(out.ls, ~(ggs(., "theta") %>% 
                         filter(Parameter %in% beta.labels$TopoClimCens)))
gg.beta.ls <- map2(gg.ls, 
                   pmap(list(beta.labels, beta.index, beta.vars), cbind), 
                   add_beta_id)
gg.beta <- bind_rows(gg.beta.ls, .id="model")
ggplot(gg.beta, aes(x=value, colour=model)) + geom_density() + 
      facet_grid(bias~var) + geom_vline(xintercept=0, linetype=3)
ggplot(gg.beta, aes(x=bias, y=value, colour=model)) + geom_violin() +
      facet_wrap(~var) + geom_hline(yintercept=0, linetype=3)

