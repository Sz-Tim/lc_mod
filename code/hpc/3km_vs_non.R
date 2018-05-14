library(stringr); library(doSNOW)

n.chain <- 3
n.iter <- 10000
p.c <- makeCluster(7); registerDoSNOW(p.c)

v <- str_remove(list.files("data_stan", "_VS"), ".Rdump")
foreach(i=seq_along(v)) %dopar% {
  library(rstan)
  options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
  mod <- paste0(v[i], "_nonspatial_vs")
  out <- stan(file="R/nonspatial_vs.stan",
              data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
              pars=c("nu", "Z_", "Z_new_"), include=F,
              iter=n.iter, chains=n.chain, refresh=50, thin=10,
              sample_file=paste0("out_stan/", mod, ".csv"))
}
stopCluster(p.c)
