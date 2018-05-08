library(stringr); library(doSNOW)

n.chain <- 3
n.iter <- 10000
p.c <- makeCluster(7); registerDoSNOW(p.c)

v <- str_remove(grep("_VS", list.files("data_stan"), inv=T, value=T), ".Rdump")
foreach(i=seq_along(v)) %dopar% {
  library(rstan)
  options(mc.cores=parallel::detectCores()); rstan_options(auto_write=TRUE)
  mod <- paste0(v[i], "_spatial")
  out <- stan(file="R/spatial.stan",
              data=read_rdump(paste0("data_stan/", v[i], ".Rdump")),
              pars=c("nu", "Z_", "Z_new_"), include=F,
              iter=n.iter, chains=n.chain, refresh=50, thin=20,
              sample_file=paste0("out_stan/", mod, ".csv"))
}
stopCluster(p.c)
