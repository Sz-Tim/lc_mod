# This script analyzes the results of LC_processVS_hpc.R

Packages <- c("rstan", "coda", "ggmcmc", "bayesplot", "tidyverse", 
              "loo", "purrr")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

setwd("C:/Users/tms1044/Desktop/lc_mod")

ggs_autocor <- function (D, family = NA, nLags = 50, greek = FALSE) 
{
  if (!is.na(family)) {
    D <- get_family(D, family = family)
  }
  nIter <- attr(D, "nIteration")
  if (nIter < nLags) {
    warning(sprintf("nLags=%d is larger than number of iterations, computing until max possible lag %d", 
                    nLags, nIter))
    nLags <- nIter
  }
  wc.ac <- D %>% dplyr::group_by(Parameter, Chain) %>% dplyr::do(ac(.$value, 
                                                                    nLags))
  f <- ggplot(wc.ac, aes(x=Lag, y=Autocorrelation, colour=as.factor(Chain))) + 
    theme(legend.position="none") + geom_hline(yintercept=0, colour="gray70") +
    geom_point(size=0.5) + geom_line(stat="identity", position="identity") + 
    scale_colour_brewer(type="qual") + facet_wrap(~Parameter) + ylim(-1, 1)
  return(f)
}

thin.seq <- seq(4,20,2)
chn.dir <- "out/"
sum.dir <- "thin_summaries/"
mods <- c("Cens", "Clim", "ClimCens")
for(m in 1:3) {
  mod <- mods[m]
  
  # GGMCMC
  beta.fls <- list.files(sum.dir, pattern="out_betas_", full.names=TRUE)
  beta.fls <- beta.fls[grepl(paste0(mod,"."), beta.fls, fixed=TRUE)]
  names(beta.fls) <- beta.fls
  beta.fls <- beta.fls[c(7:9,1:6)]
  map(beta.fls, readRDS) %>% 
    map2(., beta.fls, ~(ggs_density(ggs(.x)) + ggtitle(.y) +
                    facet_wrap(~Parameter, scales="free"))) %>% print
  map(beta.fls, readRDS) %>% 
    map2(., beta.fls, ~(ggs_traceplot(ggs(.x)) + ggtitle(.y) + 
                    facet_wrap(~Parameter, scales="free"))) %>% print
  map(beta.fls, readRDS) %>% 
    map2(., beta.fls, ~(ggs_autocor(ggs(.x)) + ggtitle(.y))) %>% print
  map(beta.fls, readRDS) %>% 
    map2(., beta.fls, ~(ggs_caterpillar(ggs(.x)) + ggtitle(.y))) %>% print
  map(beta.fls, readRDS) %>% 
    map2(., beta.fls, ~(ggs_geweke(ggs(.x)) + ggtitle(.y))) %>% print
  map(beta.fls, readRDS) %>% 
    map2(., beta.fls, ~(ggs_Rhat(ggs(.x)) + ggtitle(.y))) %>% print
  
  
  # effective sample size
  par(mfrow=c(3,3))
  map(beta.fls, readRDS) %>% 
    walk2(., thin.seq, ~hist(effectiveSize(.x), main=.y, xlab="n_eff"))
  par(mfrow=c(3,3))
  map(beta.fls, readRDS) %>% 
    walk2(., thin.seq, ~hist(effectiveSize(.x)/(5000*6/.y), main=.y,
                  xlab="n_eff/n_iter", xlim=c(0,1.5)))
  par(mfrow=c(1,1))
  
}
