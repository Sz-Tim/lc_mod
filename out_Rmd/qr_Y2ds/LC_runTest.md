---
title: "Run model tests"
author: "Tim Szewczyk"
output: html_document
---





```r
grdSz <- "01_1a"
blockSize <- 10  # block = (blockSize x blockSize) grid cells

# cell-block reference tibble
cb.i <- read_csv(paste0("data/roads_", grdSz, ".csv")) %>% 
  mutate(CellRow=1:n_distinct(top) %>% rep(n_distinct(left)),
         CellCol=1:n_distinct(left) %>% rep(each=n_distinct(top))) %>%
  filter((CellRow <= max((CellRow %/% blockSize) * blockSize)) &
           (CellCol <= max((CellCol %/% blockSize) * blockSize))) %>%
  mutate(BlockRow=((CellRow-1)%/%blockSize)+1, 
         BlockCol=((CellCol-1)%/%blockSize)+1,
         BlockID=paste(BlockCol, BlockRow) %>% factor %>% as.numeric) %>%
  select(c(CellID, CellRow, CellCol, BlockID, BlockRow, BlockCol, left, top))

# covariates summarized to blocks
pop00 <- read_csv(paste0("data/pop00_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(popTot=sum(sum)) 
hous00 <- read_csv(paste0("data/housing00_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(secHome=sum(sum)) 
rdLen <- read_csv(paste0("data/roads_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(rdLen=sum(roadLen)) 
clim <- read_csv(paste0("data/clim_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% 
  summarise(b1=mean(bio1_mean), b7=mean(bio7_mean), b12=mean(bio12_mean))
pWP <- read_csv(paste0("data/pWP_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>%
  add_blocks(cb.i=cb.i) %>% summarise(mnWP=mean(mean)/100)

# land cover summarized to blocks
grnt <- read_csv(paste0("data/out_", grdSz, "_grnt.csv")) %>% 
  mutate(CellID=1:nrow(.)) %>% add_blocks(cb.i=cb.i) %>% 
  summarise(Dev=sum(V1)/n(), Oth=sum(V2)/n(), Hwd=sum(V3)/n(), 
            WP=sum(V4)/n(), Evg=sum(V5)/n(), Mxd=sum(V6)/n()) %>%
  select(-BlockID) %>% as.matrix
nlcd <- read_csv(paste0("data/out_",grdSz,"_nlcd.csv"))  %>% 
  mutate(CellID=1:nrow(.)) %>% add_blocks(cb.i=cb.i) %>% 
  summarise(Dev=sum(V1)/n(), Oth=sum(V2)/n(), Hwd=sum(V3)/n(), 
            Evg=sum(V4)/n(), Mxd=sum(V5)/n()) %>%
  select(-BlockID) %>% as.matrix
```



```r
# small scale runs: set nCell for Y1&Y2 and Y2
set.seed(2211)
nFit <- 1200
nNew <- 800
n <- sampleCells(nFit, nNew, nrow(grnt))

# Y1 & Y2
Y1.fit <- grnt[n$fit,]
Y1.new <- grnt[n$new,]
Y2 <- nlcd[n$all,]

# covariates: bias (Dev, Oth, Hwd, Evg, Mxd)
Xd <- vector("list", 4)
Xd[[1]] <- cbind(scale(rdLen$rdLen[n$all]), scale(hous00$secHome[n$all]))
Xd[[2]] <- cbind(scale(rdLen$rdLen[n$all]), scale(pop00$popTot[n$all]))
Xd[[3]] <- matrix(scale(clim$b7[n$all]), ncol=1)
Xd[[4]] <- cbind(scale(clim$b1[n$all]), scale(clim$b12[n$all]))
nBd <- map_int(Xd, ncol)  # nBeta for each covariate

# covariates: WP|Evg
Xp <- cbind(scale(pWP$mnWP[n$all]), scale(clim$b1[n$all]), scale(clim$b12[n$all]))
nBp <- ncol(Xp)
```



```r
d <- list(n1=nFit, n2=nFit+1, n3=n$tot, L=6, nB_d=nBd, nB_p=nBp,
          Y1=Y1.fit[,-6], Y2=Y2[,-5], 
          X_d1=Xd[[1]], X_d2=Xd[[2]], X_d3=Xd[[3]], X_d4=Xd[[4]], X_p=Xp)
out <- stan(file="code/LC_mod_qrBeta.stan", data=d, init=0, thin=10, 
            iter=10000, warmup=5000, chains=8, seed=4337, refresh=500,
            include=FALSE, pars=c("Y2_ds", "nu"),
            control=list(max_treedepth=30, adapt_delta=0.85, metric="diag_e"))
```

```
## Loading required namespace: rstudioapi
```

```
## Warning: There were 26 divergent transitions after warmup. Increasing adapt_delta above 0.85 may help. See
## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

```
## Warning: There were 8 chains where the estimated Bayesian Fraction of Missing Information was low. See
## http://mc-stan.org/misc/warnings.html#bfmi-low
```

```
## Warning: Examine the pairs() plot to diagnose sampling problems
```

```r
check_treedepth(out); check_energy(out); check_div(out)
```

```
## [1] "24 of 4000 iterations saturated the maximum tree depth of 10 (0.6%)"
## [1] "Run again with max_depth set to a larger value to avoid saturation"
```

```
## [1] "Chain 1: E-BFMI = 0.063811519766332"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 2: E-BFMI = 0.0288461415045991"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 3: E-BFMI = 0.0624903319940153"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 4: E-BFMI = 0.0424209706543063"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 5: E-BFMI = 0.0787610550677791"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 6: E-BFMI = 0.0418571622779001"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 7: E-BFMI = 0.0557775213619378"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
## [1] "Chain 8: E-BFMI = 0.0521845014743615"
## [1] "E-BFMI below 0.2 indicates you may need to reparameterize your model"
```

```
## [1] "26 of 4000 iterations ended with a divergence (0.65%)"
## [1] "Try running with larger adapt_delta to remove the divergences"
```

```r
sampler_params <- get_sampler_params(out, inc_warmup=FALSE)
n_gradients <- sapply(sampler_params, function(x) sum(x[,'n_leapfrog__']))
n_gradients; sum(n_gradients)
```

```
## [1] 153100 116972 193036 193036  98700 129036 110476 259852
```

```
## [1] 1254208
```



```r
##########
## munging
##########

# Full posterior
gg.nu <- ggs(out, "n_eta") %>% arrange(Parameter, Chain, Iteration)
nGG <- attr(gg.nu, "nChains")*attr(gg.nu, "nIterations")
gg.nu %<>% mutate(Y1=t(rbind(Y1.fit, Y1.new)) %>% c %>% rep(each=nGG),
                  LC=1:6 %>% rep(each=nGG) %>% rep(times=n$tot),
                  Cell=1:n$tot %>% rep(each=nGG*6),
                  Set=c("Y1+Y2", "Y2") %>% rep(times=c(nFit, nNew)*nGG*6))

# Medians
gg.med <- gg.nu %>% group_by(Cell, LC, Set, Parameter) %>%
  summarise(Y1=first(Y1), med=median(value), 
            q05=quantile(value, 0.05), q25=quantile(value, 0.25),
            q75=quantile(value, 0.75), q95=quantile(value, 0.95)) %>%
  ungroup() %>% group_by(Cell)

# Combine WP + Evg to compare to Y2
gg.EvgComb <- gg.nu
gg.EvgComb$LC[gg.EvgComb$LC==5] <- 4
gg.EvgMed <- gg.EvgComb %>% group_by(Cell, LC, Set) %>%
  summarise(med=median(value), Y1=first(Y1), 
            q05=quantile(value, 0.05), q25=quantile(value, 0.25),
            q75=quantile(value, 0.75), q95=quantile(value, 0.95)) %>%
  ungroup %>% mutate(Y2=t(Y2) %>% c)


##########
## plots
##########

ggplot(gg.EvgMed, aes(x=Y1, y=med)) + xlim(0,1) + ylim(0,1) + 
  geom_point() + geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) 
```

![plot of chunk outNu](LC_runTest/outNu-1.png)

```r
ggplot(gg.EvgMed, aes(x=Y1, y=Y2)) + xlim(0,1) + ylim(0,1) + 
  geom_point() + geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) 
```

![plot of chunk outNu](LC_runTest/outNu-2.png)

```r
ggplot(gg.EvgMed, aes(x=Y2, y=med)) + xlim(0,1) + ylim(0,1) + 
  geom_point() + geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) 
```

![plot of chunk outNu](LC_runTest/outNu-3.png)

```r
ggplot(gg.med, aes(x=Y1, y=med, ymin=q25, ymax=q75)) + xlim(0,1) + ylim(0,1) + 
  geom_pointrange(alpha=0.5, colour="dodgerblue", fatten=1.5) + 
  geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) 
```

![plot of chunk outNu](LC_runTest/outNu-4.png)

```r
ggplot(gg.EvgMed, aes(x=Y1, xend=Y1, y=Y2, yend=med,
                      colour=abs(Y2-Y1)<abs(med-Y1))) + 
  geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) +
  scale_colour_manual(values=c("darkgreen", "red")) + xlim(0,1) + ylim(0,1) +
  geom_segment(arrow=arrow(length=unit(0.1, "cm")), alpha=0.4) + 
  labs(x="Y1", y="Y2 -> median") + theme(legend.position="none")
```

![plot of chunk outNu](LC_runTest/outNu-5.png)

```r
##########
## RMSE
##########

gg.med %>% ungroup %>% group_by(Set, LC) %>%
  summarise(rmse.mod=(med-Y1)^2 %>% mean %>% sqrt %>% round(3))
```

```
## # A tibble: 12 x 3
## # Groups:   Set [?]
##      Set    LC rmse.mod
##    <chr> <int>    <dbl>
##  1 Y1+Y2     1    0.027
##  2 Y1+Y2     2    0.042
##  3 Y1+Y2     3    0.076
##  4 Y1+Y2     4    0.046
##  5 Y1+Y2     5    0.055
##  6 Y1+Y2     6    0.081
##  7    Y2     1    0.066
##  8    Y2     2    0.097
##  9    Y2     3    0.174
## 10    Y2     4    0.093
## 11    Y2     5    0.113
## 12    Y2     6    0.169
```

```r
gg.EvgMed %>% ungroup %>% group_by(Set, LC) %>%
  summarise(rmse.mod=(med-Y1)^2 %>% mean %>% sqrt %>% round(3),
            rmse.Y2=(Y2-Y1)^2 %>% mean %>% sqrt %>% round(3),
            diff=rmse.mod-rmse.Y2, prop=(diff/rmse.Y2) %>% round(3))
```

```
## # A tibble: 10 x 6
## # Groups:   Set [?]
##      Set    LC rmse.mod rmse.Y2   diff   prop
##    <chr> <dbl>    <dbl>   <dbl>  <dbl>  <dbl>
##  1 Y1+Y2     1    0.027   0.073 -0.046 -0.630
##  2 Y1+Y2     2    0.042   0.103 -0.061 -0.592
##  3 Y1+Y2     3    0.076   0.188 -0.112 -0.596
##  4 Y1+Y2     4    0.097   0.202 -0.105 -0.520
##  5 Y1+Y2     6    0.081   0.216 -0.135 -0.625
##  6    Y2     1    0.066   0.075 -0.009 -0.120
##  7    Y2     2    0.097   0.106 -0.009 -0.085
##  8    Y2     3    0.174   0.184 -0.010 -0.054
##  9    Y2     4    0.104   0.200 -0.096 -0.480
## 10    Y2     6    0.169   0.206 -0.037 -0.180
```



```r
gg.b <- ggs(out, "beta")
ggs_caterpillar(gg.b) + geom_vline(xintercept=0)
```

![plot of chunk outBeta](LC_runTest/outBeta-1.png)

```r
ggs_traceplot(gg.b) + facet_wrap(~Parameter, scales="free")
```

![plot of chunk outBeta](LC_runTest/outBeta-2.png)

```r
gg.b %>% group_by(Parameter) %>%
  summarise(q025=quantile(value, 0.025) %>% round(3), 
            q25=quantile(value, 0.25) %>% round(3),
            med=median(value) %>% round(3),
            q75=quantile(value, 0.75) %>% round(3), 
            q975=quantile(value, 0.975) %>% round(3))
```

```
## # A tibble: 10 x 6
##    Parameter   q025    q25    med    q75   q975
##       <fctr>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
##  1 beta_d[1]  0.020  0.023  0.025  0.026  0.029
##  2 beta_d[2]  0.001  0.004  0.006  0.007  0.010
##  3 beta_d[3] -0.037 -0.033 -0.031 -0.028 -0.024
##  4 beta_d[4] -0.007 -0.002  0.001  0.004  0.009
##  5 beta_d[5] -0.055 -0.049 -0.045 -0.042 -0.035
##  6 beta_d[6] -0.084 -0.080 -0.077 -0.075 -0.069
##  7 beta_d[7] -0.028 -0.023 -0.021 -0.018 -0.013
##  8 beta_p[1] -0.080 -0.006  0.034  0.073  0.150
##  9 beta_p[2]  0.578  0.681  0.740  0.794  0.902
## 10 beta_p[3] -0.418 -0.351 -0.315 -0.281 -0.215
```



```r
stan_ess(out)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk outDiag](LC_runTest/outDiag-1.png)

```r
stan_diag(out)
```

![plot of chunk outDiag](LC_runTest/outDiag-2.png)

```r
stan_rhat(out)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 22 rows containing non-finite values (stat_bin).
```

![plot of chunk outDiag](LC_runTest/outDiag-3.png)

```r
ggs_crosscorrelation(gg.b)
```

![plot of chunk outDiag](LC_runTest/outDiag-4.png)

```r
shinystan::launch_shinystan(out)
```

```
## 
## Creating shinystan object...
```

```
## 
## Launching ShinyStan interface... for large models this  may take some time.
```

```
## Loading required package: shiny
```

```
## 
## Listening on http://127.0.0.1:3341
```

```
## Warning: Removed 208 rows containing missing values (geom_bar).
```

