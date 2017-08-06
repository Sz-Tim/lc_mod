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
  add_blocks(cb.i=cb.i) %>% summarise(popTot=sum(sum)) %>% 
  ungroup %>% mutate(popSc=scale(popTot))
hous00 <- read_csv(paste0("data/housing00_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(secHome=sum(sum)) %>% 
  ungroup %>% mutate(secHomeSc=scale(secHome))
rdLen <- read_csv(paste0("data/roads_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(rdLen=sum(roadLen)) %>% 
  ungroup %>% mutate(lenSc=scale(rdLen))
clim <- read_csv(paste0("data/clim_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>%  
  summarise(b1=mean(bio1_mean), b7=mean(bio7_mean), b12=mean(bio12_mean)) %>% 
  ungroup %>% mutate(tmean=scale(b1), tseas=scale(b7), precip=scale(b12))

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
set.seed(2121)
nFit <- 1200
nNew <- 800
n <- sampleCells(nFit, nNew, nrow(grnt))

# Y1 & Y2
Y1.fit <- grnt[n$fit,]
Y1.new <- grnt[n$new,]
Y2 <- nlcd[n$all,]

# covariates: bias (Dev, Oth, Hwd, Evg, Mxd)
Xd <- vector("list", 4)
Xd[[1]] <- cbind(rdLen$lenSc[n$all], hous00$secHomeSc[n$all], clim$tmean[n$all])
Xd[[2]] <- cbind(rdLen$lenSc[n$all], pop00$popSc[n$all])
Xd[[3]] <- matrix(clim$tseas[n$all], ncol=1)
Xd[[4]] <- cbind(clim$tmean[n$all], clim$precip[n$all], rdLen$lenSc[n$all])
nBd <- map_int(Xd, ncol)  # nBeta for each covariate

# covariates: WP|Evg
Xp <- cbind(clim$tmean[n$all], clim$precip[n$all])
nBp <- ncol(Xp)
```



```r
d <- list(n1=nFit, n2=nFit+1, n3=n$tot, L=6, nB_d=nBd, nB_p=nBp,
          Y1=Y1.fit[,-6], Y2=Y2[,-5], 
          X_d1=Xd[[1]], X_d2=Xd[[2]], X_d3=Xd[[3]], X_d4=Xd[[4]], X_p=Xp)
out <- stan(file="code/LC_mod.stan", data=d, init=0, thin=50, 
            iter=20000, warmup=15000, chains=8, seed=43337, 
            control=list(max_treedepth=15))
```

```
## Loading required namespace: rstudioapi
```

```
## Warning: There were 7 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
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
##  1 Y1+Y2     1    0.045
##  2 Y1+Y2     2    0.067
##  3 Y1+Y2     3    0.078
##  4 Y1+Y2     4    0.046
##  5 Y1+Y2     5    0.059
##  6 Y1+Y2     6    0.088
##  7    Y2     1    0.064
##  8    Y2     2    0.089
##  9    Y2     3    0.177
## 10    Y2     4    0.096
## 11    Y2     5    0.116
## 12    Y2     6    0.173
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
##  1 Y1+Y2     1    0.045   0.073 -0.028 -0.384
##  2 Y1+Y2     2    0.067   0.109 -0.042 -0.385
##  3 Y1+Y2     3    0.078   0.182 -0.104 -0.571
##  4 Y1+Y2     4    0.093   0.194 -0.101 -0.521
##  5 Y1+Y2     6    0.088   0.211 -0.123 -0.583
##  6    Y2     1    0.064   0.073 -0.009 -0.123
##  7    Y2     2    0.089   0.098 -0.009 -0.092
##  8    Y2     3    0.177   0.186 -0.009 -0.048
##  9    Y2     4    0.105   0.210 -0.105 -0.500
## 10    Y2     6    0.173   0.212 -0.039 -0.184
```



```r
gg.b <- ggs(out, "beta")
ggs_caterpillar(gg.b) + geom_vline(xintercept=0)
```

![plot of chunk outBeta](LC_runTest/outBeta-1.png)

```r
gg.b %>% group_by(Parameter) %>%
  summarise(q025=quantile(value, 0.025) %>% round(3), 
            q25=quantile(value, 0.25) %>% round(3),
            med=median(value) %>% round(3),
            q75=quantile(value, 0.75) %>% round(3), 
            q975=quantile(value, 0.975) %>% round(3))
```

```
## # A tibble: 11 x 6
##     Parameter   q025    q25    med    q75   q975
##        <fctr>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
##  1 beta_d1[1]  0.013  0.016  0.018  0.019  0.022
##  2 beta_d1[2] -0.001  0.001  0.003  0.004  0.006
##  3 beta_d1[3]  0.009  0.011  0.013  0.014  0.017
##  4 beta_d2[1] -0.038 -0.033 -0.031 -0.028 -0.024
##  5 beta_d2[2] -0.004  0.000  0.001  0.004  0.008
##  6 beta_d3[1] -0.055 -0.050 -0.046 -0.043 -0.036
##  7 beta_d4[1] -0.080 -0.075 -0.072 -0.069 -0.064
##  8 beta_d4[2] -0.024 -0.019 -0.017 -0.014 -0.009
##  9 beta_d4[3]  0.002  0.006  0.008  0.011  0.015
## 10  beta_p[1]  0.645  0.717  0.760  0.796  0.880
## 11  beta_p[2] -0.379 -0.319 -0.288 -0.251 -0.185
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
## Warning: Removed 42 rows containing non-finite values (stat_bin).
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
## Listening on http://127.0.0.1:4829
```

