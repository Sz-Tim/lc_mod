---
title: "Run model tests"
author: "Tim Szewczyk"
output: html_document
---





```r
grdSz <- "01_1a"
blockSize <- 20  # block = (blockSize x blockSize) grid cells

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
  add_blocks(cb.i=cb.i) %>% summarise(popTot=log(sum(sum)+0.001))
hous00 <- read_csv(paste0("data/housing00_", grdSz, ".csv")) %>% 
  rename(CellID=category) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(secHome=log(sum(sum)+0.001))
rdLen <- read_csv(paste0("data/roads_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% summarise(rdLen=log(sum(roadLen)+0.001)) 
clim <- read_csv(paste0("data/clim_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% 
  summarise(b1=mean(bio1_mean), b7=mean(bio7_mean), b12=mean(bio12_mean))
topo <- read_csv(paste0("data/topo_", grdSz, ".csv")) %>% 
  add_blocks(cb.i=cb.i) %>% 
  summarise(el=mean(el_mean), rugg=mean(rugg_mean))
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
set.seed(2222)
nFit <- 324
nNew <- 216
n <- sampleCells(nFit, nNew, nrow(grnt))

# Y1 & Y2
Y1.fit <- grnt[n$fit,]
Y1.new <- grnt[n$new,]
Y2 <- nlcd[n$all,]

# covariates: bias (Dev, Oth, Hwd, Evg, Mxd)
Xd <- vector("list", 4)
Xd[[1]] <- cbind(scale(rdLen$rdLen[n$all]), 
                 scale(pop00$popTot[n$all]),
                 scale(topo$el[n$all]))
Xd[[2]] <- cbind(scale(rdLen$rdLen[n$all]), 
                 scale(pop00$popTot[n$all]),
                 scale(clim$b7[n$all]))
Xd[[3]] <- cbind(scale(clim$b7[n$all]), 
                 scale(pop00$popTot[n$all]), 
                 scale(topo$el[n$all]))
Xd[[4]] <- cbind(scale(clim$b1[n$all]), 
                 scale(clim$b12[n$all]),
                 scale(pop00$popTot[n$all]))
nBd <- map_int(Xd, ncol)  # nBeta for each covariate

# covariates: WP|Evg
Xp <- cbind(scale(pWP$mnWP[n$all]), 
            scale(clim$b1[n$all]),
            scale(pop00$popTot[n$all]),
            scale(topo$rugg[n$all]))
nBp <- ncol(Xp)

Yd <- tibble(d1=c(scale(grnt[,1]-nlcd[,1])),
             d2=c(scale(grnt[,2]-nlcd[,2])),
             d3=c(scale(grnt[,3]-nlcd[,3])),
             d4=c(scale((grnt[,4] + grnt[,5])-nlcd[,4])),
             nuWP=c(scale((grnt[,4]+0.0001)/(grnt[,4] + grnt[,5] + 0.0001))),
             valWP=c(scale(pWP$mnWP)),
             rdLen=c(scale(rdLen$rdLen)),
             pop00=c(scale(pop00$popTot)),
             hous00=c(scale(hous00$secHome)),
             tmean=c(scale(clim$b1)),
             tseas=c(scale(clim$b7)),
             precip=c(scale(clim$b12)),
             el=c(scale(topo$el)),
             rugg=c(scale(topo$rugg)))
cor(Yd)
```

```
##                 d1          d2          d3          d4       nuWP
## d1      1.00000000 -0.33882794 -0.07262417 -0.19912831  0.4876080
## d2     -0.33882794  1.00000000 -0.00454309  0.01113119 -0.3094890
## d3     -0.07262417 -0.00454309  1.00000000 -0.13472275 -0.1364609
## d4     -0.19912831  0.01113119 -0.13472275  1.00000000 -0.4075154
## nuWP    0.48760804 -0.30948902 -0.13646091 -0.40751538  1.0000000
## valWP   0.47481480 -0.22331345 -0.28958844 -0.47015834  0.6797474
## rdLen   0.39382035 -0.31451645 -0.02512107 -0.26627295  0.3858912
## pop00   0.49752105 -0.30226466 -0.18959821 -0.37157838  0.5797619
## hous00  0.35997008 -0.19173991 -0.07044672 -0.31146868  0.4430759
## tmean   0.40957253 -0.22040829 -0.22428348 -0.54978630  0.6375672
## tseas   0.40317102 -0.23346639 -0.30446337 -0.44314581  0.6064141
## precip -0.25880873  0.13201364 -0.05354982  0.23043485 -0.3590041
## el     -0.43763285  0.22072614  0.24923787  0.53843591 -0.6517388
## rugg   -0.38363829  0.09114154  0.09613058  0.07773937 -0.3979684
##             valWP       rdLen      pop00      hous00      tmean      tseas
## d1      0.4748148  0.39382035  0.4975210  0.35997008  0.4095725  0.4031710
## d2     -0.2233134 -0.31451645 -0.3022647 -0.19173991 -0.2204083 -0.2334664
## d3     -0.2895884 -0.02512107 -0.1895982 -0.07044672 -0.2242835 -0.3044634
## d4     -0.4701583 -0.26627295 -0.3715784 -0.31146868 -0.5497863 -0.4431458
## nuWP    0.6797474  0.38589119  0.5797619  0.44307591  0.6375672  0.6064141
## valWP   1.0000000  0.46844349  0.6433010  0.47972409  0.8394792  0.7872503
## rdLen   0.4684435  1.00000000  0.5154575  0.45647512  0.4411659  0.3243775
## pop00   0.6433010  0.51545749  1.0000000  0.73980013  0.6158805  0.5403887
## hous00  0.4797241  0.45647512  0.7398001  1.00000000  0.4543204  0.3680021
## tmean   0.8394792  0.44116588  0.6158805  0.45432036  1.0000000  0.8141676
## tseas   0.7872503  0.32437753  0.5403887  0.36800212  0.8141676  1.0000000
## precip -0.3508613 -0.27356829 -0.3689006 -0.27441598 -0.5489150 -0.5083765
## el     -0.8748212 -0.46044833 -0.6175100 -0.44375580 -0.9801083 -0.8733614
## rugg   -0.4926853 -0.52938015 -0.4523009 -0.35244912 -0.4401005 -0.2597428
##             precip         el        rugg
## d1     -0.25880873 -0.4376328 -0.38363829
## d2      0.13201364  0.2207261  0.09114154
## d3     -0.05354982  0.2492379  0.09613058
## d4      0.23043485  0.5384359  0.07773937
## nuWP   -0.35900413 -0.6517388 -0.39796841
## valWP  -0.35086126 -0.8748212 -0.49268526
## rdLen  -0.27356829 -0.4604483 -0.52938015
## pop00  -0.36890055 -0.6175100 -0.45230085
## hous00 -0.27441598 -0.4437558 -0.35244912
## tmean  -0.54891504 -0.9801083 -0.44010050
## tseas  -0.50837652 -0.8733614 -0.25974280
## precip  1.00000000  0.4947391  0.39044291
## el      0.49473913  1.0000000  0.44505949
## rugg    0.39044291  0.4450595  1.00000000
```



```r
# block distances & adjacency
b.rows <- cb.i$BlockRow[match(rdLen$BlockID[n$all], cb.i$BlockID)]
b.cols <- cb.i$BlockCol[match(rdLen$BlockID[n$all], cb.i$BlockID)]
coords <- data.frame(b.rows, b.cols)
D <- as.matrix(dist(coords))
W <- as.matrix(dist(coords, diag=T, upper=T)) < 1.5
diag(W) <- 0

# knots
m.x <- 5
m.y <- 6
m <- m.x * m.y
coords_star <- place_knots(m.x, m.y, coords)
D_star <- as.matrix(dist(coords_star))
D_site_star <- as.matrix(dist(rbind(coords, coords_star)))[1:n$tot, 
                                                           (n$tot+1):(n$tot+m)]
```



```r
d <- list(n1=nFit, n2=nFit+1, n3=n$tot, L=6, nB_d=nBd, nB_p=nBp,
          Y1=Y1.fit[,-6], Y2=Y2[,-5], W=W, W_n=sum(W)/2,
          m=m, D_star=D_star, D_site_star=D_site_star,
          X_d1=Xd[[1]], X_d2=Xd[[2]], X_d3=Xd[[3]], X_d4=Xd[[4]], X_p=Xp)
stan_rdump(ls(d), file="code/LC_mod_examp.Rdump", envir=list2env(d))
out <- stan(file="code/LC_mod_Y2X_GPP.stan", init=0, thin=25,
            data=read_rdump("code/LC_mod_examp.Rdump"), 
            iter=5000, warmup=2500, chains=8, seed=4337, refresh=500,
            include=FALSE, pars=c("Y2_ds", "w_z", "e_z", "w", "sigma_e_tilde",
                                  "Cstar", "w_star", "inv_Cstar",
                                  "C_site_star", "C_ss_inv_Cstar"),
            control=list(max_treedepth=15, adapt_delta=0.9))
```

```
## hash mismatch so recompiling; make sure Stan code ends with a blank line
```

```
## In file included from C:/Users/tms1044/Documents/R/win-library/3.4/BH/include/boost/config.hpp:39:0,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/BH/include/boost/math/tools/config.hpp:13,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/core/var.hpp:7,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/core/gevv_vvv_vari.hpp:5,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/core.hpp:12,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/mat.hpp:4,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/StanHeaders/include/stan/math.hpp:4,
##                  from C:/Users/tms1044/Documents/R/win-library/3.4/StanHeaders/include/src/stan/model/model_header.hpp:4,
##                  from file6ca431e3001.cpp:8:
## C:/Users/tms1044/Documents/R/win-library/3.4/BH/include/boost/config/compiler/gcc.hpp:186:0: warning: "BOOST_NO_CXX11_RVALUE_REFERENCES" redefined
##  #  define BOOST_NO_CXX11_RVALUE_REFERENCES
##  ^
## <command-line>:0:0: note: this is the location of the previous definition
## cc1plus.exe: warning: unrecognized command line option "-Wno-ignored-attributes"
```

```
## Warning: There were 41 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help. See
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
## [1] "0 of 800 iterations saturated the maximum tree depth of 10 (0%)"
```

```
## [1] "41 of 800 iterations ended with a divergence (5.125%)"
## [1] "Try running with larger adapt_delta to remove the divergences"
```

```r
sampler_params <- get_sampler_params(out, inc_warmup=FALSE)
n_gradients <- sapply(sampler_params, function(x) sum(x[,'n_leapfrog__']))
n_gradients; sum(n_gradients)
```

```
## [1] 21788 12325 11843 22231 12700 24672 24440 51100
```

```
## [1] 181099
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
                  BlockID=n$all %>% rep(each=nGG*6),
                  CellID=1:length(n$all) %>% rep(each=nGG*6),
                  Set=c("Y1+Y2", "Y2") %>% rep(times=c(nFit, nNew)*nGG*6)) %>%
  mutate(BlockRow=cb.i$BlockRow[match(.$BlockID, cb.i$BlockID)], 
         BlockCol=cb.i$BlockCol[match(.$BlockID, cb.i$BlockID)])

# Medians
gg.med <- gg.nu %>% 
  group_by(CellID, BlockID, BlockRow, BlockCol, LC, Set, Parameter) %>%
  summarise(Y1=first(Y1), med=median(value), 
            q05=quantile(value, 0.05), q25=quantile(value, 0.25),
            q75=quantile(value, 0.75), q95=quantile(value, 0.95)) %>%
  ungroup() %>% group_by(BlockID)

# Combine WP + Evg to compare to Y2
gg.EvgComb <- gg.nu
gg.EvgComb$LC[gg.EvgComb$LC==5] <- 4
gg.EvgMed <- gg.EvgComb %>% group_by(CellID, BlockID, LC, Set) %>%
  summarise(med=median(value), Y1=first(Y1), 
            q05=quantile(value, 0.05), q25=quantile(value, 0.25),
            q75=quantile(value, 0.75), q95=quantile(value, 0.95)) %>%
  arrange(CellID, LC) %>%
  ungroup %>% mutate(Y2=t(Y2) %>% c)


##########
## plots
##########

ggplot(gg.EvgMed, aes(x=Y1, y=med)) + xlim(0,1) + ylim(0,1) + 
  geom_point(alpha=0.5) + facet_grid(Set~LC) + 
  geom_abline(slope=1, linetype=2, colour="red") + ggtitle("Y2_X")
```

![plot of chunk outNu](LC_runTest/outNu-1.png)

```r
ggplot(gg.EvgMed, aes(x=Y1, y=Y2)) + xlim(0,1) + ylim(0,1) + 
  geom_point(alpha=0.5) + facet_grid(Set~LC) + 
  geom_abline(slope=1, linetype=2, colour="red") + ggtitle("Y2_X") 
```

![plot of chunk outNu](LC_runTest/outNu-2.png)

```r
ggplot(gg.EvgMed, aes(x=Y2, y=med)) + xlim(0,1) + ylim(0,1) + 
  geom_point(alpha=0.5) + facet_grid(Set~LC) + 
  geom_abline(slope=1, linetype=2, colour="red")  + ggtitle("Y2_X")
```

![plot of chunk outNu](LC_runTest/outNu-3.png)

```r
ggplot(gg.med, aes(x=Y1, y=med, ymin=q25, ymax=q75)) + xlim(0,1) + ylim(0,1) + 
  geom_pointrange(alpha=0.5, colour="dodgerblue", fatten=1.5) + 
  geom_abline(slope=1, linetype=3) + facet_grid(Set~LC)  + ggtitle("Y2_X")
```

![plot of chunk outNu](LC_runTest/outNu-4.png)

```r
ggplot(gg.EvgMed, aes(x=Y1, xend=Y1, y=Y2, yend=med,
                      colour=abs(Y2-Y1)<abs(med-Y1))) + 
  geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) +
  scale_colour_manual(values=c("darkgreen", "red")) + xlim(0,1) + ylim(0,1) +
  geom_segment(arrow=arrow(length=unit(0.1, "cm")), alpha=0.4) + 
  labs(x="Y1", y="Y2 -> median") + theme(legend.position="none") +
  ggtitle("Y2_X")
```

![plot of chunk outNu](LC_runTest/outNu-5.png)

```r
ggplot(gg.med, aes(x=BlockCol, y=BlockRow, fill=med-Y1)) + 
  geom_tile() + facet_grid(Set~LC) +
  scale_fill_gradient2() + ggtitle("Y2_X")
```

![plot of chunk outNu](LC_runTest/outNu-6.png)

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
##  1 Y1+Y2     1    0.025
##  2 Y1+Y2     2    0.052
##  3 Y1+Y2     3    0.061
##  4 Y1+Y2     4    0.040
##  5 Y1+Y2     5    0.052
##  6 Y1+Y2     6    0.056
##  7    Y2     1    0.033
##  8    Y2     2    0.055
##  9    Y2     3    0.118
## 10    Y2     4    0.057
## 11    Y2     5    0.070
## 12    Y2     6    0.119
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
##  1 Y1+Y2     1    0.025   0.063 -0.038 -0.603
##  2 Y1+Y2     2    0.052   0.078 -0.026 -0.333
##  3 Y1+Y2     3    0.061   0.142 -0.081 -0.570
##  4 Y1+Y2     4    0.075   0.174 -0.099 -0.569
##  5 Y1+Y2     6    0.056   0.174 -0.118 -0.678
##  6    Y2     1    0.033   0.051 -0.018 -0.353
##  7    Y2     2    0.055   0.069 -0.014 -0.203
##  8    Y2     3    0.118   0.150 -0.032 -0.213
##  9    Y2     4    0.087   0.172 -0.085 -0.494
## 10    Y2     6    0.119   0.166 -0.047 -0.283
```



```r
gg.b <- ggs(out, "beta")
ggs_caterpillar(gg.b) + geom_vline(xintercept=0)
```

![plot of chunk outBeta](LC_runTest/outBeta-1.png)

```r
traceplot(out, pars=c("beta_d", "beta_p"))
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
## # A tibble: 16 x 6
##     Parameter   q025    q25    med    q75   q975
##        <fctr>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
##  1  beta_d[1]  0.002  0.006  0.008  0.009  0.013
##  2  beta_d[2]  0.008  0.014  0.016  0.018  0.022
##  3  beta_d[3] -0.027 -0.021 -0.017 -0.014 -0.008
##  4  beta_d[4] -0.025 -0.019 -0.016 -0.013 -0.007
##  5  beta_d[5] -0.022 -0.015 -0.012 -0.009 -0.002
##  6  beta_d[6] -0.006  0.003  0.008  0.013  0.024
##  7  beta_d[7] -0.100 -0.059 -0.041 -0.024  0.009
##  8  beta_d[8] -0.016 -0.004  0.003  0.009  0.022
##  9  beta_d[9] -0.046 -0.004  0.015  0.031  0.063
## 10 beta_d[10] -0.050 -0.032 -0.023 -0.013  0.004
## 11 beta_d[11] -0.024 -0.012 -0.007 -0.001  0.009
## 12 beta_d[12] -0.018 -0.010 -0.006 -0.001  0.006
## 13  beta_p[1] -0.136 -0.018  0.052  0.127  0.265
## 14  beta_p[2]  0.189  0.391  0.491  0.596  0.789
## 15  beta_p[3]  0.160  0.276  0.342  0.415  0.557
## 16  beta_p[4] -0.280 -0.197 -0.156 -0.114 -0.035
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
## Warning: Removed 11 rows containing non-finite values (stat_bin).
```

![plot of chunk outDiag](LC_runTest/outDiag-3.png)

```r
ggs_crosscorrelation(gg.b)
```

![plot of chunk outDiag](LC_runTest/outDiag-4.png)

```r
ggs_crosscorrelation(ggs(out, "theta"))
```

![plot of chunk outDiag](LC_runTest/outDiag-5.png)

```r
sampler_params <- get_sampler_params(out, inc_warmup=FALSE) %>% do.call(rbind, .)
samp.out <- cbind(sampler_params[,c(1,6)], extract(out, pars="lp__")[[1]], 
                  extract(out, pars="L_sigma_unif")[[1]])
colnames(samp.out) <- c("accept_stat__", "energy__", "lp__", 
                        "L_sig1[1]", "L_sig1[2]", "L_sig1[3]", 
                        "L_sig1[4]", "L_sig1[5]")
pairs(samp.out, diag.panel=panel.hist, lower.panel=panel.cor,
      upper.panel=function(...) smoothScatter(...,nrpoints=0, add=TRUE))
```

![plot of chunk outDiag](LC_runTest/outDiag-6.png)

```r
samp.out <- cbind(sampler_params[,c(1,6)], extract(out, pars="lp__")[[1]], 
                  extract(out, pars="beta_d")[[1]],
                  extract(out, pars="beta_p")[[1]])
pairs(samp.out, diag.panel=panel.hist, lower.panel=panel.cor,
      upper.panel=function(...) smoothScatter(...,nrpoints=0, add=TRUE))
```

![plot of chunk outDiag](LC_runTest/outDiag-7.png)

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
## 
## Listening on http://127.0.0.1:6415
```


