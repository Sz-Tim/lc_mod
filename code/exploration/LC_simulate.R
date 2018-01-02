library(MASS); library(sevcheck)
tr_gjam_inv <- function(w, a=0.99) {
  eta <- w[-length(w)]
  eta[eta<0] <- 0
  w.p <- sum( (eta > 0 & eta < 1)*eta + (eta > 1) )
  
  if(w.p >= a) {
    while(sum(eta) > 1) {
      D.i <- (w.p^(-1)) * (1 - ((1-a)^(w.p/a)))
      eta <- D.i*eta
    }
  }
  c(eta, 1-sum(eta)) 
}





beta.md <- gg.b %>% group_by(Parameter) %>%
  summarize(med=median(value)*100)
cov.md <- ggs(out, "L_") %>% group_by(Parameter) %>%
  summarize(med=median(value))
L_Omega_1 <- matrix(cov.md$med[1:25], byrow = TRUE, nrow=5)
L_Sigma_1 <- L_Omega_1 %*% t(L_Omega_1)
L_Omega_2 <- matrix(cov.md$med[26:50], byrow = TRUE, nrow=5)
L_Sigma_2 <- L_Omega_2 %*% t(L_Omega_2) 


nu.sim <-rbind(Y1.fit, Y1.new)
Y1.sim <- matrix(nrow=nrow(nu.sim), ncol=6)
Y2_.sim <- Y1.sim
for(i in 1:nrow(nu.sim)) {
  Y1.sim[i,1:5] <- mvrnorm(1, nu.sim[i,1:5], L_Sigma_1)
  Y2_.sim[i,1:5] <- mvrnorm(1, nu.sim[i,1:5], L_Sigma_2)
}

bias <- list(4)
bias[[1]] <- Xd[[1]] %*% beta.md$med[1:3]
bias[[2]] <- Xd[[2]] %*% beta.md$med[4:6]
bias[[3]] <- Xd[[3]] %*% beta.md$med[7:9]
bias[[4]] <- Xd[[4]] %*% beta.md$med[10:12]

Y2.sim <- matrix(nrow=nrow(nu.sim), ncol=5)
Y2.sim[,1] <- Y2_.sim[,1] + bias[[1]]
Y2.sim[,2] <- Y2_.sim[,2] + bias[[2]]
Y2.sim[,3] <- Y2_.sim[,3] + bias[[3]]
Y2.sim[,4] <- (Y2_.sim[,4] + Y2_.sim[,5]) + bias[[4]]

Y1.sim.eta <- t(apply(Y1.sim, 1, tr_gjam_inv))
Y2.sim.eta <- t(apply(Y2.sim, 1, tr_gjam_inv))



d <- list(n1=nFit, n2=nFit+1, n3=n$tot, L=6, nB_d=nBd, nB_p=nBp,
          Y1=Y1.sim.eta[1:nFit,-6], Y2=Y2.sim.eta[,-5], 
          X_d1=Xd[[1]], X_d2=Xd[[2]], X_d3=Xd[[3]], X_d4=Xd[[4]], X_p=Xp)
out <- stan(file="code/exploration/sc_theta/latent.stan", init=0, thin=20,
            data=d, 
            iter=6000, warmup=4000, chains=6, seed=43337, refresh=100,
            sample_file="out/sim_test.csv",
            include=FALSE, pars=c("Y2_", "Y2new_", "nu"))


gg.nu <- ggs(out, "n_eta") %>% arrange(Parameter, Chain, Iteration)
nGG <- attr(gg.nu, "nChains")*attr(gg.nu, "nIterations")
gg.nu %<>% mutate(nu=t(nu.sim) %>% c %>% rep(each=nGG),
                  Y1=t(Y1.sim.eta) %>% c %>% rep(each=nGG),
                  LC=1:6 %>% rep(each=nGG) %>% rep(times=n$tot),
                  BlockID=pop00$BlockID[n$all] %>% rep(each=nGG*6),
                  CellID=1:length(n$all) %>% rep(each=nGG*6),
                  Set=c("Y1+Y2", "Y2") %>% rep(times=c(nFit, nNew)*nGG*6)) %>%
  mutate(BlockRow=cb.i$BlockRow[match(.$BlockID, cb.i$BlockID)], 
         BlockCol=cb.i$BlockCol[match(.$BlockID, cb.i$BlockID)])
gg.med <- gg.nu %>% 
  group_by(CellID, BlockID, BlockRow, BlockCol, LC, Set, Parameter) %>%
  summarise(nu=first(nu), Y1=first(Y1), med=median(value), 
            q05=quantile(value, 0.05), q25=quantile(value, 0.25),
            q75=quantile(value, 0.75), q95=quantile(value, 0.95)) %>%
  ungroup()
gg.med$LC <- factor(gg.med$LC, labels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))

gg.EvgComb <- gg.nu
gg.EvgComb$LC[gg.EvgComb$LC==5] <- 4
gg.EvgComb %<>% group_by(Iteration, Chain, CellID, BlockID, BlockRow, BlockCol, LC) %>%
  mutate(Y1=sum(Y1), nu=sum(nu), med=sum(value)) %>% ungroup
gg.EvgMed <- gg.EvgComb %>% group_by(CellID, BlockID, BlockRow, BlockCol, LC, Set) %>%
  summarise(med=median(value), nu=median(nu), Y1=median(Y1), 
            q05=quantile(value, 0.05), q25=quantile(value, 0.25),
            q75=quantile(value, 0.75), q95=quantile(value, 0.95)) %>%
  arrange(CellID, LC) %>%
  ungroup %>% mutate(Y2=t(Y2.sim.eta) %>% c)
gg.EvgMed$LC <- factor(gg.EvgMed$LC, labels=c("Opn", "Oth", "Dec", "Evg", "Mxd"))



LC.cols <- c("Opn"="#B51700", "Oth"="#5E5E5E", "Dec"="#19651A",
             "WP"="#B5447F", "Evg"="#0E390F", "Mxd"="#5F8B25")
LC.labs <- c("Opn"="Open Invasible", "Oth"="Other", "Dec"="Deciduous",
             "WP"="White Pine", "Evg"="Evergreen", "Mxd"="Mixed",
             "Y1+Y2"="Y1 & Y2", "Y2"="Y2 only")
mod.cols <- c("Posterior median"="#b2182b", "Y1"="#4393c3", "Y2"="#053061")


med.p <- ggplot(gg.med, aes(x=nu, y=med, ymin=q25, ymax=q75, colour=LC)) + 
  scale_x_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_y_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_colour_manual(values=LC.cols) +
  geom_pointrange(alpha=0.1, fatten=1.5) + 
  geom_abline(slope=1, linetype=3) + 
  facet_grid(.~LC, labeller=as_labeller(LC.labs)) + 
  theme(legend.position="none", axis.text=element_text(size=16),
        strip.text = element_text(size=16))
ggsave("out/sim_med.jpg", med.p, width=14, height=3, units="in")

med.p <- ggplot(gg.med, aes(x=nu, y=med, colour=LC)) + 
  scale_x_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_y_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_colour_manual(values=LC.cols) +
  geom_point(alpha=0.1) + 
  geom_abline(slope=1, linetype=3) + 
  facet_grid(.~LC, labeller=as_labeller(LC.labs)) + 
  theme(legend.position="none", axis.text=element_text(size=16),
        strip.text = element_text(size=16))
ggsave("out/sim_med_pts.jpg", med.p, width=14, height=3, units="in")

y1.p <- ggplot(gg.med, aes(x=nu, y=Y1, colour=LC)) + 
  scale_x_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_y_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_colour_manual(values=LC.cols) +
  geom_point(alpha=0.1) + 
  geom_abline(slope=1, linetype=3) + 
  facet_grid(.~LC, labeller=as_labeller(LC.labs)) +
  theme(legend.position="none", axis.text=element_text(size=16),
        strip.text = element_text(size=16))
ggsave("out/sim_Y1.jpg", y1.p, width=14, height=3, units="in")

y2.p <- ggplot(gg.EvgMed, aes(x=nu, y=Y2, colour=LC)) + 
  scale_x_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_y_continuous("", breaks=c(0, 0.5, 1), limits=c(0,1),
                     labels=c("0", "0.5", "1")) +
  scale_colour_manual(values=LC.cols) +
  geom_point(alpha=0.1) + 
  geom_abline(slope=1, linetype=3) + 
  facet_grid(.~LC, labeller=as_labeller(LC.labs)) +
  theme(legend.position="none", axis.text=element_text(size=16),
        strip.text = element_text(size=16))
ggsave("out/sim_Y2.jpg", y2.p, width=11.667, height=3, units="in")


med.box <- gg.med %>% dplyr::select(c(1:6,8:10)) %>% 
  gather(Predictor, val, 8:9)
evgMed.box <- gg.EvgMed %>% dplyr::select(c(1:9, 14)) %>% 
  gather(Predictor, val, c(7,9,10))
med.box <- rbind(med.box, filter(evgMed.box, Predictor=="Y2"))
med.box$Predictor <- factor(med.box$Predictor, levels=c("Y2", "Y1", "med"),
                               labels=c("Y2", "Y1", "Posterior median"))

ggplot(med.box, aes(x=val-nu, colour=Predictor, fill=Predictor)) + 
  geom_density(alpha=0.2, size=1) + 
  facet_wrap(~LC, scales="free_y", labeller=as_labeller(LC.labs)) + 
  scale_colour_manual(values=mod.cols) + scale_fill_manual(values=mod.cols)
dens.p <- ggplot(med.box, aes(x=val-nu, colour=Predictor, fill=Predictor)) + 
  geom_vline(xintercept=0, linetype=3, colour="gray30") +
  geom_density(alpha=0.15, size=1) + 
  scale_x_continuous("", breaks=c(-1, 0, 1), limits=c(-1,1),
                     labels=c("-1", "0", "1")) + ylab("") +
  facet_grid(Set~LC, scales="free_y", labeller=as_labeller(LC.labs)) + 
  scale_colour_manual(values=mod.cols) + scale_fill_manual(values=mod.cols) +
  theme(legend.position=c("none"),#(0.92, 0.34),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        strip.text = element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=18),
        panel.grid = element_blank())
ggsave("out/sim_dens_set.pdf", dens.p, width=14, height=5, units="in")


ggplot(gg.EvgMed, aes(x=nu, xend=nu, y=Y2, yend=med,
                      colour=abs(Y2-nu)<abs(med-nu))) + 
  geom_abline(slope=1, linetype=3) + facet_grid(Set~LC) +
  scale_colour_manual(values=c("darkgreen", "red")) + xlim(0,1) + ylim(0,1) +
  geom_segment(arrow=arrow(length=unit(0.1, "cm")), alpha=0.4) + 
  labs(x="nu", y="Y2 -> median") + theme(legend.position="none")




gg.med %>% ungroup %>% group_by( LC) %>%
  summarise(rmse.mod=(med-nu)^2 %>% mean %>% sqrt %>% round(3),
            rmse.Y1=(Y1-nu)^2 %>% mean %>% sqrt %>% round(3),
            diff=rmse.mod-rmse.Y1, prop=(diff/rmse.Y1) %>% round(3))

gg.EvgMed %>% ungroup %>% group_by(LC) %>%
  summarise(rmse.mod=(med-nu)^2 %>% mean %>% sqrt %>% round(3),
            rmse.Y2=(Y2-nu)^2 %>% mean %>% sqrt %>% round(3),
            diff=rmse.mod-rmse.Y2, prop=(diff/rmse.Y2) %>% round(3))





ggplot(gg.med, aes(x=BlockCol, y=BlockRow)) + 
  geom_tile(aes(fill=nu, colour=Set)) + facet_grid(.~LC) +
  scale_fill_gradient(limits=c(0, 1)) + 
  scale_colour_manual(values=c("NA", "gray70"))




