# This script processes the output from procFull.R

library(tidyverse); library(stringr); theme_set(theme_bw())
load("data/lc_base.Rdata")
data_df <- data_df %>% arrange(CellID)

########
## Diagnostics
########
Rhat.ls <- readRDS("out/gelman_FULL_99.rds")
Gwk.ls <- readRDS("out/geweke_FULL_99.rds")

par(mfrow=c(1,1))
hist(Rhat.ls$psrf, main="Beta R-hats")
par(mfrow=c(4,6))
walk(Gwk.ls, ~{plot(.$z, ylim=c(-3,3)); abline(h=c(-1.96, 1.96))})

Gwk.df <- map(Gwk.ls, ~.$z) %>% 
  do.call(rbind, .) %>%
  as.data.frame %>%
  gather(., key=Parameter, value=z)
ggplot(Gwk.df, aes(x=z)) + geom_density() + facet_wrap(~Parameter) +
  geom_vline(xintercept=c(-1.96, 1.96), linetype=3) + 
  ggtitle("Beta Geweke diagnostics")


########
## Reshape model output
########


f.ls <- list.files("out", pattern="summary_FULL", full.names=TRUE)
out <- map(f.ls, ~read.csv(.)) %>% do.call(rbind, .) %>%
  rename(q025=X2.5., q25=X25., q50=X50., q75=X75., q975=X97.5.) %>%
  filter(grepl("n_eta", X))
cell.mx <- str_split_fixed(out$X, '\\.', 3)[,-1]
out <- out %>%
  mutate(CellID=as.integer(cell.mx[,1]),
         LC=cell.mx[,2],
         ID.LC=paste(CellID, LC, sep="_")) %>%
  group_by(CellID) %>% mutate(mnSc=Mean/sum(Mean)) %>% ungroup %>%
  arrange(LC, CellID) %>%
  mutate(LC=factor(LC, labels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd")),
         GRANIT=with(data_df, c(nhlc1_mean, nhlc2_mean, nhlc3_mean,
                                nhlc6_mean, nhlc5_mean, nhlc4_mean)),
         NLCD=with(data_df, c(nlcd1_mean, nlcd2_mean, nlcd3_mean,
                              nlcd5_mean*pWP_mean/100, 
                              nlcd5_mean*(1-pWP_mean/100), nlcd4_mean))) %>%
  full_join(., data_df %>% select(left, top, Set, pWP_mean, CellID))
sample_rows <- sample.int(nrow(out), size=ceiling(nrow(out)*.01))
out.thin <- out[sample_rows,]


# Maps: posterior
ggplot(out, aes(x=left, y=top, fill=mnSc)) + geom_tile() + 
  facet_wrap(~LC) + scale_fill_gradient(low="white", high="red")
ggplot(out, aes(x=left, y=top, fill=q50)) + geom_tile() + 
  facet_wrap(~LC) + scale_fill_gradient(low="white", high="red") 
ggplot(out, aes(x=left, y=top, fill=HPD_hi_0.99-HPD_lo_0.99)) + geom_tile() + 
  facet_wrap(~LC) + scale_fill_gradient(low="white", high="red")
ggplot(data_df, aes(x=left, y=top, fill=pWP_mean)) + geom_tile() + 
  scale_fill_gradient(low="white", high="red") 


# Maps: comparing datasets
ggplot(out, aes(x=left, y=top, fill=q50-NLCD)) + geom_tile() + 
  facet_wrap(~LC) + ggtitle("Posterior median - NLCD") +
  scale_fill_gradient2("", low="blue", mid="white", high="red", limits=c(-1,1)) 
ggplot(out, aes(x=left, y=top, fill=q50-GRANIT)) + geom_tile() + 
  facet_wrap(~LC) + ggtitle("Posterior median - GRANIT") +
  scale_fill_gradient2("", low="blue", mid="white", high="red", limits=c(-1,1))
ggplot(out, aes(x=left, y=top, fill=NLCD-GRANIT)) + geom_tile() + 
  facet_wrap(~LC) + ggtitle("NLCD - GRANIT") +
  scale_fill_gradient2("", low="blue", mid="white", high="red", limits=c(-1,1))


# Scatter plots
ggplot(out, aes(x=GRANIT, y=Mean)) + geom_point(alpha=0.01) + 
  geom_abline(linetype=3) + facet_wrap(~LC)
ggplot(out, aes(x=NLCD, y=Mean)) + geom_point(alpha=0.01) + 
  geom_abline(linetype=3) + facet_wrap(~LC)
ggplot(out.thin, aes(x=GRANIT, y=q50, ymin=HPD_lo_0.95, ymax=HPD_hi_0.95)) +
  geom_pointrange(aes(colour=NLCD), shape=95) + geom_abline(linetype=3) + 
  scale_colour_gradient2("NLCD", low="blue", mid="gray90", high="red", 
                         midpoint=0.5, limits=c(0,1)) + facet_wrap(~LC) + 
  ylim(0,1) + xlim(0,1) + labs(x="GRANIT", y="Posterior median ± 95% HPD")
ggplot(out.thin, aes(x=GRANIT, y=mnSc, colour=NLCD)) + facet_wrap(~LC) + 
  geom_point() + geom_abline(linetype=3) + 
  scale_colour_gradient2("NLCD", low="blue", mid="gray90", high="red", 
                         midpoint=0.5, limits=c(0,1)) + 
  ylim(0,1) + xlim(0,1) + labs(x="GRANIT", y="Posterior median ± 95% HPD")
ggplot(out.thin, aes(x=GRANIT, y=NLCD, colour=mnSc)) + facet_wrap(~LC) + 
  geom_point() + geom_abline(linetype=3) + 
  scale_colour_gradient2("Posterior", low="blue", mid="gray90", high="red", 
                         midpoint=0.5, limits=c(0,1)) + 
  ylim(0,1) + xlim(0,1) + labs(x="GRANIT", y="NLCD")


# Residuals
ggplot(out, aes(x=GRANIT, y=mnSc-NLCD)) + facet_wrap(~LC) +
  geom_point(alpha=0.01) + geom_hline(yintercept=0, linetype=3)
ggplot(out, aes(x=NLCD, y=mnSc-GRANIT)) + facet_wrap(~LC) +
  geom_point(alpha=0.01) + geom_hline(yintercept=0, linetype=3)
ggplot(out, aes(x=NLCD-GRANIT, y=mnSc-GRANIT)) + facet_wrap(~LC) +
  geom_point(alpha=0.01) + geom_abline(linetype=3) + geom_rug(alpha=0.01) + 
  xlim(-1,1) + ylim(-1,1) + labs(x="NLCD - GRANIT", y="Posterior - GRANIT") 
ggplot(out) + xlim(-1,1) + facet_wrap(~LC, scales="free_y") + 
  geom_vline(xintercept=0, linetype=3) + labs(x="X - GRANIT") +
  geom_density(aes(x=NLCD-GRANIT), colour="blue") +
  geom_density(aes(x=Mean-GRANIT), colour="red")
ggplot(out.thin, aes(x=GRANIT, xend=GRANIT, y=NLCD, yend=mnSc,
                     colour=abs(mnSc-GRANIT) < abs(NLCD-GRANIT))) + 
  geom_segment(arrow=arrow(length=unit(0.03, "npc")), alpha=0.2) + 
  facet_wrap(~LC) + geom_abline(linetype=3) +
  scale_colour_manual(values=c("red", "darkgreen"))
out %>% group_by(LC) %>% 
  summarise(NLCD.mn=mean(NLCD-GRANIT, na.rm=TRUE), 
            NLCD.sd=sd(NLCD-GRANIT, na.rm=TRUE), 
            post.mn=mean(Mean-GRANIT, na.rm=TRUE),
            post.sd=sd(Mean-GRANIT, na.rm=TRUE),
            mn.diff=abs(post.mn)-abs(NLCD.mn), 
            sd.diff=post.sd-abs(NLCD.sd),
            mn.pct=mn.diff/abs(NLCD.mn),
            sd.pct=sd.diff/NLCD.sd)



  
  
  
  