library(tidyverse)

out.20a <- read.csv("out/20a_out.csv")
out.20a$LC <- factor(out.20a$LC,
                     levels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
out.20a$Fit <- as.numeric(!out.20a$Set)
out.3km <- read.csv("out/3km_out.csv") 
out.3km$LC <- factor(out.3km$LC, 
                     levels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
out.all <- data.frame(lat=c(out.3km$lat, out.20a$top),
                      lon=c(out.3km$lon, out.20a$left),
                      Fit=factor(c(out.3km$Fit, out.20a$Fit),
                                 levels=0:1, labels=c("Predicted", "Fit")),
                      LC=factor(c(out.3km$LC, out.20a$LC),
                                labels=levels(out.3km$LC)),
                      grnt=c(out.3km$grnt, out.20a$GRANIT),
                      nlcd=c(out.3km$nlcd, out.20a$NLCD),
                      Model=rep(c("eta[Non]: coarse", "eta[CAR]: coarse", 
                                  "eta[Non]: fine"),
                                times=c(rep(nrow(out.3km)/2,2), nrow(out.20a))),
                      mn=c(out.3km$mn, out.20a$Mean),
                      q025=c(out.3km$q025, out.20a$q025),
                      q25=c(out.3km$q25, out.20a$q25),
                      q75=c(out.3km$q75, out.20a$q75),
                      q975=c(out.3km$q975, out.20a$q975))
out.all$LC_full <- lvls_revalue(out.all$LC, 
                                c("Open", "Other", "Deciduous Forest",
                                  "White Pine Forest", "Other Evergreen Forest",
                                  "Mixed Forest"))
out.long <- data.frame(lat=c(rep(out.3km$lat,2), rep(out.20a$top,3)),
                       lon=c(rep(out.3km$lon,2), rep(out.20a$left,3)),
                       LC=factor(c(rep(out.3km$LC,2), rep(out.20a$LC,3)),
                                 labels=levels(out.3km$LC)),
                       Estimate=rep(c("eta[Non]: coarse", 
                                      "eta[CAR]: coarse", "Y: coarse", 
                                      "Z[s]: coarse",
                                      "eta[Non]: fine", "Y: fine", 
                                      "Z[s]: fine"),
                                    times=c(rep(nrow(out.3km)/2,4), 
                                            rep(nrow(out.20a),3))),
                       prop=c(out.3km$mn, out.3km$grnt[out.3km$Space=="non"],
                              out.3km$nlcd[out.3km$Space=="non"],
                              out.20a$Mean, out.20a$GRANIT, out.20a$NLCD),
                       Y=c(rep(out.3km$grnt,2), rep(out.20a$GRANIT,3)))
out.long$Estimate <- lvls_reorder(out.long$Estimate, c(5,7,3,4,6,2,1))
out.long$LC_full <- lvls_revalue(out.long$LC, 
                                 c("Open", "Other", "Deciduous Forest",
                                   "White Pine Forest", "Other Evergreen Forest",
                                   "Mixed Forest"))


########
## FigPropMaps
########
mn.all <- ggplot(out.long, aes(x=lon, y=lat, fill=prop)) + theme_bw() + 
  geom_raster(data=filter(out.long, Estimate=="Y: fine")) +
  geom_raster(data=filter(out.long, Estimate=="Z[s]: fine")) +
  geom_raster(data=filter(out.long, Estimate=="eta[Non]: fine")) +
  geom_raster(data=filter(out.long, Estimate=="eta[Non]: coarse")) +
  geom_raster(data=filter(out.long, Estimate=="eta[CAR]: coarse")) +
  geom_raster(data=filter(out.long, Estimate=="Y: coarse")) +
  geom_raster(data=filter(out.long, Estimate=="Z[s]: coarse")) +
  scale_fill_gradient("Estimated\nProportion", 
                      low="white", high="red", limits=c(0,1)) +
  facet_grid(LC_full~Estimate, 
             labeller=labeller(LC_full=label_value, Estimate=label_parsed)) + 
  theme(axis.text=element_blank()) +
  labs(x="Easting", y="Northing")
ggsave("ms/figs/FigPropMaps.jpeg", mn.all, height=9, width=12)



########
## FigResidMaps
########
resid.all <- ggplot(droplevels(filter(out.long, 
                                      !Estimate %in% c("Y: fine", "Y: coarse"))), 
                    aes(x=lon, y=lat, fill=prop-Y)) + theme_bw() +
  geom_raster(data=filter(out.long, Estimate=="Z[s]: fine")) +
  geom_raster(data=filter(out.long, Estimate=="eta[Non]: fine")) +
  geom_raster(data=filter(out.long, Estimate=="eta[Non]: coarse")) +
  geom_raster(data=filter(out.long, Estimate=="eta[CAR]: coarse")) +
  geom_raster(data=filter(out.long, Estimate=="Z[s]: coarse")) +
  scale_fill_gradient2("Estimate - Y", 
                       low="blue", mid="white", high="red", limits=c(-1,1)) +
  facet_grid(LC_full~Estimate, 
             labeller=labeller(LC_full=label_value, Estimate=label_parsed)) + 
  theme(axis.text=element_blank()) +
  labs(x="Easting", y="Northing")
ggsave("ms/figs/FigResidMaps.jpeg", resid.all, height=9, width=9)



########
## FigVSPostZ
########
# This shows the difference between the model estimate and the simply processed
# NLCD estimate for each model, separated by LC and by cells used for fitting vs
# extrapolation/prediction. The main takeaways:
# - Fitted estimates are comparable across model for most LCs
# - Changes are typically < .25, centered around .5-.15 (=estimated bias)
# - Predictions are most similar to fitted values for CAR
# - Predictions are more similar to NLCD (vs fitted) if without SRE
# - All differ from NLCD, CAR differs more - particularly for predictions
box.nlcd <- ggplot(out.all, aes(x=Model, y=abs(mn-nlcd), fill=Fit)) + 
  geom_boxplot(outlier.alpha=0.1, outlier.size=0.25) + facet_grid(LC_full~.) +
  geom_hline(yintercept=0, colour="gray30", linetype=3) +
  scale_fill_brewer("Cells", type="qual", palette=2, direction=-1) +
  scale_x_discrete("", labels=c("CAR: coarse", "Non: coarse", "Non: fine")) +
  labs(y=expression(abs(bar(eta)-Z[s]))) + ylim(0,1) +
  guides(fill=guide_legend(reverse=TRUE)) +
  coord_flip() + theme_bw() + theme(legend.position=c(.825,.95))
ggsave("ms/figs/FigEtaZ.jpg", box.nlcd, width=4, height=9)



########
## FigCI95
########
ci.95 <- ggplot(out.all, aes(x=Model, y=q975-q025, fill=Fit)) + 
  geom_boxplot(outlier.alpha=0.1, outlier.size=0.25) + facet_grid(LC_full~.) +
  geom_hline(yintercept=0, colour="gray30", linetype=3) +
  scale_fill_brewer("Cells", type="qual", palette=2, direction=-1) +
  scale_x_discrete("", labels=parse(text=levels(out.all$Model))) +
  labs(x="", y="95% CI") + guides(fill=guide_legend(reverse=TRUE)) +
  coord_flip() + theme_bw() + theme(legend.position=c(.825,.95))
ggsave("ms/figs/FigCI95.jpg", ci.95, width=4, height=9)



########
## RMSE -- FULL RUN
########
# For cells used for fitting, CAR puts eta estimates very near GRANIT, while
# the models without SRE put etas between GRANIT & NLCD. Eta estimates for all
# models are nearer to GRANIT than NLCD is.
RMSE.3km <- out.3km %>% filter(Fit==1) %>%
  group_by(Space, LC) %>% 
  summarise(RMSE_nlcd=sqrt(mean( (nlcd - grnt)^2 )),
            RMSE_mod=sqrt(mean( (mn - grnt)^2 ))) %>%
  mutate(pct_diff=(RMSE_mod - RMSE_nlcd)/RMSE_nlcd*100) %>%
  rename(Model=Space) %>% ungroup
RMSE.3km$Model <- lvls_revalue(RMSE.3km$Model, 
                               c("eta[CAR]: coarse", "eta[Non]: coarse"))
RMSE.20a <- out.20a %>% filter(Fit==1) %>%
  group_by(LC) %>%
  summarise(RMSE_nlcd=sqrt(mean( (NLCD - GRANIT)^2 )),
            RMSE_mod=sqrt(mean( (Mean - GRANIT)^2 ))) %>%
  mutate(pct_diff=(RMSE_mod - RMSE_nlcd)/RMSE_nlcd*100,
         Model="eta[Non]: fine")

write_csv(rbind(RMSE.3km, RMSE.20a), "ms/figs/SuppFullRMSE.csv")


