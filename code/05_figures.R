# This script produces the figures for the manuscript. It reads in a csv of the
# model output at each resolution, combines them, and reshapes as needed. The
# figures produced here are:
# - FigPropMaps: Map of land cover proportions for each LC * estimate
# - FigResidMaps: Map of residuals (x - GRANIT) for NLCD and models
# - FigEtaZ: Density plot of abs(model - NLCD) for fitted vs predicted cells
# - Fig95CI: Boxplots of 95% credible intervals for model output
# - RMSE TABLE: RMSE (x : GRANIT) for fitted cells in full run
# - SuppScatter: Scatter plot of x vs GRANIT


library(tidyverse); theme_set(theme_bw())
out.20a <- read.csv("out/20a_out.csv")
out.20a$LC <- factor(out.20a$LC,
                     levels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
out.20a$Fit <- as.numeric(!out.20a$Set)
out.9km <- read.csv("out/9km_out.csv") 
out.9km$LC <- factor(out.9km$LC, 
                     levels=c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd"))
out.all <- data.frame(lat=c(out.9km$lat, out.20a$top),
                      lon=c(out.9km$lon, out.20a$left),
                      Fit=factor(c(out.9km$Fit, out.20a$Fit),
                                 levels=0:1, labels=c("Predict", "Fit")),
                      LC=factor(c(out.9km$LC, out.20a$LC),
                                labels=levels(out.9km$LC)),
                      grnt=c(out.9km$grnt, out.20a$GRANIT),
                      nlcd=c(out.9km$nlcd, out.20a$NLCD),
                      Model=rep(c("eta[Non]: 9km^2", "eta[CAR]: 9km^2", 
                                  "eta[Non]: 20ac"),
                                times=c(rep(nrow(out.9km)/2,2), nrow(out.20a))),
                      Res=rep(c("coarse", "fine"), 
                              times=c(nrow(out.9km), nrow(out.20a))),
                      mn=c(out.9km$mn, out.20a$Mean),
                      q025=c(out.9km$q025, out.20a$q025),
                      q25=c(out.9km$q25, out.20a$q25),
                      q75=c(out.9km$q75, out.20a$q75),
                      q975=c(out.9km$q975, out.20a$q975))
out.all$LC_full <- lvls_revalue(out.all$LC, 
                                c("Open Canopy", "Other", "Decid. Forest",
                                  "White Pine Forest", "Oth. Evg. Forest",
                                  "Mixed Forest"))
out.all$LC_full <- lvls_reorder(out.all$LC_full, c(1,3,4,5,6,2))
out.all$Model <- lvls_reorder(out.all$Model, c(1,3,2))
out.all$Mod_parse <- lvls_revalue(out.all$Model, 
                                  c("CAR~~'9 km'^2",
                                    "Non~~'9 km'^2",
                                    "Non~~'20 ac'"))
out.long <- data.frame(lat=c(rep(out.9km$lat,2), rep(out.20a$top,3)),
                       lon=c(rep(out.9km$lon,2), rep(out.20a$left,3)),
                       LC=factor(c(rep(out.9km$LC,2), rep(out.20a$LC,3)),
                                 labels=levels(out.9km$LC)),
                       Estimate=rep(c("bar(eta)[Non]~~'9 km'^2", 
                                      "bar(eta)[CAR]~~'9 km'^2", 
                                      "Y~~'9 km'^2", 
                                      "Z^s~~'9 km'^2",
                                      "bar(eta)[Non]~~'20 ac'", 
                                      "Y~~'20 ac'", 
                                      "Z^s~~'20 ac'"),
                                    times=c(rep(nrow(out.9km)/2,4), 
                                            rep(nrow(out.20a),3))),
                       prop=c(out.9km$mn, out.9km$grnt[out.9km$Space=="non"],
                              out.9km$nlcd[out.9km$Space=="non"],
                              out.20a$Mean, out.20a$GRANIT, out.20a$NLCD),
                       Y=c(rep(out.9km$grnt,2), rep(out.20a$GRANIT,3)))
out.long$Estimate <- lvls_reorder(out.long$Estimate, c(6,4,2,7,5,3,1))
out.long$LC_full <- lvls_revalue(out.long$LC, 
                                 c("Open Canopy", "Other", "Decid. Forest",
                                   "White Pine Forest", "Oth. Evg. Forest",
                                   "Mixed Forest"))
out.long$LC_full <- lvls_reorder(out.long$LC_full, c(1,3,4,5,6,2))

# Outline of Fit/Predict regions at 20 acre resolution
coord.20 <- filter(out.20a, LC=="Opn") %>% 
  select(left, top, Fit) %>% 
  raster::rasterFromXYZ() %>% 
  raster::rasterToPolygons(dissolve=T) %>%
  fortify %>% select(long, lat, group) %>%
  rename(lon=long)
coord.20.z <- coord.20.non <- coord.20
coord.20.y <- filter(coord.20, group %in% c("2.1", "2.2", "2.3"))
coord.20.z$Estimate <- "Z^s~~'20 ac'"
coord.20.non$Estimate <- "bar(eta)[Non]~~'20 ac'"
coord.20.y$Estimate <- "Y~~'20 ac'"

# Outline of Fit/Predict regions at 9 km2 resolution
coord.9 <- filter(out.9km, LC=="Opn" & Space=="non") %>% 
  select(lon, lat, Fit) %>% 
  raster::rasterFromXYZ() %>% 
  raster::rasterToPolygons(dissolve=T) %>%
  fortify %>% select(long, lat, group) %>%
  rename(lon=long)
coord.9.z <- coord.9.non <- coord.9.car <- coord.9
coord.9.y <- filter(coord.9, group == "2.1")
coord.9.z$Estimate <- "Z^s~~'9 km'^2"
coord.9.non$Estimate <- "bar(eta)[Non]~~'9 km'^2"
coord.9.car$Estimate <- "bar(eta)[CAR]~~'9 km'^2"
coord.9.y$Estimate <- "Y~~'9 km'^2"

max.coord <- rbind(coord.20.y, coord.20.z, coord.20.non, 
                   coord.9.y, coord.9.z, coord.9.non, coord.9.car)
max.coord$Estimate <- factor(max.coord$Estimate, levels=levels(out.long$Estimate))


gg_fonts <- theme(axis.title=element_text(size=12),
                  axis.text=element_text(size=10, colour="gray20"),
                  legend.title=element_text(size=10),
                  legend.text=element_text(size=9),
                  panel.grid.minor=element_blank(),
                  strip.text=element_text(size=10),
                  plot.title=element_text(size=12))
fit.col <- c(Fit="#40004b", Predict="#1b7837")
fit.fill <- c(Fit="#762a83", Predict="#5aae61")
coord.col <- c(fit.col[2], fit.fill[1])[as.numeric(substr(levels(max.coord$group),1,1))]
names(coord.col) <- levels(max.coord$group)



########
## FigPropMaps
########
mn.all <- ggplot(out.long, aes(x=lon, y=lat, fill=prop)) + theme_bw() + 
  geom_raster(data=filter(out.long, Estimate=="Z^s~~'20 ac'")) +
  geom_raster(data=filter(out.long, Estimate=="Y~~'20 ac'")) +
  geom_raster(data=filter(out.long, Estimate=="bar(eta)[Non]~~'20 ac'")) +
  geom_raster(data=filter(out.long, Estimate=="Z^s~~'9 km'^2")) +
  geom_raster(data=filter(out.long, Estimate=="Y~~'9 km'^2")) +
  geom_raster(data=filter(out.long, Estimate=="bar(eta)[Non]~~'9 km'^2")) +
  geom_raster(data=filter(out.long, Estimate=="bar(eta)[CAR]~~'9 km'^2")) +
  # geom_path(data=max.coord, inherit.aes=F, aes(x=lon, y=lat, group=group), 
  #           colour="gray20", size=0.25) + 
  geom_path(data=max.coord, inherit.aes=F, size=0.4,
            aes(x=lon, y=lat, group=group, colour=group)) + 
  scale_colour_manual(values=coord.col, guide="none") +
  scale_fill_gradient2("Proportion Cover", na.value=NA, midpoint=0.5,
                       low="#fff7fb", mid="#74a9cf", high="#023858", 
                      limits=c(0,1), breaks=c(0, 0.5, 1), guide="colourbar") +
  facet_grid(Estimate~LC_full, 
             labeller=labeller(LC_full=label_wrap_gen(width=11), 
                               Estimate=label_parsed)) + 
  gg_fonts + theme(axis.text.x=element_text(angle=300, hjust=0, vjust=1),
                   axis.text.y=element_text(angle=30, hjust=1, vjust=0.5),
                   legend.position="bottom", 
                   legend.key.height=unit(0.1, "in"),
                   panel.background=element_rect(fill="gray20"), 
                   panel.grid.major=element_blank(),
                   panel.spacing.y=unit(c(.1, .1, .4, .1, .1, .1), "lines"),
                   panel.spacing.x=unit(.1, "lines"),
                   legend.margin=margin(-.3,0,0,3.5, unit="in")) +
  guides(fill=guide_colourbar(title.position="top", title.hjust=0.5)) +
  scale_x_continuous("Easting", breaks=c(700000, 750000, 800000, 850000),
                     labels=c(expression(70^scriptscriptstyle("0000")), 
                              expression(75^scriptscriptstyle("0000")),
                              expression(80^scriptscriptstyle("0000")), 
                              expression(85^scriptscriptstyle("0000")))) +
  scale_y_continuous("Northing", breaks=c(4750000, 4800000, 4850000),
                     labels=c(expression(475^scriptscriptstyle("0000")), 
                              expression(480^scriptscriptstyle("0000")),
                              expression(485^scriptscriptstyle("0000"))))
ggsave("ms/figs/FigPropMaps.jpeg", mn.all, height=7, width=6, 
       units="in", dpi=600)



########
## FigResidMaps
########
resid.all <- ggplot(droplevels(filter(out.long, 
                                      !Estimate %in% c("Y~~'20 ac'", "Y~~'9 km'^2"))), 
                    aes(x=lon, y=lat, fill=prop-Y)) + theme_bw() +
  geom_raster(data=filter(out.long, Estimate=="Z^s~~'20 ac'")) +
  geom_raster(data=filter(out.long, Estimate=="bar(eta)[Non]~~'20 ac'")) +
  geom_raster(data=filter(out.long, Estimate=="Z^s~~'9 km'^2")) +
  geom_raster(data=filter(out.long, Estimate=="bar(eta)[Non]~~'9 km'^2")) +
  geom_raster(data=filter(out.long, Estimate=="bar(eta)[CAR]~~'9 km'^2")) +
  scale_fill_gradient2("Estimate - Y", low="#313695", mid="white", na.value=NA,
                       high="#a50026", limits=c(-1,1), breaks=c(-1, 0, 1), 
                       guide="colourbar") +
  facet_grid(Estimate~LC_full, 
             labeller=labeller(LC_full=label_wrap_gen(width=11), 
                               Estimate=label_parsed)) + 
  gg_fonts + theme(axis.text.x=element_text(angle=300, hjust=0, vjust=1),
                   axis.text.y=element_text(angle=30, hjust=1, vjust=0.5),
                   legend.position="bottom", 
                   legend.key.height=unit(0.1, "in"),
                   panel.background=element_rect(fill="gray20"), 
                   panel.grid.major=element_blank(),
                   panel.spacing.y=unit(c(.1, .4, .1, .1), "lines"),
                   panel.spacing.x=unit(.1, "lines"),
                   legend.margin=margin(-.3,0,0,3.5, unit="in")) +
  guides(fill=guide_colourbar(title.position="top", title.hjust=0.5)) +
  scale_x_continuous("Easting", breaks=c(700000, 750000, 800000, 850000),
                     labels=c(expression(70^scriptscriptstyle("0000")), 
                              expression(75^scriptscriptstyle("0000")),
                              expression(80^scriptscriptstyle("0000")), 
                              expression(85^scriptscriptstyle("0000")))) +
  scale_y_continuous("Northing", breaks=c(4750000, 4800000, 4850000),
                     labels=c(expression(475^scriptscriptstyle("0000")), 
                              expression(480^scriptscriptstyle("0000")),
                              expression(485^scriptscriptstyle("0000"))))
ggsave("ms/figs/FigResidMaps.jpeg", resid.all, height=5.5, width=6, 
       units="in", dpi=600)



########
## FigEtaZ
########
# This shows the difference between the model estimate and the simply processed
# NLCD estimate for each model, separated by LC and by cells used for fitting vs
# extrapolation/prediction. The main takeaways:
# - Fitted estimates are comparable across model for most LCs
# - Changes are typically < .25, centered around .5-.15 (=estimated bias)
# - Predictions are most similar to fitted values for CAR
# - Predictions are more similar to NLCD (vs fitted) if without SRE
# - All differ from NLCD, CAR differs more - particularly for predictions
dens.nlcd <- ggplot(out.all, aes(x=abs(mn-nlcd), 
                                 colour=fct_rev(Fit), fill=fct_rev(Fit))) + 
  geom_density(alpha=0.4) +
  geom_hline(yintercept=0, colour="gray30") +
  facet_grid(Mod_parse~LC_full, scales="free_y",
             labeller=labeller(LC_full=label_wrap_gen(width=11), 
                               Mod_parse=label_parsed)) +
  scale_colour_manual("", values=fit.col) +
  scale_fill_manual("", values=fit.fill) +
  scale_x_continuous(name=expression(abs(bar(eta)-Z^s)), limits=c(0, NA),
                     breaks=c(0, 0.5), labels=c("0", "0.5")) +
  labs(y="Density") +
  theme_bw() + theme(legend.position=c(.935,.95),
                     legend.background=element_blank(),
                     legend.key.size=unit(.15, "in"),
                     panel.grid=element_blank()) + gg_fonts
ggsave("ms/figs/FigEtaZ.jpg", dens.nlcd, width=6, height=3.5, units="in", dpi=600)



########
## FigCI95
########
ci.95 <- ggplot(out.all, aes(x=fct_rev(Model), y=q975-q025, 
                             colour=fct_rev(Fit), fill=fct_rev(Fit))) + 
  facet_grid(.~LC_full, labeller=label_wrap_gen(width=11)) +
  geom_boxplot(outlier.alpha=0.1, outlier.size=0.1, alpha=0.5) + 
  geom_hline(yintercept=0, colour="gray30", size=0.25) +
  scale_x_discrete("", labels=parse(text=c("Non~~'20 ac'",
  "Non~~'9 km'^2", "CAR~~'9 km'^2"))) +
  labs(x="", y=expression(eta~~"(95% CI)")) + 
  scale_colour_manual("", values=fit.col) +
  scale_fill_manual("", values=fit.fill) +
  theme_bw() + theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.4),
                     legend.position=c(.905,.85),
                     legend.background=element_blank(),
                     legend.key.size=unit(.15, "in"),
                     panel.grid=element_blank()) + gg_fonts
ggsave("ms/figs/FigCI95.jpg", ci.95, width=6, height=2.5, units="in", dpi=600)



########
## RMSE -- FULL RUN
########
# For cells used for fitting, CAR puts eta estimates very near GRANIT, while
# the models without SRE put etas between GRANIT & NLCD. Eta estimates for all
# models are nearer to GRANIT than NLCD is.
RMSE.9km <- out.9km %>% filter(Fit==1) %>%
  group_by(Space, LC) %>% 
  summarise(RMSE_nlcd=sqrt(mean( (nlcd - grnt)^2 )),
            RMSE_mod=sqrt(mean( (mn - grnt)^2 ))) %>%
  mutate(pct_diff=(RMSE_mod - RMSE_nlcd)/RMSE_nlcd*100) %>%
  rename(Model=Space) %>% ungroup
RMSE.9km$Model <- lvls_revalue(RMSE.9km$Model, 
                               c("eta[CAR]: coarse", "eta[Non]: coarse"))
RMSE.20a <- out.20a %>% filter(Fit==1) %>%
  group_by(LC) %>%
  summarise(RMSE_nlcd=sqrt(mean( (NLCD - GRANIT)^2 )),
            RMSE_mod=sqrt(mean( (Mean - GRANIT)^2 ))) %>%
  mutate(pct_diff=(RMSE_mod - RMSE_nlcd)/RMSE_nlcd*100,
         Model="eta[Non]: fine")

write_csv(rbind(RMSE.9km, RMSE.20a), "ms/figs/SuppFullRMSE.csv")



########
## Y vs Z: Scatter plots
########
out.all$GrdRes <- lvls_revalue(out.all$Res, c("'9 km'^2", "'20 ac'"))
scat.YZ <- ggplot(out.all, aes(x=nlcd, y=grnt)) + 
  geom_point(data=filter(out.all, Model == "eta[Non]: 9km^2"), alpha=0.1) + 
  geom_point(data=filter(out.all, Model == "eta[Non]: 20ac"), alpha=0.01) +
  geom_abline(colour="red", linetype=2) +
  facet_grid(GrdRes~LC_full, labeller=labeller(LC_full=label_wrap_gen(width=11), 
                                              GrdRes=label_parsed)) + 
  scale_x_continuous(expression(Z^s), breaks=c(0, 0.5, 1), 
                     labels=c("0", "0.5", "1")) + 
  scale_y_continuous("Y", breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1")) + 
  gg_fonts + theme_bw() +
  theme(panel.spacing=unit(.1, "lines"),
        panel.grid=element_blank())
ggsave("ms/supp/SuppScatter.jpeg", scat.YZ, height=2.75, width=6)

scat.Y_all <- ggplot(droplevels(filter(out.long, 
                         !Estimate %in% c("Y~~'20 ac'", "Y~~'9 km'^2"))), 
       aes(x=prop, y=Y)) + theme_bw() +
  geom_point(data=filter(out.long, Estimate=="Z^s~~'20 ac'"), 
             alpha=0.01) +
  geom_point(data=filter(out.long, Estimate=="bar(eta)[Non]~~'20 ac'"), 
             alpha=0.01) +
  geom_point(data=filter(out.long, Estimate=="Z^s~~'9 km'^2"), 
             alpha=0.1) +
  geom_point(data=filter(out.long, Estimate=="bar(eta)[Non]~~'9 km'^2"), 
             alpha=0.1) +
  geom_point(data=filter(out.long, Estimate=="bar(eta)[CAR]~~'9 km'^2"), 
             alpha=0.1) +
  geom_abline(colour="red", linetype=2) +
  facet_grid(Estimate~LC_full, 
             labeller=labeller(LC_full=label_wrap_gen(width=11), 
                               Estimate=label_parsed)) + 
  scale_x_continuous(expression(Z^s~~'or'~~bar(eta)~~'proportion'), 
                     breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1")) + 
  scale_y_continuous("Y", breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1")) + 
  gg_fonts + theme(panel.spacing.x=unit(.2, "lines"),
                   panel.spacing.y=unit(c(.2, .6, .2, .2), "lines"),
                   panel.grid=element_blank())
ggsave("ms/supp/SuppScatter_all.jpeg", scat.Y_all, height=5, width=6)

dens.YZ <- ggplot(out.all, aes(x=nlcd-grnt)) + 
  geom_vline(xintercept=0, linetype=3) +
  geom_density(data=filter(out.all, Model == "eta[Non]: 9km^2"), 
               colour="black") +
  geom_density(data=filter(out.all, Model == "eta[Non]: 20ac"), 
               colour="gray40") + 
  facet_wrap(~LC_full) + xlim(-1,1) + 
  labs(x=expression(Z[s]-Y))
dens.YZ


