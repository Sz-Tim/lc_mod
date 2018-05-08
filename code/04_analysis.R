library(tidyverse); library(rstan)

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
                      Model=rep(c("3km Non", "3km CAR", "20a Non"),
                                times=c(rep(nrow(out.3km)/2,2), nrow(out.20a))),
                      mn=c(out.3km$mn, out.20a$Mean),
                      q025=c(out.3km$q025, out.20a$q025),
                      q25=c(out.3km$q25, out.20a$q25),
                      q75=c(out.3km$q75, out.20a$q75),
                      q975=c(out.3km$q975, out.20a$q975))
out.all$Model <- factor(out.all$Model, levels=levels(out.all$Model)[c(1,3,2)])
out.all$LC_full <- factor(out.all$LC, 
                          labels=c("Open", "Other", "Deciduous Forest",
                                   "White Pine Forest", "Other Evergreen Forest",
                                   "Mixed Forest"))
out.long <- data.frame(lat=c(rep(out.3km$lat,2), rep(out.20a$top,3)),
                       lon=c(rep(out.3km$lon,2), rep(out.20a$left,3)),
                       LC=factor(c(rep(out.3km$LC,2), rep(out.20a$LC,3)),
                                 labels=levels(out.3km$LC)),
                       Estimate=rep(c("3km Non", "3km CAR", "3km Y", "3km Z[s]",
                                      "20a Non", "20a Y", "20a Z[s]"),
                                    times=c(rep(nrow(out.3km)/2,4), 
                                            rep(nrow(out.20a),3))),
                       prop=c(out.3km$mn, out.3km$grnt[out.3km$Space=="non"],
                              out.3km$nlcd[out.3km$Space=="non"],
                              out.20a$Mean, out.20a$GRANIT, out.20a$NLCD),
                       Y=c(rep(out.3km$grnt,2), rep(out.20a$GRANIT,3)))
out.long$LC_full <- factor(out.long$LC, 
                           labels=c("Open", "Other", "Deciduous Forest",
                                    "White Pine Forest", "Other Evergreen Forest",
                                    "Mixed Forest"))



