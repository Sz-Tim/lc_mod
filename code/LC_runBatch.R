# set up environment in LC_runTest.Rmd through data dump
# mods: vector of stan model names
mods <- paste0("sc_theta/", 
               c("direct_NNGP", "direct_NNGP", 
                 "direct_NNGP_phiSc", "direct_NNGP_phiSc",
                 "direct_NNGP_wSc", "direct_NNGP_wSc"),
               ".stan")
td <- rep(10:11, 3)
run_stan <- function(mod, td, nChain=1, iter=2, warmup=1) {
  out <- stan(file=paste0("code/", mod), 
       data=read_rdump("code/all_data.Rdump"), 
       iter=iter, warmup=warmup, chains=nChain, seed=4337, init=0,
       control=list(max_treedepth=td),
       include=FALSE, pars=c("Y2_", "Y2new_", "nu","V", "uw_dp", "um"))
  return(out)
}
out.ls <- map2(mods, td, run_stan, nChain=2, iter=2000, warmup=1000)

names(out.ls) <- paste(mods, td, sep="_")



nGG <- length(get_sampler_params(out.ls[[1]], inc_warmup=FALSE)) *
  dim(get_sampler_params(out.ls[[1]], inc_warmup=FALSE)[[1]])[1]
out.gg <- map(out.ls, ~(ggs(., "n_eta") %>% 
                          arrange(Parameter, Chain, Iteration) %>% 
                          mutate(Y1=t(rbind(Y1.fit, Y1.new)) %>% 
                                   c %>% rep(each=nGG),
                                 LC=1:6 %>% rep(each=nGG) %>% rep(times=n$tot),
                                 BlockID=pop00$BlockID[n$all] %>% 
                                   rep(each=nGG*6),
                                 CellID=1:length(n$all) %>% rep(each=nGG*6),
                                 Set=c("Y1+Y2", "Y2") %>% 
                                   rep(times=c(nFit, nNew)*nGG*6)) %>%
                          mutate(BlockRow=cb.i$BlockRow[match(.$BlockID, 
                                                              cb.i$BlockID)], 
                                 BlockCol=cb.i$BlockCol[match(.$BlockID, 
                                                              cb.i$BlockID)]) %>%
                          group_by(CellID, BlockID, BlockRow, 
                                   BlockCol, LC, Set, Parameter) %>%
                          summarise(Y1=first(Y1), med=median(value), 
                                    q05=quantile(value, 0.05),
                                    q25=quantile(value, 0.25),
                                    q75=quantile(value, 0.75), 
                                    q95=quantile(value, 0.95)) %>%
                          ungroup() %>% group_by(BlockID)))

map_df(out.ls, get_elapsed_time) %>% t %>% cbind(rowSums(.))
map_dbl(out.ls, ~(sum(get_elapsed_time(.))/summary(., pars="lp__")$summary[9]))
map(out.ls, check_treedepth)
map(out.ls, check_div)

map_dfc(out.gg, ~(ungroup(.) %>% group_by(Set, LC) %>% 
                    summarise(rmse.mod=(med-Y1)^2 %>% 
                                mean %>% sqrt %>% round(3)) %>% 
                    ungroup %>% select(rmse.mod))) %>%
  add_column(LC=rep(1:6, 2), Set=rep(c("Y1+Y2", "Y2"), each=6)) 

map(out.gg, ~(ggplot(., aes(x=Y1, y=med)) + geom_point(alpha=0.5) + 
                xlim(0,1) + ylim(0,1) + facet_grid(Set~LC) + 
                geom_abline(slope=1)))

map(out.gg, ~(ggplot(., aes(x=BlockCol, y=BlockRow)) + 
                geom_tile(aes(fill=med-Y1, colour=Set)) + facet_grid(.~LC) +
                scale_fill_gradient2(limits=c(-1, 1)) + 
                scale_colour_manual(values=c("NA", "gray30"))))

map(out.ls, ~plot(., pars="beta_d"))
map(out.ls, stan_ess)
map(out.ls, stan_diag)
map(out.ls, stan_rhat)

