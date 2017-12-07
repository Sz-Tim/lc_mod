# This script performs stratified random sampling to create the variable 
# selection grid. The full grid takes too long to run, so variable selection is
# done on a subset of the cells, with each variable represented proportionally
# in the subset as it is found on the landscape. As the variable selection is 
# only possible for cells with both GRANIT & NLCD, it is restricted to cells 
# with data_df$Set == FALSE.

library(tidyverse); library(purrr); library(stringr)
load("data/lc_base.Rdata")

# extract covariates
NH_df <- data_df %>% filter(!Set) %>%
  arrange(left, top) %>% 
  mutate(id=row_number(), 
         x=as.numeric(as.factor(str_pad(left, 3, "left", "0"))),
         y=as.numeric(as.factor(str_pad(top, 3, "left", "0"))))
X <- NH_df %>% 
  select(pWP_mean, bio1_mean, bio7_mean, bio12_mean,
         el_mean, rugg_mean, rdLen_, pop, hms)

# calculate quartiles
X.q <- X %>% mutate_all(ntile, n=5)

n <- floor(nrow(X)*0.15)

samp.id <- sample(1:nrow(X), n, replace=FALSE)
x.sum <- sapply(X.q[samp.id,], table)
x.sum

n.q <- colSums(x.sum)[1]/n_distinct(X.q$el_mean)

i <- 0
while( sum(abs(x.sum - n.q) > 0.01*n.q) != 0 ) {
  samp.id <- sample(1:nrow(X), n, replace=FALSE)
  x.sum <- sapply(X.q[samp.id,], table)
  i <- i+1
  if((i %% 100) == 0) cat("Still trying...", i, "attempts\n")
}
x.sum

X.q.samp <- X.q[samp.id,]
names(X.q.samp) <- paste0(names(X.q.samp), "_q")
X_df <- cbind(NH_df[samp.id,], X.q.samp)

write_csv(data.frame(id=samp.id), "data/stratified_sample_15pct_20a_rowID.csv")


theme_set(theme_bw())
ggplot(X_df, aes(x=x, y=y, colour=factor(pWP_mean_q))) + 
  geom_point(size=0.1) + scale_colour_brewer(palette=8)
ggplot(X_df, aes(x=x, y=y, fill=factor(pWP_mean_q))) + 
  geom_tile() + scale_fill_brewer(palette=8)
