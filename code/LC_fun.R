# Landcover Model
# Tim Szewczyk
# Accessory functions

##########
## transformations
##########

tr_gjam_inv <- function(w, a=0.99) {
  # based on Clark et al 2017
  # transformation for unconstrained w to compositional eta
  # reference class (last) is not modeled, but = 1 - sum(eta[-last])
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


tr_gjam <- function(eta, a=0.99) {
  # based on Clark et al 2017
  # transformation for compositional eta to unconstrained w
  # all eta = 0 are drawn randomly with w < 0
  require(msm)
  w <- eta
  eta.p <- sum(eta[-length(eta)])
  if(eta.p > a) {
    C.i <- (eta.p^(-1)) * a * (log(1-eta.p)/log(1-a))
    w <- C.i*eta
  }
  w[w==0] <- rtnorm(sum(w==0), mean=0, sd=0.5, upper=0)
  w
}




##########
## data munging
##########

add_blocks <- function(x, cb.i=cb.i) {
  # adds Block IDs to cells for aggregating 1-acre cells to larger blocks
  # cb.i is a reference dataframe identifying which cells belong to which blocks
  require(magrittr)
  x %>%
    mutate(BlockID=cb.i$BlockID[match(.$CellID, cb.i$CellID)]) %>%
    filter(!is.na(BlockID)) %>% 
    group_by(BlockID)
}



sampleCells <- function(nFit, nNew, nPop) {
  # samples cells for running the model at smaller extents
  totCell <- nFit + nNew
  all <- sample(1:nPop, totCell)
  fit <- all[1:nFit]
  new <- all[(nFit+1):totCell]
  return(list(all=all, fit=fit, new=new, tot=totCell))
}




##########
## placing knots for GPP
##########

place_knots <- function(mx, my, coords) {
  # from https://github.com/mbjoseph/gpp-speed-test
  xcoords <- seq(min(coords[, 2]), max(coords[, 2]), length.out = mx + 1)
  ycoords <- seq(min(coords[, 1]), max(coords[, 1]), length.out = my + 1)
  x_offset <- diff(xcoords)[1] / 2
  y_offset <- diff(ycoords)[1] / 2
  expand.grid(b.rows = ycoords[-c(my + 1)] + y_offset,
              b.cols = xcoords[-c(mx + 1)] + x_offset)
}




