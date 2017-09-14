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



sampleCells <- function(nFit, nNew, nPop, partition=FALSE) {
  # samples cells for running the model at smaller extents
  # partition controls whether the cells are randomly assigned or partitioned
  totCell <- nFit + nNew
  if(partition) {
    all <- 1:totCell
  } else {
    all <- sample(1:nPop, totCell)
  }
    fit <- all[1:nFit]
    new <- all[(nFit+1):totCell]
  return(list(all=all, fit=fit, new=new, tot=totCell))
}




##########
## placing knots for GPP
##########

place_knots <- function(mx, my, coords) {
  # from https://github.com/mbjoseph/gpp-speed-test
  xcoords <- seq(min(coords[, 1]), max(coords[, 1]), length.out = mx + 1)
  ycoords <- seq(min(coords[, 2]), max(coords[, 2]), length.out = my + 1)
  x_offset <- diff(xcoords)[1] / 2
  y_offset <- diff(ycoords)[1] / 2
  expand.grid(b.cols = xcoords[-c(mx + 1)] + x_offset,
              b.rows = ycoords[-c(my + 1)] + y_offset)
}




##########
## NNGP setup: Lu Zhang example
##########
# 
# i_index <- function(i, s, M) {
#   require(fields)
#   if(M >= (i - 1)) {im <- 1:(i - 1)}
#   else 	{
#     dist <- rdist(s[c(1,i),], s[c(1:(i-1)), ])[-1,]
#     im <- sort(order(dist)[1:M])
#   }
#   return(im)
# }
# 
# 
# 
# i_dist <- function(i, neighbor_index, s) {
#   dist(s[c(i, neighbor_index[[i - 1]]), ])
# }
# 
# 
# 
# get_index_dist <- function(s, M) {
#   
#   n = nrow(s)
#   M = min(M, n-1)
#   # get index of neighborhood
#   neighbor_index <- sapply(2:n, i_index, s, M)
#   # get distance matrix for each i
#   neighbor_dist <- sapply(2:n, i_dist, neighbor_index, s)
#   
#   return(list(i = neighbor_index, d = neighbor_dist))
# }
# 
# 
# 
# get_neardistM <- function (ind, ind_distM_d) {
#   if (ind < M ){l = ind } else {l = M}; 
#   M_i <- rep(0, M * (M - 1) / 2);
#   if (l == 1) {}
#   else{
#     M_i[1: (l * (l - 1) / 2)] <- 
#       c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
#   }
#   return(M_i)
# }
# 
# 
# 
# get_neardist <- function (ind, ind_distM_d) {
#   if (ind < M ){l = ind } else {l = M}; 
#   D_i <- rep(0, M);
#   D_i[1:l]<-c(ind_distM_d[[ind]])[1:l]
#   return( D_i)
# }
# 
# 
# 
# get_nearind <- function (ind, ind_distM_i) {
#   if (ind < M ){l = ind } else {l = M}; 
#   D_i <- rep(0, M);
#   D_i[1:l]<-c(ind_distM_i[[ind]])[1:l]
#   return( D_i)
# }



##########
## NNGP setup: gridded pixels
##########

i_index <- function(i, s, M_r) {
  ### generates list of nearest neighbor indices for pixel i
  # i: row index within coordinates s
  # s: coordinates as pixel column[,1] & row[,2]
  # M_r: neighborhood size as ± M_r rows & cols
  require(dplyr)
  x_nn <- seq.int(max(s[i,1]-M_r, min(s[,1])), min(s[i,1]+M_r, max(s[,1])))
  y_nn <- seq.int(max(s[i,2]-M_r, min(s[,2])), min(s[i,2]+M_r, max(s[,2])))
  nn <- filter(expand.grid(nn.col=x_nn, nn.row=y_nn), 
               !(nn.col==s[i,1] & nn.row==s[i,2]))
  im <- which(paste(s[,1], s[,2]) %in% paste(nn[,1], nn[,2]))
  return(im)
}

i_dist <- function(i, nn_ind, s) {
  ### calculates nearest neighbor distances for pixel i
  dist(s[c(i, nn_ind[[i]]), ])
}

get_index_dist <- function(s, M_r) {
  ### generates nn indices and distance matrix for each pixel in s
  # s: coordinates as pixel column[,1] & row[,2]
  # M_r: neighborhood size as ± M_r rows & cols
  require(purrr)
  # generate indices, order, & distances
  nn_index <- sapply(1:nrow(s), i_index, s, M_r)
  nn_order <- order(purrr::map_int(nn_index, length))
  nn_dist <- sapply(1:nrow(s), i_dist, nn_index, s)[nn_order]
  # reorder
  nn_index <- nn_index[nn_order]
  s_nn <- s[nn_order,]
  # create index map for nn vs all else
  n_M <- purrr::map_int(nn_index, length)
  nn_YX <- cbind(XY_id=nn_order,
                 nn_grp=as.numeric(as.factor(n_M)))
  # reference for indexing nn operations by neighborhood size
  nn_M <- unique(n_M)
  n_M_ref <- matrix(nrow=length(nn_M), ncol=3)
  colnames(n_M_ref) <- c("nn_M", "i_start", "i_end")
  for(r in 1:length(nn_M)) {
    n_M_ref[r,] <- c(nn_M[r], range(which(n_M == nn_M[r])))
  }
  
  return(list(i=nn_index, d=nn_dist, 
              s_nn=s_nn, n_M_ref=n_M_ref, 
              nn_YX=nn_YX))
}

get_neardistM <- function (i, id_d_M, M_r) {
  ### generates vector of distance matrix for pixel i
  # i: row index within coordinates s
  # id_d_M: output from get_index_dist with i=nn_index & d=nn_distance
  max_nn <- (1 + 2*M_r)^2 - 1
  n_nn <- length(id_d_M$i[[i]])
  M_i <- rep(0, max_nn*(max_nn-1)/2)
  M_i[1:(n_nn*(n_nn-1)/2)] <- c(id_d_M$d[[i]])[(n_nn+1):(n_nn*(n_nn+1)/2)]
  return(M_i)
}

get_neardist <- function (i, id_d_M, M_r) {
  ### generates vector of distances from pixel i
  # i: row index within coordinates s
  # id_d_M: output from get_index_dist with i=nn_index & d=nn_distance
  # M_r: neighborhood size as ± M_r rows & cols
  max_nn <- (1 + 2*M_r)^2 - 1
  n_nn <- length(id_d_M$i[[i]])
  D_i <- rep(0, max_nn)
  D_i[1:n_nn] <- c(id_d_M$d[[i]])[1:n_nn]
  return( D_i)
}

get_nearind <- function (i, id_d_M, M_r) {
  ### generates vector of indices to pair with get_neardist
  # i: row index within coordinates s
  # id_d_M: output from get_index_dist with i=nn_index & d=nn_distance
  # M_r: neighborhood size as ± M_r rows & cols
  max_nn <- (1 + 2*M_r)^2 - 1
  n_nn <- length(id_d_M$i[[i]])
  D_i <- rep(0, max_nn)
  D_i[1:n_nn] <- c(id_d_M$i[[i]])[1:n_nn]
  return( D_i)
}


