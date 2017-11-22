data {
  //counts and indices
  int n;  //number of grid cells
  //landcover: observed
  vector<lower=0, upper=1>[4] Y1[n];  //Land cover proportions -- no bias
  matrix<lower=0, upper=1>[n,3] Y2;  //Land cover proportions -- with bias
  //covariates
  matrix[n,3] X_p;  //pr(WP|Evg) covariates
  matrix[n,3] X_d1;  //bias covariates: LC1
  matrix[n,2] X_d2;  //bias covariates: LC2
  matrix[n,3] X_d3;  //bias covariates: LC3
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[4] L_Omega[2];
  vector<lower=0>[4] L_sigma[2];
  //landcover: latent 
  vector[4] nu[n];
  //betas: slopes for p(WP|Evg) & bias correction
  vector[3] beta_p;  //pr(WP|Evg) betas
  vector[8] beta_d;  //bias betas
}

transformed parameters {
  //Y2 de-biased and split
  vector[4] Y2_[n];  
  
  //correct bias & split Evg (Y2[,3]) into WP (Y2_[,3]) & OthEvg (Y2_[,4])
  Y2_[,1] = to_array_1d(Y2[,1] + (X_d1 * beta_d[1:3]));
  Y2_[,2] = to_array_1d(Y2[,2] + (X_d2 * beta_d[4:5]));
  Y2_[,3] = to_array_1d((Y2[,3] + (X_d3 * beta_d[6:8])) 
      .* inv_logit(X_p * beta_p));
  Y2_[,4] = to_array_1d((Y2[,3] + (X_d3 * beta_d[6:8])) 
      .* (1-inv_logit(X_p * beta_p)));
}

model {
  //covariance priors
  for(j in 1:2) {
    L_Omega[j] ~ lkj_corr_cholesky(8);
    L_sigma[j] ~ normal(0, 1);
  }
 
  //nu priors
  for(l in 1:4) {
    nu[,l] ~ normal(0.5, 1);
  }
  
  //beta priors
  beta_p ~ normal(0, 1);
  beta_d ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, diag_pre_multiply(L_sigma[1], L_Omega[1]));
   Y2_ ~ multi_normal_cholesky(nu, diag_pre_multiply(L_sigma[2], L_Omega[2]));
}

