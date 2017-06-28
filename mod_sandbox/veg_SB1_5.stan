data {
  int N;  //number of grid cells
  int L;  //number of land cover classes
  int nB_d;  //number of bias covariates
  row_vector[L] Y1[N];  //GRANIT proportions
  row_vector[L-1] Y2[N];  //NLCD proportions
  matrix[N,nB_d] X_d;  //bias covariates
  vector<lower=0, upper=1>[N] p;  //pr(WP|Evg)
}
parameters {
  cholesky_factor_corr[L] L_Omega;  //covariance matrix for Y1 & Y2.ds from nu
  vector<lower=0>[L] L_sigma;  //covariance matrix for Y1 & Y2.ds from nu
  row_vector<lower=0, upper=1>[L] nu[N];  //latent LC proportions
  matrix[nB_d,L-1] beta_d;
}
transformed parameters {
  row_vector[L-1] Y2_d[N];  //de-biased NLCD proportions
  row_vector[L] Y2_ds[N];  //de-biased and split NLCD proportions
  matrix<lower=-1, upper=1>[N,L-1] d; //bias between NLCD and nu
  
  //estimate bias
  d = X_d * beta_d;
  
  for(n in 1:N) {                       // POSSIBLE TO VECTORIZE? //
    //de-bias NLCD proportions
    Y2_d[n] = Y2[n] + d[n];
    //split WP to [,5] and Evg to [,6]
    Y2_ds[n,5] = Y2_d[n,5]*(1-p[n]);
    Y2_ds[n,6] = Y2_d[n,5]*p[n];
  }
  Y2_ds[,1:4] = Y2_d[,1:4];
}
model {
  //priors
  matrix[L,L] L_Sigma;
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  for(n in 1:N) {
    nu[n] ~ uniform(0,1);
  }
  to_vector(beta_d) ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma);
}
