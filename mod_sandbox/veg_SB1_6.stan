data {
  int N;  //number of grid cells
  int L;  //number of land cover classes
  int nB_d[L-1];  //number of bias covariates for each LC
  row_vector[L] Y1[N];  //GRANIT proportions
  row_vector[L-1] Y2[N];  //NLCD proportions
  matrix[N,nB_d[1]] X_d1;  //bias covariates
  matrix[N,nB_d[2]] X_d2;
  matrix[N,nB_d[3]] X_d3;
  matrix[N,nB_d[4]] X_d4;
  matrix[N,nB_d[5]] X_d5;
  vector<lower=0, upper=1>[N] p;  //pr(WP|Evg)
}
parameters {
  cholesky_factor_corr[L] L_Omega[2]; //covariance matrix for Y1 & Y2 from nu
  vector<lower=0>[L] L_sigma[2];  //covariance matrix for Y1 & Y2 from nu
  row_vector<lower=0, upper=1>[L] nu[N];  //latent LC proportions
  vector[nB_d[1]] beta_d1;
  vector[nB_d[2]] beta_d2;
  vector[nB_d[3]] beta_d3;
  vector[nB_d[4]] beta_d4;
  vector[nB_d[5]] beta_d5;
}
transformed parameters {
  row_vector[L-1] Y2_d[N];  //de-biased NLCD proportions
  row_vector[L] Y2_ds[N];  //de-biased and split NLCD proportions
  matrix<lower=-1, upper=1>[N,L-1] d; //bias between NLCD and nu
  
  //estimate bias
  d[,1] = X_d1 * beta_d1;
  d[,2] = X_d2 * beta_d2;
  d[,3] = X_d3 * beta_d3;
  d[,4] = X_d4 * beta_d4;
  d[,5] = X_d5 * beta_d5;
  
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
  matrix[L,L] L_Sigma[2];
  for(j in 1:2) {
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
    L_Omega[j] ~ lkj_corr_cholesky(4);
    L_sigma[j] ~ cauchy(0, 2.5);
  }
  for(n in 1:N) {
    nu[n] ~ uniform(0,1);
  }
  to_vector(beta_d1) ~ normal(0, 1);
  to_vector(beta_d2) ~ normal(0, 1);
  to_vector(beta_d3) ~ normal(0, 1);
  to_vector(beta_d4) ~ normal(0, 1);
  to_vector(beta_d5) ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma[2]);
}
