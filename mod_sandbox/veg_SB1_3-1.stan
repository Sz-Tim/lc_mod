data {
  int N;  //number of grid cells
  int L;  //number of land cover classes
  row_vector[L] Y1[N];  //GRANIT proportions
  row_vector[L-1] Y2[N];  //NLCD proportions
  vector<lower=0, upper=1>[N] p;  //pr(WP|Evg)
}
parameters {
  cholesky_factor_corr[L] L_Omega;  //covariance matrix for Y1 & Y2.ds from nu
  vector<lower=0>[L] L_sigma;  //covariance matrix for Y1 & Y2.ds from nu
  row_vector[L] nu[N];  //latent LC proportions
}
transformed parameters {
  row_vector[L] Y2_ds[N];  //de-biased and split NLCD proportions
  {
    Y2_ds[,1:4] = Y2[,1:4];
    for(n in 1:N) {
      Y2_ds[n,5] = Y2[n,5]*(1-p[n]);
      Y2_ds[n,6] = Y2[n,5]*p[n];
    }
  }
}
model {
  //priors
  matrix[L,L] L_Sigma;
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma);
}
