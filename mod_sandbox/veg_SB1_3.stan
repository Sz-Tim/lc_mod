data {
  int N; //number of grid cells
  int L; //number of land cover classes
  row_vector[L] Y1[N]; //GRANIT proportions
  row_vector[L] Y2[N]; //NLCD proportions
}
parameters {
  cholesky_factor_corr[L] L_Omega[2]; //covariance matrix for Y1 & Y2 from nu
  vector<lower=0>[L] L_sigma[2];
  row_vector[L] nu[N];
}
model {
  //priors
  matrix[L,L] L_Sigma[2];
  for(j in 1:2) {
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
    L_Omega[j] ~ lkj_corr_cholesky(4);
    L_sigma[j] ~ cauchy(0, 2.5);
  }
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma[j]);
   Y2 ~ multi_normal_cholesky(nu, L_Sigma[j]);
}
