data {
  int N; //number of grid cells
  int L; //number of land cover classes
  row_vector[L] Y1_raw[N]; //GRANIT proportions
  row_vector[L] Y2_raw[N]; //NLCD proportions
}
// transformed data {
//   simplex[L] Y1[N];
//   simplex[L] Y2[N];
//   for(n in 1:N) {
//     Y1[n] = softmax(Y1_raw[n]');
//     Y2[n] = softmax(Y2_raw[n]');
//   }
// }
parameters {
  cholesky_factor_corr[L] L_Omega; //covariance matrix for Y1 & Y2 from nu
  vector<lower=0>[L] L_sigma;
  row_vector[L] nu_raw[N];
}
transformed parameters {
  simplex[L] nu_s[N];
  for(n in 1:N) 
    nu_s[n] = softmax(nu_raw[n]');
}
model {
  //priors
  matrix[L,L] L_Sigma;
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  
  //likelihood
   Y1_raw ~ multi_normal_cholesky(nu_raw, L_Sigma);
   Y2_raw ~ multi_normal_cholesky(nu_raw, L_Sigma);
}
