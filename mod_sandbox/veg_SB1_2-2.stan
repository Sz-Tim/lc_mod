data {
  int N; //number of grid cells
  int L; //number of land cover classes
  simplex[L] Y1[N]; //GRANIT proportions
  simplex[L] Y2[N]; //NLCD proportions
  vector<lower=0, upper=1>[L] alpha;
}
parameters {
  // cholesky_factor_corr[L] L_Omega[2]; //covariance matrix for Y1 & Y2 from nu
  // vector<lower=0>[L] L_sigma[2];
  simplex[L] nu[N];
}
// transformed parameters {
//   simplex[L] nu[N];
//   for(n in 1:N)
//     nu[n] = softmax(raw_nu[n]');
// }
model {
  //priors
  // matrix[L,L] L_Sigma[2];
  // for(j in 1:2) {
  //   L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
  //   L_Omega[j] ~ lkj_corr_cholesky(4);
  //   L_sigma[j] ~ cauchy(0, 2.5);
  // }
  for(n in 1:N) {
    nu[n] ~ dirichlet(alpha);
    Y1[n] ~ dirichlet(nu[n]);
    Y2[n] ~ dirichlet(nu[n]);
  }
  
  //likelihood
   // Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   // Y2 ~ multi_normal_cholesky(nu, L_Sigma[2]);
}
