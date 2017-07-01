functions {
  row_vector tr_gjam_rng(row_vector eta) {
    row_vector[6] w;
    real eta_p;
    real C_i;
    
    eta_p = sum(eta[1:5]);
    
    if(eta_p > 0.99) {
      C_i = (eta_p^(-1)) * 0.99 * (log(1-eta_p)/log(0.01));
      w = C_i*eta;
    } else {
      w = eta;
    }
    
    for(i in 1:num_elements(w)) {
      if(w[i]==0) {
        while(w[i]>=0) {
          w[i] = normal_rng(0, 0.5);
          // w[i] = -0.5;
        }
      }
    }
    return w;
  }
}
data {
  int N; //number of grid cells
  int L; //number of land cover classes
  row_vector[L] Y1[N]; //GRANIT proportions
  row_vector[L] Y2[N]; //NLCD proportions
}
transformed data {
  row_vector[L] Z1[N];
  row_vector[L] Z2[N];
  for(n in 1:N) {
    Z1[n] = tr_gjam_rng(Y1[n]);
    Z2[n] = tr_gjam_rng(Y2[n]);
  }
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
   Z1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   Z2 ~ multi_normal_cholesky(nu, L_Sigma[2]);
}
