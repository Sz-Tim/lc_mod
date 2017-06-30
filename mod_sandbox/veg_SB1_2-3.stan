functions {
  row_vector tr_gjam(row_vector w) {
    int L;
    row_vector[L-1] eta;
    real w_pred;
    real D_i;
    row_vector[L] eta_full;
    real l1_sum;
    real g1_sum;
    
    L = num_elements(w);
    l1_sum = 0;
    g1_sum = 0;
    
    eta = w[1:(L-1)];
    for(i in 1:num_elements(eta)) {
      if(eta[i] < 0) {
        eta[i] = 0;
      }
      if(eta[i] < 1) {
        l1_sum = l1_sum + eta[i];
      } else {
        g1_sum = g1_sum + eta[i];
      }
    }
    w_pred = l1_sum + g1_sum;
    if(w_pred > 0.95) {
      D_i = (1 - (0.05)^(w_pred/0.95)) / w_pred;
      eta = D_i * eta;
    }
    eta_full[1:(L-1)] = eta;
    eta_full[L] = 1 - sum(eta);
    return eta_full;
  }
}
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
   Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   Y2 ~ multi_normal_cholesky(nu, L_Sigma[2]);
}
