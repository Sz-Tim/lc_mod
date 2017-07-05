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
    
    // for(i in 1:num_elements(w)) {
    //   if(w[i]==0) {
    //     while(w[i]>=0) {
    //       w[i] = normal_rng(0, 0.25);
    //     }
    //   }
    // }
    return w;
  }
  
  row_vector tr_gjam_inv(row_vector w) {
    row_vector[6] eta;
    real w_p;
    real D_i;
    real ls_1;
    int gr_1;
    
    eta = w;
    ls_1 = 0;
    gr_1 = 0;
    
    for(l in 1:5) {
      if(eta[l] < 0) {
        eta[l] = 0;
      }
      if(eta[l] < 1) {
        ls_1 = ls_1 + eta[l];
      } else {
        gr_1 = gr_1 + 1;
      }
    }
    w_p = ls_1 + gr_1;
    
    if(w_p >= 0.99) {
      row_vector[5] tmp;
      D_i = (w_p^(-1)) * (1 - (0.01)^(w_p/0.99));
      while(sum(eta[1:5]) > 0.99) {
        tmp = D_i * eta[1:5];
        eta[1:5] = tmp;
      }
    }
    eta[6] = 1 - sum(eta[1:5]);
    return eta;
  }
}
data {
  int N; //number of grid cells
  int L; //number of land cover classes
  simplex[L] Y1[N]; //GRANIT proportions
  simplex[L] Y2[N]; //NLCD proportions
}
transformed data {
  row_vector[L] Z1[N];
  row_vector[L] Z2[N];
  for(n in 1:N) {
    Z1[n] = tr_gjam_rng(Y1[n]');
    Z2[n] = tr_gjam_rng(Y2[n]');
  }
}
parameters {
  cholesky_factor_corr[L] L_Omega[2]; //covariance matrix for Y1 & Y2 from nu
  vector<lower=0>[L] L_sigma[2];
  row_vector[L] nu[N];
}
transformed parameters {
  simplex[L] eta[N];
  for(n in 1:N) {
    eta[n] = tr_gjam_inv(nu[n])';
  }
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
