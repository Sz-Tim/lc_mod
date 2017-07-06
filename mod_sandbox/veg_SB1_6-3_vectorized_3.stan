functions {
  row_vector tr_gjam_inv(row_vector w) {
    row_vector[6] eta;
    real w_p;
    real D_i;
    real ls_1;
    int gr_1;
    
    eta = w;
    ls_1 = 0;
    gr_1 = 0;
    
    for(l in 1:5) {  //replace with step functions?
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
  //counts and indices
  int N;  //number of grid cells
  int L;  //number of land cover classes
  int nB_d[L-1];  //number of bias covariates for each LC
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  row_vector[L] Y1[N];  //GRANIT proportions
  matrix[N,L-1] Y2;  //NLCD proportions
  //covariates
  matrix[N,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[N,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[N,nB_d[2]] X_d2;  //bias covariates: Oth
  matrix[N,nB_d[4]] X_d4;  //bias covariates: Evg
}

parameters {
  //covariance
  cholesky_factor_corr[L] L_Omega[2]; //covariance matrix for Y1 & Y2
  vector<lower=0>[L] L_sigma[2];  //covariance matrix for Y1 & Y2
  //landcover: latent non-constrained
  row_vector[L] nu[N];  //latent LC proportions
  //betas
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
  vector[nB_d[1]] beta_d1;  //bias betas: Dev
  vector[nB_d[2]] beta_d2;  //bias betas: Oth
  vector[nB_d[4]] beta_d4;  //bias betas: Evg
}

transformed parameters {
////////// definitions //////////
  //NLCD de-biasing and splitting
  vector[N] d_Evg;  //bias between NLCD and nu
  vector<lower=0, upper=1>[N] p;  //pr(WP|Evg)
  row_vector[L] Y2_ds[N];  //unbiased, split NLCD
  //landcover: latent compositional
  simplex[L] n_eta[N];  //gjam transformed nu
  
  
////////// operations //////////
  
  //pr(WP|Evg)
  p = inv_logit(X_p * beta_p);
  
  //store bias in Evg
  d_Evg = X_d4 * beta_d4;
  
  //add bias & split WP to [,4] and Evg to [,5]
  Y2_ds[,1] = to_array_1d(Y2[,1] + (X_d1 * beta_d1));
  Y2_ds[,2] = to_array_1d(Y2[,1] + (X_d2 * beta_d2));
  Y2_ds[,3] = to_array_1d(Y2[,1]);
  Y2_ds[,4] = to_array_1d((Y2[,4] + d_Evg) .* p);
  Y2_ds[,5] = to_array_1d((Y2[,4] + d_Evg) .* (1-p));
  Y2_ds[,6] = to_array_1d(Y2[,5]);
  
  //nu to eta
  for(n in 1:N) {
    n_eta[n] = tr_gjam_inv(nu[n])';
  }
}

model {
  matrix[L,L] L_Sigma[2];
  
  //covariance priors
  for(j in 1:2) {
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
    L_Omega[j] ~ lkj_corr_cholesky(4);
    L_sigma[j] ~ cauchy(0, 2.5);
  }
 
  //nu priors
  for(l in 1:L) {
    nu[,l] ~ normal(0.5, 1);
  }
  
  //beta priors
  beta_p ~ normal(0, 1);
  to_vector(beta_d1) ~ normal(0, 1);
  to_vector(beta_d2) ~ normal(0, 1);
  to_vector(beta_d4) ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma[2]);
}
