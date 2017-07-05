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
  int N;  //number of grid cells
  int L;  //number of land cover classes
  int nB_d[L-1];  //number of bias covariates for each LC
  int nB_p;  //number of covariates for pr(WP|Evg)
  row_vector[L] Y1[N];  //GRANIT proportions
  row_vector[L-1] Y2[N];  //NLCD proportions
  matrix[N,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[N,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[N,nB_d[2]] X_d2;  //bias covariates: Oth
  // matrix[N,nB_d[3]] X_d3;  //bias covariates: Hwd
  matrix[N,nB_d[4]] X_d4;  //bias covariates: Evg
  // matrix[N,nB_d[5]] X_d5;  //bias covariates: Mxd
}
parameters {
  cholesky_factor_corr[L] L_Omega[2]; //covariance matrix for Y1 & Y2
  vector<lower=0>[L] L_sigma[2];  //covariance matrix for Y1 & Y2
  row_vector[L] nu[N];  //latent LC proportions
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
  vector[nB_d[1]] beta_d1;  //bias betas: Dev
  vector[nB_d[2]] beta_d2;  //bias betas: Oth
  // vector[nB_d[3]] beta_d3;  //bias betas: Hwd
  vector[nB_d[4]] beta_d4;  //bias betas: Evg
  // vector[nB_d[5]] beta_d5;  //bias betas: Mxd
}
transformed parameters {
  vector<lower=0, upper=1>[N] p;  //pr(WP|Evg)
  row_vector[L-1] Y2_d[N];  //de-biased NLCD proportions
  row_vector[L] Y2_ds[N];  //de-biased and split NLCD proportions
  matrix<lower=-1, upper=1>[N,L-1] d; //bias between NLCD and nu
  simplex[L] n_eta[N];
  
  //pr(WP|Evg)
  p = inv_logit(X_p * beta_p);
  
  //estimate bias
  d[,1] = X_d1 * beta_d1;
  d[,2] = X_d2 * beta_d2;
  // d[,3] = X_d3 * beta_d3;
  d[,4] = X_d4 * beta_d4;
  // d[,5] = X_d5 * beta_d5;
  
  //bias terms set to 0
  d[,3] = to_vector(rep_array(0, N));
  d[,5] = to_vector(rep_array(0, N));
  
  for(n in 1:N) {                       // POSSIBLE TO VECTORIZE? //
    //nu to eta
    n_eta[n] = tr_gjam_inv(nu[n])';
    //de-bias NLCD proportions
    Y2_d[n] = Y2[n] + d[n];
    //split WP to [,4] and Evg to [,5]
    Y2_ds[n,4] = Y2_d[n,4]*(p[n]);
    Y2_ds[n,5] = Y2_d[n,4]*(1-p[n]);
  }
  Y2_ds[,1:3] = Y2_d[,1:3];
  Y2_ds[,6] = Y2_d[,5];
}
model {
  //priors
  matrix[L,L] L_Sigma[2];
  for(j in 1:2) {
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
    L_Omega[j] ~ lkj_corr_cholesky(4);
    L_sigma[j] ~ cauchy(0, 2.5);
  }
 
  for(l in 1:L) {
    nu[,l] ~ normal(0.5, 1);
  }
  
  beta_p ~ normal(0, 1);
  to_vector(beta_d1) ~ normal(0, 1);
  to_vector(beta_d2) ~ normal(0, 1);
  // to_vector(beta_d3) ~ normal(0, 1);
  to_vector(beta_d4) ~ normal(0, 1);
  // to_vector(beta_d5) ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma[2]);
}
