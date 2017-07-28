functions {
  row_vector tr_gjam_inv(row_vector w) {
    row_vector[6] eta;
    real w_p;
    real D_i;
    real ls_1;
    int gr_1;
    
    eta[1:5] = w;
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
  //counts and indices
  int n1;  //number of cells for GRANIT
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells for NLCD + covariates
  int L;  //number of land cover classes
  int nB_d[L-2];  //number of bias covariates for each LC
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  row_vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n3,L-2] Y2;  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[n3,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[n3,nB_d[2]] X_d2;  //bias covariates: Oth
  matrix[n3,nB_d[3]] X_d3;  //bias covariates: Hwd
  matrix[n3,nB_d[4]] X_d4;  //bias covariates: Evg
}

parameters {
  //covariance
  cholesky_factor_corr[L-1] L_Omega[2]; //covariance matrix for Y1 & Y2
  vector<lower=0, upper=pi()/2>[L-1] L_sigma_unif[2];  //covariance
  //landcover: latent non-constrained
  row_vector<lower=-1, upper=2>[L-1] nu[n3];  //latent LC proportions
  //betas
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
  vector[nB_d[1]] beta_d1;  //bias betas: Dev
  vector[nB_d[2]] beta_d2;  //bias betas: Oth
  vector[nB_d[3]] beta_d3;  //bias betas: Hwd
  vector[nB_d[4]] beta_d4;  //bias betas: Evg
}

transformed parameters {
  //covariance
  vector<lower=0>[L-1] L_sigma[2];  //covariance matrix for Y1 & Y2
  matrix[L-1,L-1] L_Sigma[2];
  //NLCD de-biasing and splitting
  row_vector[L-1] Y2_ds[n3];  //unbiased, split NLCD
  //landcover: latent compositional
  simplex[L] n_eta[n3];  //gjam transformed nu
  
  
  // L_sigma ~ cauchy(0, 2.5) --- reparameterized for speed
  for(j in 1:2) {
    L_sigma[j] = 2.5 * tan(L_sigma_unif[j]);
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
  }
  
  
  //correct bias & split WP to [,4] and Evg to [,5]
  ////fit betas using cells with Y1 & Y2
  Y2_ds[1:n1,1] = to_array_1d(Y2[1:n1,1] + (X_d1[1:n1,] * beta_d1));
  Y2_ds[1:n1,2] = to_array_1d(Y2[1:n1,2] + (X_d2[1:n1,] * beta_d2));
  Y2_ds[1:n1,3] = to_array_1d(Y2[1:n1,3] + (X_d3[1:n1,] * beta_d3));
  Y2_ds[1:n1,4] = to_array_1d((Y2[1:n1,4] + (X_d4[1:n1,] * beta_d4)) 
      .* inv_logit(X_p[1:n1,] * beta_p));
  Y2_ds[1:n1,5] = to_array_1d((Y2[1:n1,4] + (X_d4[1:n1,] * beta_d4)) 
      .* (1-inv_logit(X_p[1:n1,] * beta_p)));
  ////correct for bias in cells with only Y2 using fit betas
  {
    vector[nB_d[1]] b_d1;
    vector[nB_d[2]] b_d2;
    vector[nB_d[3]] b_d3;
    vector[nB_d[4]] b_d4;
    vector[nB_p] b_p;
    b_d1 = beta_d1;
    b_d2 = beta_d2;
    b_d3 = beta_d3;
    b_d4 = beta_d4;
    b_p = beta_p;
    Y2_ds[n2:n3,1] = to_array_1d(Y2[n2:n3,1] + (X_d1[n2:n3,] * b_d1));
    Y2_ds[n2:n3,2] = to_array_1d(Y2[n2:n3,2] + (X_d2[n2:n3,] * b_d2));
    Y2_ds[n2:n3,3] = to_array_1d(Y2[n2:n3,3] + (X_d3[n2:n3,] * b_d3));
    Y2_ds[n2:n3,4] = to_array_1d((Y2[n2:n3,4] + (X_d4[n2:n3,] * b_d4)) 
      .* inv_logit(X_p[n2:n3,] * b_p));
    Y2_ds[n2:n3,5] = to_array_1d((Y2[n2:n3,4] + (X_d4[n2:n3,] * b_d4)) 
      .* (1-inv_logit(X_p[n2:n3,] * b_p)));
  }
  
  //nu to eta
  for(n in 1:(n3)) {
    n_eta[n] = tr_gjam_inv(nu[n])';
  }
}

model {
  
  //covariance priors
  for(j in 1:2) {
    L_Omega[j] ~ lkj_corr_cholesky(8);
  }
 
  //nu priors
  for(l in 1:(L-1)) {
    nu[,l] ~ uniform(-1, 2);
  }
  
  //beta priors
  beta_p ~ normal(0, 1);
  beta_d1 ~ normal(0, 1);
  beta_d2 ~ normal(0, 1);
  beta_d4 ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu[1:n1], L_Sigma[1]);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma[2]);
}

