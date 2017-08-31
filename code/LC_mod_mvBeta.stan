functions {
  vector tr_gjam_inv(vector w) {
    vector[6] eta;
    vector[5] w_p_;
    real w_p;
    real D_i;
    
    eta[1:5] = w;
    
    for(l in 1:5) {  
      if(eta[l] < 0)  eta[l] = 0; 
      if(eta[l] < 1)  w_p_[l] = eta[l];
      else  w_p_[l] = 1;
    }
    w_p = sum(w_p_);
    
    if(w_p >= 0.99) {
      D_i = (w_p^(-1)) * (1 - (0.01)^(w_p/0.99));
      while(sum(eta[1:5]) > 0.99) {
        vector[5] tmp;
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
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n3,L-2] Y2;  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[n3,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[n3,nB_d[2]] X_d2;  //bias covariates: Oth
  matrix[n3,nB_d[3]] X_d3;  //bias covariates: Hwd
  matrix[n3,nB_d[4]] X_d4;  //bias covariates: Evg
}
transformed data {
  int n_beta_d;  //total number of beta_ds
  // indexes for bias betas
  int d1_2;  //last d1 beta
  int d2_1;  //first d2 beta
  int d2_2;  //last d2 beta
  int d3_1;  //first d3 beta
  int d3_2;  //last d3 beta
  int d4_1;  //first d4 beta
  int d4_2;  //last d4 beta
  
  n_beta_d = sum(nB_d);
  d1_2 = nB_d[1];
  d2_1 = d1_2 + 1;
  d2_2 = d1_2 + nB_d[2];
  d3_1 = d2_2 + 1;
  d3_2 = d2_2 + nB_d[3];
  d4_1 = d3_2 + 1;
  d4_2 = d3_2 + nB_d[4];
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega[2]; //covariance matrix for Y1 & Y2
  vector<lower=0, upper=pi()/2>[L-1] L_sigma_unif[2];  //covariance
  //landcover: latent non-constrained
  vector<lower=-1, upper=2>[L-1] nu[n3];  //latent LC proportions
  //betas
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
  vector[nB_p] mu_p;
  cov_matrix[nB_p] sig_p;
  vector[n_beta_d] beta_d;  //bias betas
  vector[n_beta_d] mu_d;
  cov_matrix[n_beta_d] sig_d;
}

transformed parameters {
  //covariance
  vector<lower=0>[L-1] L_sigma[2];  //covariance matrix for Y1 & Y2
  matrix[L-1,L-1] L_Sigma[2];
  //NLCD de-biasing and splitting
  vector[L-1] Y2_ds[n3];  //unbiased, split NLCD
  
  
  // L_sigma ~ cauchy(0, 2.5) --- reparameterized for speed
  for(j in 1:2) {
    L_sigma[j] = 2.5 * tan(L_sigma_unif[j]);
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
  }
  
  
  //correct bias & split WP to [,4] and Evg to [,5]
  ////fit betas using cells with Y1 & Y2
  Y2_ds[1:n1,1] = to_array_1d(Y2[1:n1,1] 
      + (X_d1[1:n1,] * beta_d[1:d1_2]));
  Y2_ds[1:n1,2] = to_array_1d(Y2[1:n1,2] 
      + (X_d2[1:n1,] * beta_d[d2_1:d2_2]));
  Y2_ds[1:n1,3] = to_array_1d(Y2[1:n1,3] 
      + (X_d3[1:n1,] * beta_d[d3_1:d3_2]));
  Y2_ds[1:n1,4] = to_array_1d((Y2[1:n1,4] 
      + (X_d4[1:n1,] * beta_d[d4_1:d4_2])) 
      .* inv_logit(X_p[1:n1,] * beta_p));
  Y2_ds[1:n1,5] = to_array_1d((Y2[1:n1,4] 
      + (X_d4[1:n1,] * beta_d[d4_1:d4_2])) 
      .* (1-inv_logit(X_p[1:n1,] * beta_p)));
  ////correct for bias in cells with only Y2 using fit betas
  {
    vector[n_beta_d] b_d;
    vector[nB_p] b_p;
    b_d = beta_d;
    b_p = beta_p;
    Y2_ds[n2:n3,1] = to_array_1d(Y2[n2:n3,1] 
        + (X_d1[n2:n3,] * b_d[1:d1_2]));
    Y2_ds[n2:n3,2] = to_array_1d(Y2[n2:n3,2] 
        + (X_d2[n2:n3,] * b_d[d2_1:d2_2]));
    Y2_ds[n2:n3,3] = to_array_1d(Y2[n2:n3,3] 
        + (X_d3[n2:n3,] * b_d[d3_1:d3_2]));
    Y2_ds[n2:n3,4] = to_array_1d((Y2[n2:n3,4] 
        + (X_d4[n2:n3,] * b_d[d4_1:d4_2])) 
        .* inv_logit(X_p[n2:n3,] * b_p));
    Y2_ds[n2:n3,5] = to_array_1d((Y2[n2:n3,4] 
        + (X_d4[n2:n3,] * b_d[d4_1:d4_2])) 
        .* (1-inv_logit(X_p[n2:n3,] * b_p)));
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
  beta_p ~ multi_normal(mu_p, sig_p);
  beta_d ~ multi_normal(mu_d, sig_d);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu[1:n1], L_Sigma[1]);
   Y2_ds ~ multi_normal_cholesky(nu, L_Sigma[2]);
}

generated quantities {
  //landcover: latent compositional
  simplex[L] n_eta[n3];  //gjam transformed nu
    
  //nu to eta
  for(n in 1:(n3)) {
    n_eta[n] = tr_gjam_inv(nu[n]);
  }
}
