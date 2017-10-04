functions {
  vector tr_gjam_inv(vector w) {
    vector[4] eta;
    vector[3] w_p_;
    real w_p;
    real D_i;
    
    eta[1:3] = w;
    
    for(l in 1:3) {  
      if(eta[l] < 0)  eta[l] = 0; 
      if(eta[l] < 1)  w_p_[l] = eta[l];
      else  w_p_[l] = 1;
    }
    w_p = sum(w_p_);
    
    if(w_p >= 0.99) {
      D_i = (w_p^(-1)) * (1 - (0.01)^(w_p/0.99));
      while(sum(eta[1:3]) > 0.99) {
        vector[3] tmp;
        tmp = D_i * eta[1:3];
        eta[1:3] = tmp;
      }
    }
    eta[4] = 1 - sum(eta[1:3]);
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
  vector<lower=0, upper=1>[L-2] Y2[n3];  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[n3,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[n3,nB_d[2]] X_d2;  //bias covariates: Oth
}

transformed data {
  int n_beta_d = sum(nB_d);  //total number of beta_ds
  // indexes for bias betas
  int d1_2 = nB_d[1];         //last d1 beta
  int d2_1 = d1_2 + 1;        //first d2 beta
  int d2_2 = d1_2 + nB_d[2];  //last d2 beta
  //QR decomposition for covariates
  matrix[n3,nB_p] Q_p = qr_Q(X_p)[,1:nB_p] * sqrt(n3-1);
  matrix[nB_p,nB_p] R_p = qr_R(X_p)[1:nB_p,] / sqrt(n3-1);
  matrix[nB_p,nB_p] R_inv_p = inverse(R_p);
  matrix[n3,nB_d[1]] Q_d1 = qr_Q(X_d1)[,1:nB_d[1]] * sqrt(n3-1);
  matrix[nB_d[1],nB_d[1]] R_d1 = qr_R(X_d1)[1:nB_d[1],] / sqrt(n3-1);
  matrix[nB_d[1],nB_d[1]] R_inv_d1 = inverse(R_d1);
  matrix[n3,nB_d[2]] Q_d2 = qr_Q(X_d2)[,1:nB_d[2]] * sqrt(n3-1);
  matrix[nB_d[2],nB_d[2]] R_d2 = qr_R(X_d2)[1:nB_d[2],] / sqrt(n3-1);
  matrix[nB_d[2],nB_d[2]] R_inv_d2 = inverse(R_d2);
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega; 
  vector<lower=0>[L-1] L_sigma;  
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d;  //bias betas (QR decomposition)
  vector[n_beta_d] beta_d_sig;  //bias beta scale terms
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];
  //betas
  vector[nB_p] beta_p;
  vector[n_beta_d] beta_d_raw;
  vector[n_beta_d] beta_d;

  
  //QR decompositions
  beta_p = R_inv_p * theta_p;
  beta_d_raw[1:d1_2] = R_inv_d1 * theta_d[1:d1_2];
  beta_d_raw[d2_1:d2_2] = R_inv_d2 * theta_d[d2_1:d2_2];
  beta_d = beta_d_raw .* beta_d_sig;
  

  Y2_[,1] = to_array_1d(to_vector(Y2[1:n1,1]) 
      + (Q_d1[1:n1] * theta_d[1:d1_2]));
  Y2_[,2] = to_array_1d((to_vector(Y2[1:n1,2]) 
      + (Q_d2[1:n1] * theta_d[d2_1:d2_2]))
        .* inv_logit(Q_p[1:n1] * theta_p));
  Y2_[,3] = to_array_1d((to_vector(Y2[1:n1,2]) 
      + (Q_d2[1:n1] * theta_d[d2_1:d2_2]))
        .* (1 - inv_logit(Q_p[1:n1] * theta_p)));  
}

model {
  //covariance priors
  L_Omega ~ lkj_corr_cholesky(8);
  L_sigma ~ normal(0, 1);
  
  //beta priors
  beta_p ~ normal(0, 1);
  beta_d_raw ~ normal(0, 1);
  beta_d_sig ~ cauchy(0, 0.5);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(Y2_, diag_pre_multiply(L_sigma, L_Omega));
}

generated quantities {
  //landcover: latent compositional
  vector[L-1] Y2new_[n3-n1];
  simplex[L] n_eta[n3];
  
  Y2new_[,1] = to_array_1d(to_vector(Y2[n2:n3,1]) 
      + (Q_d1[n2:n3] * theta_d[1:d1_2]));
  Y2new_[,2] = to_array_1d((to_vector(Y2[n2:n3,2]) 
      + (Q_d2[n2:n3] * theta_d[d2_1:d2_2]))
        .* inv_logit(Q_p[n2:n3] * theta_p));
  Y2new_[,3] = to_array_1d((to_vector(Y2[n2:n3,2]) 
      + (Q_d2[n2:n3] * theta_d[d2_1:d2_2]))
        .* (1 - inv_logit(Q_p[n2:n3] * theta_p))); 
        
  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(Y2_[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2new_[n-n1]);
  }
}
