functions {
  //GJAM transformation to enforce compositional restrictions
  vector tr_gjam_inv(vector nu_i, int D, int d) {
    vector[D] eta_i;
    vector[d] nu_i_p_;
    real nu_i_p;
    real T_i;
    
    eta_i[1:d] = nu_i;
    
    for(j in 1:d) {  
      if(eta_i[j] < 0)  eta_i[j] = 0; 
      if(eta_i[j] < 1)  nu_i_p_[j] = eta_i[j];
      else  nu_i_p_[j] = 1;
    }
    nu_i_p = sum(nu_i_p_);
    
    if(nu_i_p >= 0.99) {
      T_i = (nu_i_p^(-1)) * (1 - (0.01)^(nu_i_p/0.99));
      while(sum(eta_i[1:d]) > 0.99) {
        vector[5] tmp;
        tmp = T_i * eta_i[1:d];
        eta_i[1:d] = tmp;
      }
    }
    eta_i[D] = 1 - sum(eta_i[1:5]);
    return eta_i;
  }
}

data {
  int n1;  //number of cells for GRANIT
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells for NLCD + covariates
  int D;  //number of land cover classes
  int n_b;  //number of covariates (beta) for rho for each LC
  int ri[2*(D-2)];  //beta indexes identifying which LC
  int n_t;  //number of covariates (theta) for p
  vector<lower=0, upper=1>[D-1] Y[n1];  //Y proportions
  matrix<lower=0, upper=1>[n3,D-2] Z;  //Z proportions
  matrix[n1,n_t] X;  //covariates: fitting
  matrix[n3-n1,n_t] X_new;  //covariates: predicting
}

transformed data {
  int d = D-1;
  int tot_b = n_b*(d-1);  //total number of betas
  //QR transformation for X
  real qr_n1 = n1-1;  //avoids cmdstan ambiguity with sqrt(int)
  matrix[n1,n_t] Q = qr_Q(X)[,1:n_t] * sqrt(qr_n1);
  matrix[n_t,n_t] R = qr_R(X)[1:n_t,] / sqrt(qr_n1);
  matrix[n_t,n_t] R_inv = inverse(R);
  //QR transformation for X_new
  real qr_n3 = (n3-n1)-1;  //avoids cmdstan ambiguity with sqrt(int)
  matrix[n3-n1,n_t] Q_new = qr_Q(X_new)[,1:n_t] * sqrt(qr_n3);
  matrix[n_t,n_t] R_new = qr_R(X_new)[1:n_t,] / sqrt(qr_n3);
  matrix[n_t,n_t] R_new_inv = inverse(R_new);
}

parameters {
  cholesky_factor_corr[d] L_Omega[2];  //covariance
  vector<lower=0>[d] L_sigma[2];  //covariance
  vector[d] nu[n1];  //landcover: latent, unconstrained
  vector[n_t] theta_qr;  //thetas (QR decomposition)
  vector[tot_b] beta_z;  //betas (QR decomposition, non-centered)
  real<lower=0> beta_scale;  //betas (QR decomposition, non-centered)
}

transformed parameters {
  matrix[d,d] Sigma[2];  //covariance
  vector[d] Z_[n1];  //Z split & unbiased  
  vector[tot_b] beta_qr = beta_z * beta_scale;  //implies beta ~ N(0,scale)
  
  for(k in 1:2) Sigma[k] = diag_pre_multiply(L_sigma[k], L_Omega[k]);
  {
    vector[n1] p = inv_logit(Q[1:n1,] * theta_qr);
    matrix[n1,d-1] rho;
    for(j in 1:(d-2)) {
      rho[,j] = Q[1:n1,2:n_t] * beta_qr[ri[j+j-1]:ri[j+j]];
      Z_[,j] = to_array_1d(Z[1:n1,j] + rho[,j]);
    }
    rho[,d-1] = Q[,2:n_t] * beta_qr[ri[d+d-3]:ri[d+d-2]];
    Z_[,d-1] = to_array_1d((Z[1:n1,d-1] + rho[,d-1]) .* p);
    Z_[,d] = to_array_1d((Z[1:n1,d-1] + rho[,d-1]) .* (1-p));
  }
}

model {
  for(k in 1:2) {
    L_Omega[k] ~ lkj_corr_cholesky(8);
    L_sigma[k] ~ normal(0, 1);
  }
  for(j in 1:d) nu[,j] ~ normal(0.5, 1);  //will be near observations [0,1]
  theta_qr ~ normal(0, 1);
  beta_z ~ normal(0, 1);
  beta_scale ~ normal(0, 1);
  Y ~ multi_normal_cholesky(nu, Sigma[1]);
  Z_ ~ multi_normal_cholesky(nu, Sigma[2]);
}

generated quantities {
  vector[n1] log_lik;  //log likelihood for model comparison
  simplex[D] eta[n3];  //gjam transformed nu
  vector[d] Z_new_[n3-n1];  //unbiased, split Z
  vector[n_t] theta = R_inv * theta_qr;
  vector[tot_b] beta;
  for(j in 1:(d-1)) {
    beta[ri[j+j-1]:ri[j+j]] = R_inv[2:n_t,2:n_t] * beta_qr[ri[j+j-1]:ri[j+j]];
  }

  {
    vector[n3-n1] p_new = inv_logit(Q_new * theta_qr);
    matrix[n3-n1,d-1] rho_new;
    for(j in 1:(d-2)) {
      rho_new[,j] = Q_new[,2:n_t] * beta_qr[ri[j+j-1]:ri[j+j]];
      Z_new_[,j] = to_array_1d(Z[n2:n3,j] + rho_new[,j]);
    }
    rho_new[,d-1] = Q_new[,2:n_t] * beta_qr[ri[d+d-3]:ri[d+d-2]];
    Z_new_[,d-1] = to_array_1d((Z[n2:n3,d-1] + rho_new[,d-1]) .* p_new);
    Z_new_[,d] = to_array_1d((Z[n2:n3,d-1] + rho_new[,d-1]) .* (1-p_new));
  }
  
  for(n in 1:n1) {
    eta[n] = tr_gjam_inv(nu[n], D, d);
    log_lik[n] = multi_normal_cholesky_lpdf(Y[n] | nu[n], Sigma[1]) +
                 multi_normal_cholesky_lpdf(Z_[n] | nu[n], Sigma[2]);
  }
  for(n in n2:n3) {
    eta[n] = tr_gjam_inv(Z_new_[n-n1], D, d);
  }
}
