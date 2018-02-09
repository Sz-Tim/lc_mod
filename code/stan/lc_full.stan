functions {
  //GJAM transformation to enforce compositional restrictions
  vector tr_gjam_inv(vector w, int D, int d) {
    vector[D] eta;
    vector[d] w_p_;
    real w_p;
    real D_i;
    
    eta[1:d] = w;
    
    for(l in 1:d) {  
      if(eta[l] < 0)  eta[l] = 0; 
      if(eta[l] < 1)  w_p_[l] = eta[l];
      else  w_p_[l] = 1;
    }
    w_p = sum(w_p_);
    
    if(w_p >= 0.99) {
      D_i = (w_p^(-1)) * (1 - (0.01)^(w_p/0.99));
      while(sum(eta[1:d]) > 0.99) {
        vector[d] tmp = D_i * eta[1:d];
        eta[1:d] = tmp;
      }
    }
    eta[D] = 1 - sum(eta[1:d]);
    return eta;
  }
}

data {
  int n1;  //number of cells for Y1
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells for Y2, X, and Z
  int D;  //number of categories
  int nB;  //number of bias covariates for each category
  int ri[2*(D-2)];  //bias beta indexes
  int nTh;  //number of covariates for pr(Y2'_{d-1} | Y2_{d-1})
  vector<lower=0, upper=1>[D-1] Y1[n1];  //Y1 proportions
  matrix<lower=0, upper=1>[n3,D-2] Y2;  //Y2 proportions
  matrix[n3,nTh] Z;  //covariates
}

transformed data {
  int d = D-1;
  int n_beta = nB*(d-1);  //total number of betas
  real qr_n = n3-1;  //avoids cmdstan ambiguity with sqrt
  matrix[n3,nTh] QZ = qr_Q(Z)[,1:nTh] * sqrt(qr_n);
  matrix[nTh,nTh] RZ = qr_R(Z)[1:nTh,] / sqrt(qr_n);
  matrix[nTh,nTh] RZ_inv = inverse(RZ);
}

parameters {
  cholesky_factor_corr[d] L_Omega[2];  //landcover: covariance
  vector<lower=0>[d] L_sigma[2];  //landcover: covariance
  vector[d] nu[n1];  //landcover: latent, unconstrained
  vector[nTh] theta_qr;  //p thetas (QR decomposition)
  vector[n_beta] beta_qr_std;  //rho betas N(0,1): non-centered parameterization
  vector[n_beta] beta_qr_ctr;  //rho betas (QR decomposition, centered)
  real<lower=0> beta_qr_scale;  //rho betas (QR decomposition, scale)
}

transformed parameters {
  vector[d] Y2_[n1];  //Y2 split & unbiased  
  vector[n_beta] beta_qr = beta_qr_ctr + beta_qr_std * beta_qr_scale;
  {
    vector[n1] p = inv_logit(QZ[1:n1,] * theta_qr);
    matrix[n1,d-1] rho;
    for(i in 1:(d-2)) {
      rho[,i] = QZ[1:n1,2:nTh] * beta_qr[ri[i+i-1]:ri[i+i]];
      Y2_[,i] = to_array_1d(Y2[1:n1,i] + rho[,i]);
    }
    rho[,d-1] = QZ[1:n1,2:nTh] * beta_qr[ri[d+d-3]:ri[d+d-2]];
    Y2_[,d-1] = to_array_1d((Y2[1:n1,d-1] + rho[,d-1]) .* p);
    Y2_[,d] = to_array_1d((Y2[1:n1,d-1] + rho[,d-1]) .* (1-p));
  }
}

model {
  matrix[d,d] L_Sigma[2];
  for(j in 1:2) {
    L_Omega[j] ~ lkj_corr_cholesky(4);
    L_sigma[j] ~ normal(0, 1);
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
  }
  for(l in 1:d) {
    nu[,l] ~ normal(0.5, 1);
  }
  theta_qr ~ normal(0, 1);
  beta_qr_std ~ normal(0, 1); 
  beta_qr_ctr ~ normal(0, 1);
  beta_qr_scale ~ normal(0, 1);
  Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
  Y2_ ~ multi_normal_cholesky(nu, L_Sigma[2]);
}

generated quantities {
  simplex[D] n_eta[n3];  //gjam transformed nu
  vector[d] Y2_pred_[n3-n1];  //unbiased, split Y2
  vector[nTh] theta = RZ_inv * theta_qr;
  vector[n_beta] beta;
  vector[n1] log_lik;

  for(i in 1:(d-1)) {
    beta[ri[i+i-1]:ri[i+i]] = RZ_inv[2:nTh,2:nTh] * beta_qr[ri[i+i-1]:ri[i+i]];
  }
  {
    vector[n3-n1] p_pred = inv_logit(QZ[n2:n3,] * theta_qr);
    matrix[n3-n1,d-1] rho_pred;
    for(i in 1:(d-2)) {
      rho_pred[,i] = QZ[n2:n3,2:nTh] * beta_qr[ri[i+i-1]:ri[i+i]];
      Y2_pred_[,i] = to_array_1d(Y2[n2:n3,i] + rho_pred[,i]);
    }
    rho_pred[,d-1] = QZ[n2:n3,2:nTh] * beta_qr[ri[d+d-3]:ri[d+d-2]];
    Y2_pred_[,d-1] = to_array_1d((Y2[n2:n3,d-1] + rho_pred[,d-1]) .* p_pred);
    Y2_pred_[,d] = to_array_1d((Y2[n2:n3,d-1] + rho_pred[,d-1]) .* (1-p_pred));
  }
  
  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(nu[n], D, d);
    log_lik[n] = multi_normal_cholesky_lpdf(Y1[n] | nu[n], L_Sigma[1]) +
                 multi_normal_cholesky_lpdf(Y2_[n] | nu[n], L_Sigma[2]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2_pred_[n-n1], D, d);
  }
}
