functions {
  //GJAM transformation to enforce compositional restrictions
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
        vector[5] tmp = D_i * eta[1:5];
        eta[1:5] = tmp;
      }
    }
    eta[6] = 1 - sum(eta[1:5]);
    return eta;
  }
}

data {
  int n1;  //number of cells for GRANIT
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells for NLCD + covariates
  int L;  //number of land cover classes
  int nB_d;  //number of bias covariates for each LC
  int di[2*(L-2)];  //bias beta indexes
  int nB_p;  //number of covariates for pr(WP|Evg)
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n3,L-2] Y2;  //NLCD proportions
  matrix[n1,nB_p] X;  //covariates
  matrix[n3-n1,nB_p] X_p;  //covariates
}

transformed data {
  int n_beta_d = nB_d*(L-2);  //total number of beta_ds
  real qr_n1 = n1-1;  //avoids cmdstan ambiguity with sqrt
  matrix[n1,nB_p] Q = qr_Q(X)[,1:nB_p] * sqrt(qr_n1);
  matrix[nB_p,nB_p] R = qr_R(X)[1:nB_p,] / sqrt(qr_n1);
  matrix[nB_p,nB_p] R_inv = inverse(R);
  real qr_n3 = (n3-n1)-1;  //avoids cmdstan ambiguity with sqrt
  matrix[n3-n1,nB_p] Q_p = qr_Q(X_p)[,1:nB_p] * sqrt(qr_n3);
  matrix[nB_p,nB_p] R_p = qr_R(X_p)[1:nB_p,] / sqrt(qr_n3);
  matrix[nB_p,nB_p] R_p_inv = inverse(R_p);
}

parameters {
  cholesky_factor_corr[L-1] L_Omega[2];  //landcover: covariance
  vector<lower=0>[L-1] L_sigma[2];  //landcover: covariance
  vector[L-1] nu[n1];  //landcover: latent, unconstrained
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d_std;  //bias betas (QR decomposition, non-centered)
  vector[n_beta_d] theta_d_z;  //bias betas (QR decomposition, non-centered)
  real<lower=0> theta_d_scale;  //bias betas (QR decomposition, non-centered)
}

transformed parameters {
  vector[L-1] Y2_[n1];  //NLCD split & unbiased  
  vector[n_beta_d] theta_d = theta_d_z + theta_d_std * theta_d_scale;
  {
    vector[n1] pWP = inv_logit(Q[1:n1,] * theta_p);
    matrix[n1,L-2] bias;
    for(i in 1:4) {
      bias[,i] = Q[1:n1,2:nB_p] * theta_d[di[i+i-1]:di[i+i]];
    }
    Y2_[,1] = to_array_1d(Y2[1:n1,1] + bias[,1]);
    Y2_[,2] = to_array_1d(Y2[1:n1,2] + bias[,2]);
    Y2_[,3] = to_array_1d(Y2[1:n1,3] + bias[,3]);
    Y2_[,4] = to_array_1d((Y2[1:n1,4] + bias[,4]) .* pWP);
    Y2_[,5] = to_array_1d((Y2[1:n1,4] + bias[,4]) .* (1-pWP));
  }
}

model {
  matrix[L-1,L-1] L_Sigma[2];
  for(j in 1:2) {
    L_Omega[j] ~ lkj_corr_cholesky(8);
    L_sigma[j] ~ normal(0, 1);
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
  }
  for(l in 1:(L-1)) {
    nu[,l] ~ normal(0.5, 1);
  }
  theta_p ~ normal(0, 1);
  theta_d_std ~ normal(0, 1);
  theta_d_z ~ normal(0, 1);
  theta_d_scale ~ normal(0, 1);
  Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
  Y2_ ~ multi_normal_cholesky(nu, L_Sigma[2]);
}

generated quantities {
  simplex[L] n_eta[n3];  //gjam transformed nu
  vector[L-1] Y2_p_[n3-n1];  //unbiased, split NLCD
  vector[nB_p] beta_p = R_inv * theta_p;
  vector[n_beta_d] beta_d;
  beta_d[di[1]:di[2]] = R_p_inv[2:nB_p,2:nB_p] * theta_d[di[1]:di[2]];
  beta_d[di[3]:di[4]] = R_p_inv[2:nB_p,2:nB_p] * theta_d[di[3]:di[4]];
  beta_d[di[5]:di[6]] = R_p_inv[2:nB_p,2:nB_p] * theta_d[di[5]:di[6]];
  beta_d[di[7]:di[8]] = R_p_inv[2:nB_p,2:nB_p] * theta_d[di[7]:di[8]];

  {
    vector[n3-n1] pWP_p = inv_logit(Q_p * theta_p);
    matrix[n3-n1,L-2] bias_p;
    for(i in 1:4) {
      bias_p[,i] = Q_p[,2:nB_p] * theta_d[di[i+i-1]:di[i+i]];
    }
    Y2_p_[,1] = to_array_1d(Y2[n2:n3,1] + bias_p[,1]);
    Y2_p_[,2] = to_array_1d(Y2[n2:n3,2] + bias_p[,2]);
    Y2_p_[,3] = to_array_1d(Y2[n2:n3,3] + bias_p[,3]);
    Y2_p_[,4] = to_array_1d((Y2[n2:n3,4] + bias_p[,4]) .* pWP_p);
    Y2_p_[,5] = to_array_1d((Y2[n2:n3,4] + bias_p[,4]) .* (1-pWP_p));
  }
  
  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(nu[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2_p_[n-n1]);
  }
}
