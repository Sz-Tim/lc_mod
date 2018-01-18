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
      D_i = (1 - (0.01)^(w_p/0.99)) / w_p;
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
  int L;  //number of land cover classes
  int nB_d;  //number of bias covariates for each LC
  int di[2*(L-2)];  //bias beta indexes
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n1,L-2] Y2;  //NLCD proportions
  //covariates
  matrix[n1,nB_p] X;  //pr(WP|Evg) covariates
}

transformed data {
  int n_beta_d = nB_d*(L-2);  //total number of beta_ds
  //QR decomposition for covariates
  real qr_n = n1-1;  //avoids cmdstan ambiguity with sqrt
  matrix[n1,nB_p] Q = qr_Q(X)[,1:nB_p] * sqrt(qr_n);
  matrix[nB_p,nB_p] R = qr_R(X)[1:nB_p,] / sqrt(qr_n);
  matrix[nB_p,nB_p] R_inv = inverse(R);
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega[2];
  vector<lower=0>[L-1] L_sigma[2];
  //landcover: latent not constrained to be compositional
  vector[L-1] nu[n1];
  //thetas: QR decomposition slopes
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d_z;  //bias betas (QR decomposition)
  real<lower=0> theta_d_scale;  //bias beta scale (so theta_d_z ~ N(0,1))
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];  
  vector[n_beta_d] theta_d;
  matrix[L-1,L-1] L_Sigma[2];
  
  for(j in 1:2) {
    L_Sigma[j] = diag_pre_multiply(L_sigma[j], L_Omega[j]);
  }
  
  theta_d = theta_d_z * theta_d_scale;
  
  //correct bias and split WP to [,4] and Evg to [,5]
  {
    //temporary storage for pWP & bias term
    vector[n1] pWP = inv_logit(Q * theta_p);
    matrix[n1,L-2] bias;
    for(i in 1:4) {
      bias[,i] = Q[,2:nB_p] * theta_d[di[i+i-1]:di[i+i]];
    }
    Y2_[,1] = to_array_1d(Y2[,1] + bias[,1]);
    Y2_[,2] = to_array_1d(Y2[,2] + bias[,2]);
    Y2_[,3] = to_array_1d(Y2[,3] + bias[,3]);
    Y2_[,4] = to_array_1d((Y2[,4] + bias[,4]) .* pWP);
    Y2_[,5] = to_array_1d((Y2[,4] + bias[,4]) .* (1-pWP));
  }
}

model {
  //covariance priors
  for(j in 1:2) {
    L_Omega[j] ~ lkj_corr_cholesky(8);
    L_sigma[j] ~ normal(0, 1);
  }
 
  //nu priors
  for(l in 1:(L-1)) {
    nu[,l] ~ normal(0.5, 1);
  }
  
  //beta priors
  theta_p ~ normal(0, 1);
  theta_d_z ~ normal(0, 1);
  theta_d_scale ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, L_Sigma[1]);
   Y2_ ~ multi_normal_cholesky(nu, L_Sigma[2]);
}

generated quantities {
  //landcover: latent compositional
  simplex[L] n_eta[n1];  //gjam transformed nu
  //betas
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
  vector[n_beta_d] beta_d;  //bias betas
  //log likelihood for model comparison
  vector[n1] log_lik;

  //QR decopmositions
  beta_p = R_inv * theta_p;
  beta_d[di[1]:di[2]] = R_inv[2:nB_p,2:nB_p] * theta_d[di[1]:di[2]];
  beta_d[di[3]:di[4]] = R_inv[2:nB_p,2:nB_p] * theta_d[di[3]:di[4]];
  beta_d[di[5]:di[6]] = R_inv[2:nB_p,2:nB_p] * theta_d[di[5]:di[6]];
  beta_d[di[7]:di[8]] = R_inv[2:nB_p,2:nB_p] * theta_d[di[7]:di[8]];
  
  //enforce compositional constraints
  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(nu[n]);
    log_lik[n] = multi_normal_cholesky_lpdf(Y1[n] | nu[n], L_Sigma[1]) +
                 multi_normal_cholesky_lpdf(Y2_[n] | nu[n], L_Sigma[2]);
  }
}
