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
  int nB_d;  //number of bias covariates for each LC
  int di[2*(L-2)];  //bias beta indexes
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n3,L-2] Y2;  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X;  //pr(WP|Evg) covariates
}

transformed data {
  int n_beta_d = nB_d*(L-2);  //total number of beta_ds
  //QR decomposition for covariates
  matrix[n3,nB_p] Q = qr_Q(X)[,1:nB_p] * sqrt(n3-1);
  matrix[nB_p,nB_p] R = qr_R(X)[1:nB_p,] / sqrt(n3-1);
  matrix[nB_p,nB_p] R_inv = inverse(R);
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega;
  vector<lower=0>[L-1] L_sigma;
  //thetas: QR decomposition slopes
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d_z;  //bias betas (QR decomposition)
  real<lower=0> theta_d_scale;  //bias beta scale (so theta_d_z ~ N(0,1))
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];  
  vector[n_beta_d] theta_d;
  
  theta_d = theta_d_z * theta_d_scale;
  
  //correct bias and split WP to [,4] and Evg to [,5]
  {
    //temporary storage for pWP & bias term
    vector[n1] pWP = inv_logit(Q[1:n1,] * theta_p);
    matrix[n1,L-2] bias;
    for(i in 1:4) {
      bias[,i] = Q[1:n1,2:nB_p] * theta_d[di[i+i-1]:di[i+i]];
    }
    Y2_[1:n1,1] = to_array_1d(Y2[1:n1,1] + bias[,1]);
    Y2_[1:n1,2] = to_array_1d(Y2[1:n1,2] + bias[,2]);
    Y2_[1:n1,3] = to_array_1d(Y2[1:n1,3] + bias[,3]);
    Y2_[1:n1,4] = to_array_1d((Y2[1:n1,4] + bias[,4]) .* pWP);
    Y2_[1:n1,5] = to_array_1d((Y2[1:n1,4] + bias[,4]) .* (1-pWP));
  }
}

model {
  //covariance priors
  L_Omega ~ lkj_corr_cholesky(8);
  L_sigma ~ normal(0, 1);
  
  //beta priors
  theta_p ~ normal(0, 1);
  theta_d_z ~ normal(0, 1);
  theta_d_scale ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(Y2_, diag_pre_multiply(L_sigma, L_Omega));
}

generated quantities {
  //landcover: latent compositional
  simplex[L] n_eta[n3];  //gjam transformed nu
  //betas
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
  vector[n_beta_d] beta_d;  //bias betas
  //log likelihood for model comparison
  vector[n1] log_lik;

  //QR decopmositions
  beta_p = R_inv * theta_p;
  beta_d[di[1]:di[2]] = R_inv[2:nB_p, 2:nB_p] * theta_d[di[1]:di[2]];
  beta_d[di[3]:di[4]] = R_inv[2:nB_p, 2:nB_p] * theta_d[di[3]:di[4]];
  beta_d[di[5]:di[6]] = R_inv[2:nB_p, 2:nB_p] * theta_d[di[5]:di[6]];
  beta_d[di[7]:di[8]] = R_inv[2:nB_p, 2:nB_p] * theta_d[di[7]:di[8]];
  
  //calculate bias and pWP
  {
	vector[L-1] Y2new_[n3-n1];  //unbiased, split NLCD
	matrix[n3-n1,L-2] bias_new;
    vector[n3-n1] pWP_new;
    for(i in 1:4) {
      bias_new[,i] = Q[n2:n3,2:nB_p] * theta_d[di[i+i-1]:di[i+i]];
    }
    pWP_new = inv_logit(Q[n2:n3,] * theta_p);
    
    //correct bias and split WP to [,4] and Evg to [,5]
    Y2new_[,1] = to_array_1d(Y2[n2:n3,1] + bias_new[,1]);
    Y2new_[,2] = to_array_1d(Y2[n2:n3,2] + bias_new[,2]);
    Y2new_[,3] = to_array_1d(Y2[n2:n3,3] + bias_new[,3]);
    Y2new_[,4] = to_array_1d((Y2[n2:n3,4] + bias_new[,4]) .* pWP_new);
    Y2new_[,5] = to_array_1d((Y2[n2:n3,4] + bias_new[,4]) .* (1-pWP_new));
  
	//enforce compositional constraints
	for(n in 1:n1) {
		n_eta[n] = tr_gjam_inv(Y2_[n]);
		log_lik[n] = multi_normal_cholesky_lpdf(Y1[n] | Y2_[n], 
                                    diag_pre_multiply(L_sigma, L_Omega));
	}
	for(n in n2:n3) {
		n_eta[n] = tr_gjam_inv(Y2new_[n-n1]);
	}
  }
}
