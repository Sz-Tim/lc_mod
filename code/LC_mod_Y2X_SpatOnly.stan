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
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  vector<lower=0, upper=1>[L-2] Y2[n3];  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  //GPP spatial random effects
  int<lower=1, upper=n3> m;  //number of knots
  matrix[m,m] D_star;  //distance matrix for knots
  matrix[n3,m] D_site_star;  //distance matrix for points to knots
}

transformed data {
  //QR decomposition for covariates
  matrix[n3,nB_p] Q_p = qr_Q(X_p)[,1:nB_p] * sqrt(n3-1);
  matrix[nB_p,nB_p] R_p = qr_R(X_p)[1:nB_p,] / sqrt(n3-1);
  matrix[nB_p,nB_p] R_inv_p = inverse(R_p);
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega; 
  vector<lower=0>[L-1] L_sigma; 
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  //GPP
  real<lower=0> eta[L-1];  //sqrt(GPP variance)
  real<lower=0> sigma[L-1];  //sqrt(nugget)
  real<lower=0> phi[L-1];  //decay rate
  vector[m] w_z[L-1];  //spatial effects
  vector[n3] e_z[L-1];  //
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];
  vector[nB_p] beta_p;
  //GPP
  vector[n3] w[L-1];
  vector[m] w_star[L-1];
  matrix[m,m] Cstar[L-1];
  matrix[n3,m] C_site_star[L-1];
  matrix[n3,m] C_site_inv_Cstar[L-1];
  real eta_sq[L-1];
  real sig_sq[L-1];


  //latent gp at knots
  for(l in 1:(L-1)) {
    eta_sq[l] = pow(eta[l], 2);
    sig_sq[l] = pow(sigma[l], 2);
    Cstar[l] = diag_matrix(rep_vector(eta_sq[l] + sig_sq[l], m));
    for(i in 1:(m-1)) {
	    for(j in (i+1):m) {
		    Cstar[l,i,j] = eta_sq[l] * exp(-D_star[i,j] * phi[l]);
		    Cstar[l,j,i] = Cstar[l,i,j];
	    }
    }
    w_star[l] = cholesky_decompose(Cstar[l]) * w_z[l];
  
    //latent gp at sample locations + bias adjustment from Finley et al. 2009
    C_site_star[l] = eta_sq[l] * exp(-D_site_star * phi[l]);
    // efficient and arithmetically stable for: C_site_star * inverse(Cstar)
    C_site_inv_Cstar[l] = mdivide_right_spd(C_site_star[l], Cstar[l]);
    w[l] = C_site_inv_Cstar[l] * w_star[l] 
          + (e_z[l] 
            .* sqrt(eta_sq[l] + sig_sq[l]
                    - rows_dot_product(C_site_inv_Cstar[l], C_site_star[l])));
  }  
  
  //QR decompositions
  beta_p = R_inv_p * theta_p;

  //split and de-bias Y2
  Y2_[,1] = to_array_1d(to_vector(Y2[1:n1,1]) 
      + w[1,1:n1]);
  Y2_[,2] = to_array_1d(to_vector(Y2[1:n1,2]) 
      + w[2,1:n1]);
  Y2_[,3] = to_array_1d(to_vector(Y2[1:n1,3]) 
      + w[3,1:n1]);
  Y2_[,4] = to_array_1d(to_vector(Y2[1:n1,4]) 
        .* inv_logit(Q_p[1:n1] * theta_p)
      + w[4,1:n1]);
  Y2_[,5] = to_array_1d(to_vector(Y2[1:n1,4]) 
        .* (1 - inv_logit(Q_p[1:n1] * theta_p))
      + w[5,1:n1]);  
}

model {
  //gpp
  sigma ~ normal(0, 1);
  eta ~ normal(0, 1);
  phi ~ normal(0, 1);
  for(l in 1:(L-1)) {
    w_z[l] ~ normal(0, 1);
    e_z[l] ~ normal(0, 1);
  }
  
  //covariance priors
  L_Omega ~ lkj_corr_cholesky(8);
  L_sigma ~ normal(0, 1);
  
  //beta priors
  beta_p ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(Y2_, 
                diag_pre_multiply(L_sigma, L_Omega));
}

generated quantities {
  //landcover: latent compositional
  vector[L-1] Y2new_[n3-n1];
  simplex[L] n_eta[n3];

  Y2new_[,1] = to_array_1d(to_vector(Y2[n2:n3,1])
      + w[1,n2:n3]);
  Y2new_[,2] = to_array_1d(to_vector(Y2[n2:n3,2])
      + w[2,n2:n3]);
  Y2new_[,3] = to_array_1d(to_vector(Y2[n2:n3,3])
      + w[3,n2:n3]);
  Y2new_[,4] = to_array_1d(to_vector(Y2[n2:n3,4])
        .* inv_logit(Q_p[n2:n3] * theta_p)
      + w[4,n2:n3]);
  Y2new_[,5] = to_array_1d(to_vector(Y2[n2:n3,4])
        .* (1 - inv_logit(Q_p[n2:n3] * theta_p))
      + w[5,n2:n3]);

  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(Y2_[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2new_[n-n1]);
  }
}
