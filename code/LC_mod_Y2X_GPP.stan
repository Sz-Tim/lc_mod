functions {
  vector tr_gjam_inv(vector w) {
    vector[6] eta;
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
      vector[5] tmp;
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
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  vector<lower=0, upper=1>[L-2] Y2[n3];  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[n3,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[n3,nB_d[2]] X_d2;  //bias covariates: Oth
  matrix[n3,nB_d[3]] X_d3;  //bias covariates: Hwd
  matrix[n3,nB_d[4]] X_d4;  //bias covariates: Evg
  //GPP spatial random effects
  int<lower=1, upper=n1> m;
  matrix[m,m] D_star;
  matrix[n3,m] D_site_star;
}

transformed data {
  int n_beta_d = sum(nB_d);  //total number of beta_ds
  // indexes for bias betas
  int d1_2 = nB_d[1];         //last d1 beta
  int d2_1 = d1_2 + 1;        //first d2 beta
  int d2_2 = d1_2 + nB_d[2];  //last d2 beta
  int d3_1 = d2_2 + 1;        //first d3 beta
  int d3_2 = d2_2 + nB_d[3];  //last d3 beta
  int d4_1 = d3_2 + 1;        //first d4 beta
  int d4_2 = d3_2 + nB_d[4];  //last d4 beta
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
  matrix[n3,nB_d[3]] Q_d3 = qr_Q(X_d3)[,1:nB_d[3]] * sqrt(n3-1);
  matrix[nB_d[3],nB_d[3]] R_d3 = qr_R(X_d3)[1:nB_d[3],] / sqrt(n3-1);
  matrix[nB_d[3],nB_d[3]] R_inv_d3 = inverse(R_d3);
  matrix[n3,nB_d[4]] Q_d4 = qr_Q(X_d4)[,1:nB_d[4]] * sqrt(n3-1);
  matrix[nB_d[4],nB_d[4]] R_d4 = qr_R(X_d4)[1:nB_d[4],] / sqrt(n3-1);
  matrix[nB_d[4],nB_d[4]] R_inv_d4 = inverse(R_d4);
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega; 
  vector<lower=0, upper=pi()/2>[L-1] L_sigma_unif;  
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d;  //bias betas (QR decomposition)
  //GPP
  real<lower=0> eta[L-2];  //sqrt(GPP variance)
  real<lower=0> sigma[L-2];  //sqrt(nugget)
  real<lower=0> phi[L-2];  //decay rate
  vector[m] w_z[L-2];  //spatial effects
  vector[n3] e_z[L-2];  //
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_ds[n1];
  //betas
  vector[nB_p] beta_p;
  vector[n_beta_d] beta_d;
  //GPP
  vector[n3] w[L-2];
  vector[n3] sigma_e_tilde[L-2];
  matrix[m,m] Cstar[L-2];
  vector[m] w_star[L-2];
  matrix[m,m] inv_Cstar[L-2];
  matrix[n3,m] C_site_star[L-2];
  matrix[n3,m] C_ss_inv_Cstar[L-2];
  real eta_sq[L-2];
  real sig_sq[L-2];


  //latent gp at knots
  for(l in 1:(L-2)) {
    eta_sq[l] = pow(eta[l], 2);
    sig_sq[l] = pow(sigma[l], 2);
    for(i in 1:(m-1)) {
	    for(j in (i+1):m) {
		    Cstar[l,i,j] = eta_sq[l] * exp(-D_star[i,j] * phi[l]);
		    Cstar[l,j,i] = Cstar[l,i,j];
	    }
    }
    for(k in 1:m) {
	    Cstar[l,k,k] = eta_sq[l] + sig_sq[l];
    }
    inv_Cstar[l] = inverse(Cstar[l]);
    w_star[l] = cholesky_decompose(Cstar[l]) * w_z[l];
  
    //latent gp at sample locations
    C_site_star[l] = eta_sq[l] * exp(-D_site_star * phi[l]);
    C_ss_inv_Cstar[l] = C_site_star[l] * inv_Cstar[l];
    w[l] = C_site_star[l] * inv_Cstar[l] * w_star[l];
  
    //bias adjustment from Finley et al. 2009
    sigma_e_tilde[l] = eta_sq[l] + sig_sq[l]
          - rows_dot_product(C_ss_inv_Cstar[l], C_site_star[l]);
    for(i in 1:n3) {
  	  w[l,i] = w[l,i] + e_z[l,i] * sqrt(sigma_e_tilde[l,i]);
    }
  }  
  //QR decompositions
  beta_p = R_inv_p * theta_p;
  beta_d[1:d1_2] = R_inv_d1 * theta_d[1:d1_2];
  beta_d[d2_1:d2_2] = R_inv_d2 * theta_d[d2_1:d2_2];
  beta_d[d3_1:d3_2] = R_inv_d3 * theta_d[d3_1:d3_2];
  beta_d[d4_1:d4_2] = R_inv_d4 * theta_d[d4_1:d4_2];
  

  Y2_ds[,1] = to_array_1d(to_vector(Y2[1:n1,1]) 
      + (Q_d1[1:n1] * theta_d[1:d1_2])
      + w[1,1:n1]);
  Y2_ds[,2] = to_array_1d(to_vector(Y2[1:n1,2]) 
      + (Q_d2[1:n1] * theta_d[d2_1:d2_2])
      + w[2,1:n1]);
  Y2_ds[,3] = to_array_1d(to_vector(Y2[1:n1,3]) 
      + (Q_d3[1:n1] * theta_d[d3_1:d3_2]) 
      + w[3,1:n1]);
  Y2_ds[,4] = to_array_1d((to_vector(Y2[1:n1,4]) 
      + (Q_d4[1:n1] * theta_d[d4_1:d4_2]))
        .* inv_logit(Q_p[1:n1] * theta_p)
      + w[4,1:n1]);
  Y2_ds[,5] = to_array_1d((to_vector(Y2[1:n1,4]) 
      + (Q_d4[1:n1] * theta_d[d4_1:d4_2]))
        .* (1 - inv_logit(Q_p[1:n1] * theta_p))
      + w[4,1:n1]);  
}

model {
  //gpp
  sigma ~ normal(0, 1);
  eta ~ normal(0, 1);
  phi ~ normal(0, 5);
  for(l in 1:(L-2)) {
    w_z[l] ~ normal(0, 1);
    e_z[l] ~ normal(0, 1);
  }
  
  //covariance priors
  L_Omega ~ lkj_corr_cholesky(8);
  
  //beta priors
  beta_p ~ normal(0, 1);
  beta_d ~ normal(0, 0.1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(Y2_ds, 
                diag_pre_multiply(2.5 * tan(L_sigma_unif), L_Omega));
}

generated quantities {
  //landcover: latent compositional
  vector[L-1] Y2_ds_new[n3-n1];
  simplex[L] n_eta[n3];

  Y2_ds_new[,1] = to_array_1d(to_vector(Y2[n2:n3,1])
      + (Q_d1[n2:n3] * theta_d[1:d1_2])
      + w[1,n2:n3]);
  Y2_ds_new[,2] = to_array_1d(to_vector(Y2[n2:n3,2])
      + (Q_d2[n2:n3] * theta_d[d2_1:d2_2])
      + w[2,n2:n3]);
  Y2_ds_new[,3] = to_array_1d(to_vector(Y2[n2:n3,3])
      + (Q_d3[n2:n3] * theta_d[d3_1:d3_2])
      + w[3,n2:n3]);
  Y2_ds_new[,4] = to_array_1d((to_vector(Y2[n2:n3,4])
      + (Q_d4[n2:n3] * theta_d[d4_1:d4_2]))
        .* inv_logit(Q_p[n2:n3] * theta_p)
      + w[4,n2:n3]);
  Y2_ds_new[,5] = to_array_1d((to_vector(Y2[n2:n3,4])
      + (Q_d4[n2:n3] * theta_d[d4_1:d4_2]))
        .* (1 - inv_logit(Q_p[n2:n3] * theta_p))
      + w[4,n2:n3]);

  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(Y2_ds[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2_ds_new[n-n1]);
  }
}
