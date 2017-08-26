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
  cholesky_factor_corr[L-1] L_Omega_Y1; 
  vector<lower=0, upper=pi()/2>[L-1] L_sigma_unif_Y1;  
  cholesky_factor_corr[L-2] L_Omega_Y2; 
  vector<lower=0, upper=pi()/2>[L-2] L_sigma_unif_Y2;  
  //landcover: latent unconstrained
  vector<lower=-1, upper=2>[L-2] nu_Y2[n3];  
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d;  //bias betas (QR decomposition)
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] nu_Y1[n3];
  vector[n3] pX;
  vector[n3] nu4_d;
  //betas
  vector[nB_p] beta_p;
  vector[n_beta_d] beta_d;

  
  //QR decompositions
  beta_p = R_inv_p * theta_p;
  beta_d[1:d1_2] = R_inv_d1 * theta_d[1:d1_2];
  beta_d[d2_1:d2_2] = R_inv_d2 * theta_d[d2_1:d2_2];
  beta_d[d3_1:d3_2] = R_inv_d3 * theta_d[d3_1:d3_2];
  beta_d[d4_1:d4_2] = R_inv_d4 * theta_d[d4_1:d4_2];
  
  
  //estimate bias & split WP [,4] and Evg [,5]
  ////fit betas using cells with Y1 & Y2
  pX[1:n1] = inv_logit(Q_p[1:n1,] * theta_p);
  nu4_d[1:n1] = to_vector(nu_Y2[1:n1,4]) + (Q_d4[1:n1,] * theta_d[d4_1:d4_2]);
  nu_Y1[1:n1,1] = to_array_1d(to_vector(nu_Y2[1:n1,1]) 
        + (Q_d1[1:n1,] * theta_d[1:d1_2]));
  nu_Y1[1:n1,2] = to_array_1d(to_vector(nu_Y2[1:n1,2])
        + (Q_d2[1:n1,] * theta_d[d2_1:d2_2]));
  nu_Y1[1:n1,3] = to_array_1d(to_vector(nu_Y2[1:n1,3]) 
        + (Q_d3[1:n1,] * theta_d[d3_1:d3_2]));
  nu_Y1[1:n1,4] = to_array_1d(nu4_d[1:n1] ./ pX[1:n1]);
  nu_Y1[1:n1,5] = to_array_1d(nu4_d[1:n1] ./ (1 - pX[1:n1]));
  ////predict bias in cells with only Y2 using fit betas
  {
    vector[nB_p] b_p;
    vector[n_beta_d] b_d;
    b_p = theta_p;
    b_d = theta_d;
    pX[n2:n3] = inv_logit(Q_p[n2:n3,] * b_p);
    nu4_d[n2:n3] = to_vector(nu_Y2[n2:n3,4]) + (Q_d4[n2:n3,] * b_d[d4_1:d4_2]);
    nu_Y1[n2:n3,1] = to_array_1d(to_vector(nu_Y2[n2:n3,1])
          + (Q_d1[n2:n3,] * b_d[1:d1_2]));
    nu_Y1[n2:n3,2] = to_array_1d(to_vector(nu_Y2[n2:n3,2])
          + (Q_d2[n2:n3,] * b_d[d2_1:d2_2]));
    nu_Y1[n2:n3,3] = to_array_1d(to_vector(nu_Y2[n2:n3,3])
          + (Q_d3[n2:n3,] * b_d[d3_1:d3_2]));
    nu_Y1[n2:n3,4] = to_array_1d(nu4_d[n2:n3] ./ pX[n2:n3]);
    nu_Y1[n2:n3,5] = to_array_1d(nu4_d[n2:n3] ./ (1 - pX[n2:n3]));
  }
}

model {
  //covariance priors
  L_Omega_Y1 ~ lkj_corr_cholesky(8);
  L_Omega_Y2 ~ lkj_corr_cholesky(8);
 
  //nu priors
  for(l in 1:(L-2)) {
    nu_Y2[,l] ~ uniform(-1, 2);
  }
  
  //beta priors
  beta_p ~ normal(0, 1);
  beta_d ~ normal(0, 0.1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu_Y1[1:n1], 
                diag_pre_multiply(2.5 * tan(L_sigma_unif_Y1), L_Omega_Y1));
   Y2 ~ multi_normal_cholesky(nu_Y2,
                diag_pre_multiply(2.5 * tan(L_sigma_unif_Y2), L_Omega_Y2));
}

generated quantities {
  //landcover: latent compositional
  simplex[L] n_eta[n3];
  for(n in 1:(n3)) {
    n_eta[n] = tr_gjam_inv(nu_Y1[n]);
  }
}