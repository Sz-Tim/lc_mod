functions {
  //GJAM transformation to enforce compositional restrictions
  vector tr_gjam_inv(vector PSI_i_, int K, int k) {
    // Final compositions are constrained to [0,1], and the sum across all 
    //categories must equal 1. This transformation takes the unconstrained 
    //values in PSI_ and imposes these constraints by truncating negative values 
    //to 0, truncating values >1 to 1, and compressing the proportions to the 
    //degree necessary to sum to 1. The input contains k categories, and so the
    //final category K is calculated as 1 - sum(categories 1:k). This function
    //is adapted from Clark et al. 2017.
    vector[K] PSI_i; //final compositional output
    vector[k] PSI_i_trunc; //truncated but uncompressed proportions
    real PSI_i_sum; //sum of truncated but unconstrained proportions
    real T_i; //compression factor if necessary
    
    PSI_i[1:k] = PSI_i_;
    //truncate values to [0,1]
    for(j in 1:k) {  
      //if unconstrained category j is negative, set to 0
      if(PSI_i[j] < 0)  PSI_i[j] = 0; 
      //if unconstrained category j is 0-1, store value
      if(PSI_i[j] < 1)  PSI_i_trunc[j] = PSI_i[j];
      //if unconstrained category j is >1, set to 1
      else  PSI_i_trunc[j] = 1;
    }
    //calculate sum of truncated proportions
    PSI_i_sum = sum(PSI_i_trunc);
    //compress if necessary
    if(PSI_i_sum >= 0.99) {
      T_i = (PSI_i_sum^(-1)) * (1 - (0.01)^(PSI_i_sum/0.99));
      while(sum(PSI_i[1:k]) > 0.99) {
        vector[k] tmp;
        tmp = T_i * PSI_i[1:k];
        PSI_i[1:k] = tmp;
      }
    }
    //calculate final category K
    PSI_i[K] = 1 - sum(PSI_i[1:k]);
    return PSI_i;
  }
}


data {
  int n1;  //number of cells for Y
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells for Z and covariates
  int K;  //number of land cover classes
  int nX;  //number of delta covariates for each LC
  int di[2*(K-2)];  //beta indexes identifying each LC
  int nV;  //number of rho covariates
  vector<lower=0, upper=1>[K-1] Y[n1];  //Y proportions: training
  matrix<lower=0, upper=1>[n3,K-2] Z;  //Z proportions
  vector[K] new_Y[n3-n1];  //Y proportions: testing
  matrix[n1,nV] V;  //covariates: fitting
  matrix[n3-n1,nV] new_V;  //covariates: predicting
}


transformed data {
  int k = K-1;
  int tot_b = nX*(k-1);  //total number of betas
  
  //QR decomposition of V (for efficiency)
  real qr_n1 = n1-1;  
  matrix[n1,nV] Q = qr_Q(V)[,1:nV] * sqrt(qr_n1);
  matrix[nV,nV] R = qr_R(V)[1:nV,] / sqrt(qr_n1);
  matrix[nV,nV] R_inv = inverse(R);
  
  //QR decomposition of new_V (for efficiency)
  real qr_n3 = (n3-n1)-1;  
  matrix[n3-n1,nV] new_Q = qr_Q(new_V)[,1:nV] * sqrt(qr_n3);
  matrix[nV,nV] new_R = qr_R(new_V)[1:nV,] / sqrt(qr_n3);
  matrix[nV,nV] new_R_inv = inverse(new_R);
}


parameters {
  cholesky_factor_corr[k] L_Omega[2];  //covariance
  vector<lower=0>[k] L_sigma[2];  //covariance
  vector[k] PSI_[n1];  //landcover: latent, unconstrained (PSI_ = PSI' in ms)
  vector[nV] theta_qr;  //thetas (QR decomposition)
  vector[tot_b] beta_z;  //betas (QR decomposition, non-centered)
  real<lower=0> beta_scale;  //betas (QR decomposition, non-centered)
}


transformed parameters {
  matrix[k,k] Sigma[2];  //covariance matrices
  vector[k] Z_[n1];  //split & unbiased Z for cells with Y  (Z_ = Z' in ms)
  vector[tot_b] beta_qr = beta_z * beta_scale;  //implies beta ~ N(0,scale)
  
  //Covariance matrix for Y, Z
  for(l in 1:2) Sigma[l] = diag_pre_multiply(L_sigma[l], L_Omega[l]);
  
  {
    //rho: pr(WP|Evg) = inv_logit(V * theta) -- QR decomposition for efficiency
    //delta: discrepancy btw Y, Z = X * beta -- QR decomposition for efficiency
    //Z' = Z + delta   [Open, Other, Deciduous]
    //Z' = (Z + delta) * rho   [WP]
    //Z' = (Z + delta) * (1 - rho)   [Evg]
    vector[n1] rho = inv_logit(Q[,1:nV] * theta_qr);
    matrix[n1,k-1] delta;
    for(j in 1:(k-2)) {
      delta[,j] = Q[,1:nX] * beta_qr[di[j+j-1]:di[j+j]];
      Z_[,j] = to_array_1d(Z[1:n1,j] + delta[,j]);
    }
    delta[,k-1] = Q[,1:nX] * beta_qr[di[k+k-3]:di[k+k-2]];
    Z_[,k-1] = to_array_1d((Z[1:n1,k-1] + delta[,k-1]) .* rho);
    Z_[,k] = to_array_1d((Z[1:n1,k-1] + delta[,k-1]) .* (1-rho));
  }
}


model {
  for(l in 1:2) {
    L_Omega[l] ~ lkj_corr_cholesky(8);
    L_sigma[l] ~ normal(0, 1);
  }
  for(j in 1:k) PSI_[,j] ~ normal(0.5, 1);  //will be near observations [0,1]
  theta_qr ~ normal(0, 1);
  beta_z ~ normal(0, 1);
  beta_scale ~ normal(0, 1);
  Y ~ multi_normal_cholesky(PSI_, Sigma[1]);
  Z_ ~ multi_normal_cholesky(PSI_, Sigma[2]);
}


generated quantities {
  vector[n3-n1] log_lik;  //log likelihood for model comparison
  vector[k] new_PSI_[n3-n1]; //latent, unconstrained for predictions
  simplex[K] PSI[n3];  //transformed, compositional PSI_ 
  vector[k] new_Z_[n3-n1];  //unbiased, split Z for cells without Y
  vector[nV] theta = R_inv * theta_qr; //theta on natural scale -- reverse QR decomposition
  vector[tot_b] beta; //beta on natural scale -- reverse QR decomposition
  vector[K] PSI_Y_error[n3-n1];  //prediction error for testing set
  
  for(j in 1:(k-1)) {
    beta[di[j+j-1]:di[j+j]] = R_inv[1:nX,1:nX] * beta_qr[di[j+j-1]:di[j+j]];
  }

  //calculate Z' for cells without Y using fitted betas & thetas
  //rho, delta, Z' are calculated as above in "transformed parameters" block
  {
    vector[n3-n1] new_rho = inv_logit(new_Q[,1:nV] * theta_qr);
    matrix[n3-n1,k-1] new_delta;
    for(j in 1:(k-2)) {
      new_delta[,j] = new_Q[,1:nX] * beta_qr[di[j+j-1]:di[j+j]];
      new_Z_[,j] = to_array_1d(Z[n2:n3,j] + new_delta[,j]);
    }
    new_delta[,k-1] = new_Q[,1:nX] * beta_qr[di[k+k-3]:di[k+k-2]];
    new_Z_[,k-1] = to_array_1d((Z[n2:n3,k-1] + new_delta[,k-1]) .* new_rho);
    new_Z_[,k] = to_array_1d((Z[n2:n3,k-1] + new_delta[,k-1]) .* (1-new_rho));
  }
  
  //calculate PSI from PSI' for cells with Y and Z
  for(n in 1:n1) PSI[n] = tr_gjam_inv(PSI_[n], K, k);
  
  //calculate PSI from PSI' for cells without Y
  for(n in n2:n3) {
    new_PSI_[n-n1] = multi_normal_cholesky_rng(new_Z_[n-n1], Sigma[2]);
    PSI[n] = tr_gjam_inv(new_PSI_[n-n1], K, k);
    PSI_Y_error[n-n1] = PSI[n] - new_Y[n-n1];
    log_lik[n-n1] = multi_normal_cholesky_lpdf(new_Y[n-n1][1:k] | new_PSI_[n-n1], Sigma[1]);
  }
}
