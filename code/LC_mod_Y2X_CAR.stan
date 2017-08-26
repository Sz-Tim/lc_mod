functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
  
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
  //adjacency matrix
  matrix<lower=0, upper=1>[n3,n3] W;  //adjacency matrix: Y1 + Y2
  int W_n;  // number of adjacent pairs: Y1 + Y2
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
  //CAR
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[n3] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n3] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n3 - 1)) {
      for (j in (i + 1):n3) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:n3) D_sparse[i] = sum(W[i]);
  {
    vector[n3] invsqrtD;  
    for (i in 1:n3) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega; 
  vector<lower=0, upper=pi()/2>[L-1] L_sigma_unif;  
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d;  //bias betas (QR decomposition)
  //CAR
  matrix[n3, L-1] phi; 
  real<lower=0> tau[L-1];
  real<lower=0, upper=1> alpha[L-1];
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_ds[n1];
  //betas
  vector[nB_p] beta_p;
  vector[n_beta_d] beta_d;

  
  //QR decompositions
  beta_p = R_inv_p * theta_p;
  beta_d[1:d1_2] = R_inv_d1 * theta_d[1:d1_2];
  beta_d[d2_1:d2_2] = R_inv_d2 * theta_d[d2_1:d2_2];
  beta_d[d3_1:d3_2] = R_inv_d3 * theta_d[d3_1:d3_2];
  beta_d[d4_1:d4_2] = R_inv_d4 * theta_d[d4_1:d4_2];
  

  Y2_ds[,1] = to_array_1d(to_vector(Y2[1:n1,1]) 
      + (Q_d1[1:n1] * theta_d[1:d1_2]) 
      + phi[1:n1,1]);
  Y2_ds[,2] = to_array_1d(to_vector(Y2[1:n1,2]) 
      + (Q_d2[1:n1] * theta_d[d2_1:d2_2]) 
      + phi[1:n1,2]);
  Y2_ds[,3] = to_array_1d(to_vector(Y2[1:n1,3]) 
      + (Q_d3[1:n1] * theta_d[d3_1:d3_2]) 
      + phi[1:n1,3]);
  Y2_ds[,4] = to_array_1d((to_vector(Y2[1:n1,4]) 
      + (Q_d4[1:n1] * theta_d[d4_1:d4_2]))
        .* inv_logit(Q_p[1:n1] * theta_p) 
      + phi[1:n1,4]);
  Y2_ds[,5] = to_array_1d((to_vector(Y2[1:n1,4]) 
      + (Q_d4[1:n1] * theta_d[d4_1:d4_2]))
        .* (1 - inv_logit(Q_p[1:n1] * theta_p)) 
      + phi[1:n1,5]);  
}

model {
  //CAR
  for(i in 1:(L-1)) {
    phi[, i] ~ sparse_car(tau[i], alpha[i], 
                              W_sparse, D_sparse, lambda, n3, W_n);
  }
  tau ~ gamma(2,2);
  
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
      + phi[n2:n3,1]);
  Y2_ds_new[,2] = to_array_1d(to_vector(Y2[n2:n3,2])
      + (Q_d2[n2:n3] * theta_d[d2_1:d2_2]) 
      + phi[n2:n3,2]);
  Y2_ds_new[,3] = to_array_1d(to_vector(Y2[n2:n3,3])
      + (Q_d3[n2:n3] * theta_d[d3_1:d3_2]) 
      + phi[n2:n3,3]);
  Y2_ds_new[,4] = to_array_1d((to_vector(Y2[n2:n3,4])
      + (Q_d4[n2:n3] * theta_d[d4_1:d4_2]))
        .* inv_logit(Q_p[n2:n3] * theta_p) 
      + phi[n2:n3,4]);
  Y2_ds_new[,5] = to_array_1d((to_vector(Y2[n2:n3,4])
      + (Q_d4[n2:n3] * theta_d[d4_1:d4_2]))
        .* (1 - inv_logit(Q_p[n2:n3] * theta_p)) 
      + phi[n2:n3,5]);

  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(Y2_ds[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2_ds_new[n-n1]);
  }
}