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
  //NNGP spatial random effects
  int nn_YX[n3,3];  //index for nn order vs Y&X order
  int<lower=1, upper=n3> M;  //number of neighbors
  int nn_id[M, n3];  //neighbor indices
  matrix[M, n3] nn_d;  //neighbor distances
  matrix[(M*(M-1)/2), n3] nn_dM;  //neighbor distance matrix
  int nn_dim[n3];  //number of neighbors for each cell
  int dim_r;  //number of unique neighborhood distance matrices
  int dim_i[dim_r, 4];  //reference for neighborhood matrix index start & end
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
  //NNGP neighborhood size indices
  int M_i[dim_r] = dim_i[,1];  //neighborhood sizes: number of neighbors
  int M_mx[dim_r] = dim_i[,2];  //neighborhood sizes: unique dist mx's
  int M_1[dim_r] = dim_i[,3];  //neighborhood sizes: index starts
  int M_2[dim_r] = dim_i[,4];  //neighborhood sizes: index ends
  int n_i[dim_r];  //number of cells w/neighborhood dist mx i
  matrix[M,M] nn_dM_i[dim_r];  //unique distance matrices
  for(r in 1:dim_r) {
    int h = 0;
    n_i[r] = M_2[r] - M_1[r] + 1;
    nn_dM_i[r] = diag_matrix(rep_vector(1, M));
    for(j in 1:(nn_dim[M_1[r]]-1)) {
      for(k in (j+1):nn_dim[M_1[r]]) {
        h = h + 1;
        nn_dM_i[r,j,k] = nn_dM[h,M_1[r]];
        nn_dM_i[r,k,j] = nn_dM[h,M_1[r]];
      }
    }
  }
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega; 
  vector<lower=0>[L-1] L_sigma; 
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  vector[n_beta_d] theta_d;  //bias thetas
  //NNGP
  real<lower=0> sigma[L-1];  //sqrt(nugget)
  real<lower=0> phi[L-1];  //decay rate
  vector[n3] w[L-1];  //spatial effects
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];
  vector[nB_p] beta_p;
  vector[n_beta_d] beta_d;
  //NNGP parameters
  matrix[M,n3] um[L-1];
  vector[n3] V[L-1];
  vector[n3] uw_dp[L-1];
  real<lower=0> sig2[L-1];
  
  //NNGP calculations
  for(l in 1:(L-1)) {
    sig2[l] = square(sigma[l]);
    um[l] = exp(-phi[l] * nn_d);
    for(r in 1:dim_r) {
      int gy[n_i[r]] = nn_YX[M_1[r]:M_2[r],1];
      matrix[M_i[r],M_i[r]] exp_nn_dM;
      row_vector[M_i[r]] L_u;
      vector[M_i[r]] v_L;
      matrix[M_i[r],M_i[r]] L_nn;
      
      exp_nn_dM = exp(-phi[l] * block(nn_dM_i[r],1,1,M_i[r],M_i[r]));
      for(j in 1:M_i[r]) exp_nn_dM[j,j] = 1;
      L_nn = cholesky_decompose(exp_nn_dM);
      L_u = mdivide_left_tri_low(L_nn, sub_col(um[l],1,M_1[r],M_i[r]))';
      v_L = mdivide_right_tri_low(L_u, L_nn)';
      V[l,gy] = rep_vector(1 - dot_self(L_u), n_i[r]);
      uw_dp[l,gy] = to_matrix(w[l,to_array_1d(nn_id[1:M_i[r],M_1[r]:M_2[r]])], 
                    n_i[r], M_i[r]) * v_L;
    }
  }
  
  //QR decompositions
  beta_p = R_inv_p * theta_p;
  beta_d[1:d1_2] = R_inv_d1 * theta_d[1:d1_2];
  beta_d[d2_1:d2_2] = R_inv_d2 * theta_d[d2_1:d2_2];
  beta_d[d3_1:d3_2] = R_inv_d3 * theta_d[d3_1:d3_2];
  beta_d[d4_1:d4_2] = R_inv_d4 * theta_d[d4_1:d4_2];

  //split and de-bias Y2
  Y2_[,1] = to_array_1d(to_vector(Y2[1:n1,1]) 
      + (Q_d1[1:n1,] * theta_d[1:d1_2]) 
      + w[1,1:n1]);
  Y2_[,2] = to_array_1d(to_vector(Y2[1:n1,2]) 
      + (Q_d2[1:n1,] * theta_d[d2_1:d2_2])
      + w[2,1:n1]);
  Y2_[,3] = to_array_1d(to_vector(Y2[1:n1,3]) 
      + (Q_d3[1:n1,] * theta_d[d3_1:d3_2])
      + w[3,1:n1]);
  Y2_[,4] = to_array_1d((to_vector(Y2[1:n1,4]) 
          + (Q_d4[1:n1,] * theta_d[d4_1:d4_2])) 
        .* inv_logit(Q_p[1:n1] * theta_p)
      + w[4,1:n1]);
  Y2_[,5] = to_array_1d((to_vector(Y2[1:n1,4]) 
          + (Q_d4[1:n1,] * theta_d[d4_1:d4_2])) 
        .* (1 - inv_logit(Q_p[1:n1] * theta_p))
      + w[5,1:n1]);  
}

model {
  //NNGP
  sigma ~ normal(0, 1);
  phi ~ normal(0, 1);
  for(l in 1:(L-1)) {
    target += -0.5 * (n3*log(sig2[l]) + sum(log(V[l])) 
                      + sum(square(w[l]-uw_dp[l]) ./ V[l]) / sig2[l]);
  }
  
  //covariance priors
  L_Omega ~ lkj_corr_cholesky(8);
  L_sigma ~ normal(0, 1);
  
  //beta priors
  beta_p ~ normal(0, 1);
  beta_d ~ normal(0, 0.1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(Y2_, diag_pre_multiply(L_sigma, L_Omega));
}

generated quantities {
  //landcover: latent compositional
  vector[L-1] Y2new_[n3-n1];
  simplex[L] n_eta[n3];

  Y2new_[,1] = to_array_1d(to_vector(Y2[n2:n3,1])
      + (Q_d1[n2:n3,] * theta_d[1:d1_2])
      + w[1,n2:n3]);
  Y2new_[,2] = to_array_1d(to_vector(Y2[n2:n3,2])
      + (Q_d2[n2:n3,] * theta_d[d2_1:d2_2])
      + w[2,n2:n3]);
  Y2new_[,3] = to_array_1d(to_vector(Y2[n2:n3,3])
      + (Q_d3[n2:n3,] * theta_d[d3_1:d3_2])
      + w[3,n2:n3]);
  Y2new_[,4] = to_array_1d((to_vector(Y2[n2:n3,4])
          + (Q_d4[n2:n3,] * theta_d[d4_1:d4_2])) 
        .* inv_logit(Q_p[n2:n3] * theta_p)
      + w[4,n2:n3]);
  Y2new_[,5] = to_array_1d((to_vector(Y2[n2:n3,4])
          + (Q_d4[n2:n3,] * theta_d[d4_1:d4_2])) 
        .* (1 - inv_logit(Q_p[n2:n3] * theta_p))
      + w[5,n2:n3]);

  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(Y2_[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2new_[n-n1]);
  }
}
