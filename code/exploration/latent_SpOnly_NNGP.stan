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
  int n1;  //number of cells with Y1 & Y2
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells with only Y2
  int L;  //number of land cover classes
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n3,L-2] Y2;  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
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

  //QR decomposition for covariates
  matrix[n3,nB_p] Q_p = qr_Q(X_p)[,1:nB_p] * sqrt(n3-1);
  matrix[nB_p,nB_p] R_p = qr_R(X_p)[1:nB_p,] / sqrt(n3-1);
  matrix[nB_p,nB_p] R_inv_p = inverse(R_p);
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
  cholesky_factor_corr[L-1] L_Omega[2];
  vector<lower=0>[L-1] L_sigma[2];
  //landcover: latent not constrained to be compositional
  vector[L-1] nu[n1];
  //thetas: QR decomposition slopes
  vector[nB_p] theta_p;  //pr(WP|Evg) thetas
  //NNGP
  real<lower=0> sigma[L-1];  //sqrt(nugget)
  real<lower=0> phi[L-1];  //decay rate
  vector[n3] w[L-1];  //spatial effects
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];  
  //betas
  vector[nB_p] beta_p;  //pr(WP|Evg) betas
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
  
  //QR decopmositions
  beta_p = R_inv_p * theta_p;
  
  //split and de-bias Y2
  Y2_[1:n1,1] = to_array_1d(Y2[1:n1,1] 
      + w[1,1:n1]);
  Y2_[1:n1,2] = to_array_1d(Y2[1:n1,2] 
      + w[2,1:n1]);
  Y2_[1:n1,3] = to_array_1d(Y2[1:n1,3] 
      + w[3,1:n1]);
  Y2_[1:n1,4] = to_array_1d((Y2[1:n1,4]) 
        .* inv_logit(Q_p[1:n1,] * theta_p)
      + w[4,1:n1]);
  Y2_[1:n1,5] = to_array_1d((Y2[1:n1,4]) 
        .* (1-inv_logit(Q_p[1:n1,] * theta_p))
      + w[5,1:n1]);
}

model {
  //NNGP
  sigma ~ normal(0, 1);
  phi ~ normal(0, 1);
  for(l in 1:(L-1)) {
    target += sum(-0.5*log(V[l]) - 0.5/sig2[l]*(square(w[l]-uw_dp[l]) ./ V[l])) 
            - 0.5*n3*log(sig2[l]);
  }
  
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
  beta_p ~ normal(0, 1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, diag_pre_multiply(L_sigma[1], L_Omega[1]));
   Y2_ ~ multi_normal_cholesky(nu, diag_pre_multiply(L_sigma[2], L_Omega[2]));
}

generated quantities {
  //landcover: latent compositional
  simplex[L] n_eta[n3];  //gjam transformed nu
  vector[L-1] Y2new_[n3-n1];  //unbiased, split NLCD
  
  Y2new_[,1] = to_array_1d(Y2[n2:n3,1] 
      + w[1,n2:n3]);
  Y2new_[,2] = to_array_1d(Y2[n2:n3,2] 
      + w[1,n2:n3]);
  Y2new_[,3] = to_array_1d(Y2[n2:n3,3] 
      + w[1,n2:n3]);
  Y2new_[,4] = to_array_1d((Y2[n2:n3,4]) 
        .* inv_logit(Q_p[n2:n3,] * theta_p)
      + w[1,n2:n3]);
  Y2new_[,5] = to_array_1d((Y2[n2:n3,4]) 
        .* (1-inv_logit(Q_p[n2:n3,] * theta_p))
      + w[1,n2:n3]);
  
  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(nu[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2new_[n-n1]);
  }
}
