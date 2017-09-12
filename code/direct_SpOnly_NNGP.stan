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
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  vector<lower=0, upper=1>[L-2] Y2[n3];  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  //NNGP spatial random effects
  int<lower=1, upper=n3> M;  //number of neighbors
  int nn_id[M, n3];  //neighbor indices
  matrix[M, n3] nn_d;  //neighbor distances
  matrix[(M*(M-1)/2), n3] nn_dM;  //neighbor distance matrix
  int nn_dim[n3];  //number of neighbors for each cell
  int dim_r;  //number of unique neighborhood sizes
  int dim_i[dim_r, 3];  //reference for neighborhood size index start & end
}

transformed data {
  //QR decomposition for covariates
  matrix[n3,nB_p] Q_p = qr_Q(X_p)[,1:nB_p] * sqrt(n3-1);
  matrix[nB_p,nB_p] R_p = qr_R(X_p)[1:nB_p,] / sqrt(n3-1);
  matrix[nB_p,nB_p] R_inv_p = inverse(R_p);
  //NNGP neighborhood size indices
  int M_i[dim_r] = dim_i[,1];  //neighborhood sizes: unique(nn_dim)
  int M_1[dim_r] = dim_i[,2];  //neighborhood sizes: index starts
  int M_2[dim_r] = dim_i[,3];  //neighborhood sizes: index ends
  int n_i[dim_r];  //number of cells w/neighborhood size i
  matrix[M,M] nn_dM_i[n3];
  for(r in 1:dim_r) n_i[r] = M_2[r] - M_1[r] + 1;
  //NNGP sub-matrix distances
  for(i in 1:n3) {
    int h = 0;
    nn_dM_i[i] = diag_matrix(rep_vector(1, M));
    for(j in 1:(nn_dim[i]-1)) {
      for(k in (j+1):nn_dim[i]) {
        h = h + 1;
        nn_dM_i[i,j,k] = nn_dM[h,i];
        nn_dM_i[i,k,j] = nn_dM[h,i];
      } //close k
    } //close j
  }
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega; 
  vector<lower=0>[L-1] L_sigma; 
  //betas
  vector[nB_p] theta_p;  //pr(WP|Evg) betas (QR decomposition)
  //NNGP
  real<lower=0> sigma[L-1];  //sqrt(nugget)
  real<lower=0> phi[L-1];  //decay rate
  vector[n3] w[L-1];  //spatial effects
}

transformed parameters {
  //NLCD de-biasing and splitting
  vector[L-1] Y2_[n1];
  vector[nB_p] beta_p;
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
      matrix[M_i[r],M_i[r]] exp_nn_dM[n_i[r]];
      matrix[M_i[r],M_i[r]] L_nn[n_i[r]];
      matrix[n_i[r],M_i[r]] L_u;
      matrix[n_i[r],M_i[r]] v_L;
      
      for(i in 1:n_i[r]) {
        int i_g = i + M_1[r] - 1;  //global index
        exp_nn_dM[i] = exp(-phi[l] * block(nn_dM_i[i_g],1,1,M_i[r],M_i[r]));
        for(j in 1:M_i[r]) exp_nn_dM[i,j,j] = 1;
        L_nn[i] = cholesky_decompose(exp_nn_dM[i]);
        L_u[i,] = mdivide_left_tri_low(L_nn[i], sub_col(um[l],1,i_g,M_i[r]))';
        v_L[i,] = mdivide_right_tri_low(L_u[i,], L_nn[i]);
      }
      V[l,M_1[r]:M_2[r]] = 1 - rows_dot_self(L_u);
      uw_dp[l,M_1[r]:M_2[r]] = rows_dot_product(v_L, 
          to_matrix(w[l, to_array_1d(nn_id[1:M_i[r], M_1[r]:M_2[r]])], 
                    n_i[r], M_i[r]));
    }
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
  //NNGP
  sigma ~ normal(0, 1);
  phi ~ normal(0, 1);
  for(l in 1:(L-1)) {
    target += sum(-0.5*log(V[l]) - 0.5/sig2[l]*(square(w[l]-uw_dp[l]) ./ V[l])) 
            - 0.5*n3*log(sig2[l]);
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
