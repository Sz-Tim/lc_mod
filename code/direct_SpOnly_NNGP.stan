functions {
  //NNGP lpdf based on Lu Zhang's work
  real nngp_w_lpdf(vector w, 
                   real sigmasq, 
                   real phi, 
                   matrix nn_d,
                   matrix nn_dM,
                   matrix[] nn_dM_i,
                   int[,] nn_id, 
                   int N,
                   int M,
                   int[] nn_dim) {
    vector[N] V;
    vector[N] Uw;
    vector[N] out_;
    real out;
    
    for(i in 1:N) {
      matrix[nn_dim[i], nn_dim[i]] exp_nn_dM;
      matrix[nn_dim[i], nn_dim[i]] L;
      vector[nn_dim[i]] v;
      row_vector[nn_dim[i]] v_t;
      
      exp_nn_dM = exp(-phi * block(nn_dM_i[i], 1, 1, nn_dim[i], nn_dim[i]));
      for(j in 1:nn_dim[i]) exp_nn_dM[j,j] = 1;
      
      L = cholesky_decompose(exp_nn_dM);
      v = mdivide_left_tri_low(L, exp(-phi * sub_col(nn_d, 1, i, nn_dim[i])));
      v_t = v';
      V[i] = 1 - (v_t * v);
      Uw[i] = w[i] - dot_product(mdivide_right_tri_low(v_t, L), 
                                     w[nn_id[i, 1:nn_dim[i]]]);
    }
    out_ = -0.5*log(V) - 0.5/sigmasq*(square(Uw) ./ V);
    out = sum(out_) - 0.5*N*log(sigmasq);
    return out;
  }
  
  
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
  int nn_id[n3, M];  //neighbor indices
  matrix[M, n3] nn_d;  //neighbor distances
  matrix[n3, (M*(M-1)/2)] nn_dM;  //neighbor distance matrix
  int nn_dim[n3];  //number of neighbors for each cell
}

transformed data {
  //QR decomposition for covariates
  matrix[n3,nB_p] Q_p = qr_Q(X_p)[,1:nB_p] * sqrt(n3-1);
  matrix[nB_p,nB_p] R_p = qr_R(X_p)[1:nB_p,] / sqrt(n3-1);
  matrix[nB_p,nB_p] R_inv_p = inverse(R_p);
  //NNGP sub-matrix distances
  matrix[M,M] nn_dM_i[n3];
  for(i in 1:n3) {
    int h = 0;
    nn_dM_i[i] = diag_matrix(rep_vector(1, M));
    for(j in 1:(nn_dim[i]-1)) {
      for(k in (j+1):nn_dim[i]) {
        h = h + 1;
        nn_dM_i[i,j,k] = nn_dM[i,h];
        nn_dM_i[i,k,j] = nn_dM[i,h];
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
    w[l] ~ nngp_w(sigma[l]^2, phi[l], nn_d, nn_dM, nn_dM_i,
                  nn_id, n3, M, nn_dim);
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
