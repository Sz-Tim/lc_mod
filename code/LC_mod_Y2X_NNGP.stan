functions {
  //NNGP lpdf based on Lu Zhang's work
  real nngp_w_lpdf(vector w, 
                   real sigmasq, 
                   real phi, 
                   matrix neardist,
                   matrix neardistM,
                   int[,] nearind, 
                   int N,
                   int M) {
    vector[N] V;
    vector[N] Uw;
    int dim;
    int h;
    real out;
    
    // removed this section, remove beta1, just set Uw = w_r
    Uw = w;
    
    for(i in 2:N) {
      matrix[i < (M+1) ? (i-1) : M, i < (M+1) ? (i-1) : M] temp_neardistM;
      matrix[i < (M+1) ? (i-1) : M, i < (M+1) ? (i-1) : M] L;
      vector[i < (M+1) ? (i-1) : M] u;
      vector[i < (M+1) ? (i-1) : M] v;
      row_vector[i < (M+1) ? (i-1) : M] v2;
      
      dim = (i < (M+1)) ? (i-1) : M;
      
      //get exp(-phi * neardistM)
      if(dim == 1) temp_neardistM[1,1] = 1;
      else {
        h = 0;
        for(j in 1:(dim-1)) {
          for(k in (j+1):dim) {
            h = h + 1;
            temp_neardistM[j,k] = exp(-phi * neardistM[(i-1),h]);
            temp_neardistM[k,j] = temp_neardistM[j,k];
          }
          temp_neardistM[j,j] = 1;
        }
        temp_neardistM[dim, dim] = 1;
      }
      
      L = cholesky_decompose(temp_neardistM);
      
      for(j in 1:dim) {
        u[j] = exp(-phi * neardist[(i-1),j]);
      }
      
      v = mdivide_left_tri_low(L, u);
      V[i] = 1 - (v' * v);
      
      v2 = mdivide_right_tri_low(v', L);
      
      for(j in 1:dim) {
        Uw[i] = Uw[i] - v2[j] * w[nearind[(i-1), j]];
      }
    }
    
    V[1] = 1;
    out = 0.0;
    for(i in 1:N) {
      out = out - 0.5 * log(V[i]) - 0.5 / sigmasq * (Uw[i] * Uw[i] / V[i]);
    }
    out = out - 0.5 * N * log(sigmasq);
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
  int nearind[n3-1, M];
  matrix[n3-1, M] neardist;
  matrix[n3-1, (M*(M-1)/2)] neardistM;
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
  //gpp
  sigma ~ normal(0, 1);
  phi ~ normal(0, 1);
  for(l in 1:(L-1)) {
    w[l] ~ nngp_w(sigma[l]^2, phi[l], neardist, neardistM, nearind, n3, M);
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
