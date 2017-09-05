functions {
  //NNGP lpdf based on Lu Zhang's work
  real nngp_w_lpdf(vector w, 
                   real sigmasq, 
                   real phi, 
                   matrix neardist,
                   matrix neardistM,
                   int[,] nearind, 
                   int N,
                   int M,
                   int[] dim) {
    vector[N] V;
    vector[N] Uw;
    int h;
    vector[N] out_;
    real out;
    
    Uw = w;
    V[1] = 1;
    out_[1] = - 0.5 * log(V[1]) - 0.5 / sigmasq * (Uw[1] * Uw[1] / V[1]);
    
    for(i in 2:N) {
      matrix[dim[i-1], dim[i-1]] temp_neardistM;
      matrix[dim[i-1], dim[i-1]] L;
      vector[dim[i-1]] u;
      vector[dim[i-1]] v;
      row_vector[dim[i-1]] v2;
      
      if(dim[i-1] == 1) temp_neardistM[1,1] = 1;
      else {
        h = 0;
        for(j in 1:(dim[i-1]-1)) {
          for(k in (j+1):dim[i-1]) {
            h = h + 1;
            temp_neardistM[j,k] = exp(-phi * neardistM[(i-1),h]);
            temp_neardistM[k,j] = temp_neardistM[j,k];
          } //close k
          temp_neardistM[j,j] = 1;
        } //close j
        temp_neardistM[dim[i-1], dim[i-1]] = 1;
      } //close else
      
      L = cholesky_decompose(temp_neardistM);
      for(j in 1:dim[i-1]) {
        u[j] = exp(-phi * neardist[(i-1),j]);
      } //close j
      v = mdivide_left_tri_low(L, u);
      V[i] = 1 - (v' * v);
      v2 = mdivide_right_tri_low(v', L);
      for(j in 1:dim[i-1]) {
        Uw[i] = Uw[i] - v2[j] * w[nearind[(i-1), j]];
      } //close j
      
      out_[i] = - 0.5 * log(V[i]) - 0.5 / sigmasq * (Uw[i] * Uw[i] / V[i]);
    } //close i
    out = sum(out_) - 0.5 * N * log(sigmasq);
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
  int n1;  //number of cells with Y1 & Y2
  int n2;  //n1 + 1 (for indexing)
  int n3;  //number of cells with only Y2
  int L;  //number of land cover classes
  int nB_d[L-2];  //number of bias covariates for each LC
  int nB_p;  //number of covariates for pr(WP|Evg)
  //landcover: observed
  vector<lower=0, upper=1>[L-1] Y1[n1];  //GRANIT proportions
  matrix<lower=0, upper=1>[n3,L-2] Y2;  //NLCD proportions
  //covariates
  matrix[n3,nB_p] X_p;  //pr(WP|Evg) covariates
  matrix[n3,nB_d[1]] X_d1;  //bias covariates: Dev
  matrix[n3,nB_d[2]] X_d2;  //bias covariates: Oth
  matrix[n3,nB_d[3]] X_d3;  //bias covariates: Hwd
  matrix[n3,nB_d[4]] X_d4;  //bias covariates: Evg
  //NNGP spatial random effects
  int<lower=1, upper=n3> M;  //number of neighbors
  int nearind[n3-1, M];
  matrix[n3-1, M] neardist;
  matrix[n3-1, (M*(M-1)/2)] neardistM;
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
  //NNGP matrix sizes
  int dim[n3-1];
  for(i in 2:n3) dim[i-1] = i < (M+1) ? (i-1) : M;
}

parameters {
  //landcover: covariance
  cholesky_factor_corr[L-1] L_Omega[2];
  vector<lower=0>[L-1] L_sigma[2];
  //landcover: latent not constrained to be compositional
  vector[L-1] nu[n1];
  //thetas: QR decomposition slopes
  vector[nB_p] theta_p;  //pr(WP|Evg) thetas
  vector[n_beta_d] theta_d;  //bias thetas
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
  vector[n_beta_d] beta_d;  //bias betas

  
  //QR decopmositions
  beta_p = R_inv_p * theta_p;
  beta_d[1:d1_2] = R_inv_d1 * theta_d[1:d1_2];
  beta_d[d2_1:d2_2] = R_inv_d2 * theta_d[d2_1:d2_2];
  beta_d[d3_1:d3_2] = R_inv_d3 * theta_d[d3_1:d3_2];
  beta_d[d4_1:d4_2] = R_inv_d4 * theta_d[d4_1:d4_2];
  
  
  //correct bias & split WP to [,4] and Evg to [,5]
  ////fit betas using cells with Y1 & Y2
  Y2_[1:n1,1] = to_array_1d(Y2[1:n1,1] 
      + (Q_d1[1:n1,] * theta_d[1:d1_2]) 
      + w[1,1:n1]);
  Y2_[1:n1,2] = to_array_1d(Y2[1:n1,2] 
      + (Q_d2[1:n1,] * theta_d[d2_1:d2_2])
      + w[2,1:n1]);
  Y2_[1:n1,3] = to_array_1d(Y2[1:n1,3] 
      + (Q_d3[1:n1,] * theta_d[d3_1:d3_2])
      + w[3,1:n1]);
  Y2_[1:n1,4] = to_array_1d((Y2[1:n1,4] 
          + (Q_d4[1:n1,] * theta_d[d4_1:d4_2])) 
        .* inv_logit(Q_p[1:n1,] * theta_p)
      + w[4,1:n1]);
  Y2_[1:n1,5] = to_array_1d((Y2[1:n1,4] 
          + (Q_d4[1:n1,] * theta_d[d4_1:d4_2])) 
        .* (1-inv_logit(Q_p[1:n1,] * theta_p))
      + w[5,1:n1]);
}

model {
  //NNGP
  sigma ~ normal(0, 1);
  phi ~ normal(0, 1);
  for(l in 1:(L-1)) {
    w[l] ~ nngp_w(sigma[l]^2, phi[l], neardist, neardistM, nearind, n3, M, dim);
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
  beta_d ~ normal(0, 0.1);
  
  //likelihood
   Y1 ~ multi_normal_cholesky(nu, diag_pre_multiply(L_sigma[1], L_Omega[1]));
   Y2_ ~ multi_normal_cholesky(nu, diag_pre_multiply(L_sigma[2], L_Omega[2]));
}

generated quantities {
  //landcover: latent compositional
  simplex[L] n_eta[n3];  //gjam transformed nu
  vector[L-1] Y2new_[n3-n1];  //unbiased, split NLCD
  
  Y2new_[,1] = to_array_1d(Y2[n2:n3,1] 
      + (Q_d1[n2:n3,] * theta_d[1:d1_2])
      + w[1,n2:n3]);
  Y2new_[,2] = to_array_1d(Y2[n2:n3,2] 
      + (Q_d2[n2:n3,] * theta_d[d2_1:d2_2])
      + w[1,n2:n3]);
  Y2new_[,3] = to_array_1d(Y2[n2:n3,3] 
      + (Q_d3[n2:n3,] * theta_d[d3_1:d3_2])
      + w[1,n2:n3]);
  Y2new_[,4] = to_array_1d((Y2[n2:n3,4] 
          + (Q_d4[n2:n3,] * theta_d[d4_1:d4_2])) 
        .* inv_logit(Q_p[n2:n3,] * theta_p)
      + w[1,n2:n3]);
  Y2new_[,5] = to_array_1d((Y2[n2:n3,4] 
          + (X_d4[n2:n3,] * theta_d[d4_1:d4_2])) 
        .* (1-inv_logit(Q_p[n2:n3,] * theta_p))
      + w[1,n2:n3]);
  
  for(n in 1:n1) {
    n_eta[n] = tr_gjam_inv(nu[n]);
  }
  for(n in n2:n3) {
    n_eta[n] = tr_gjam_inv(Y2new_[n-n1]);
  }
}
