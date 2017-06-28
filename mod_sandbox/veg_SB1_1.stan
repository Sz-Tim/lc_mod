data {
  int N; //number of grid cells
  vector[N] Y1;
  vector[N] Y2;
}
parameters {
  vector<lower=0, upper=1>[N] nu;
  real<lower=0> sig;
}
model {
  //priors
  nu ~ beta(0.6, 0.75);
  sig ~ uniform(0, 2);
  
  //likelihood
  Y1 ~ normal(nu, sig);
  Y2 ~ normal(nu, sig);
}



