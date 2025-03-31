data {
  int<lower=0> N;       // number of previous patients
  vector[N] log_AUC;    // log AUC
  matrix[N,2] doses;    // log doses + intercept
  real beta0_mean;      // beta0 distribution mean
  real beta01_sd;       // beta0 and beta1 distributions standard deviation
}

parameters {
  vector[2] b; 
  real<lower=0, upper=1> nu;
}

model {
  log_AUC ~ normal(doses*b, nu);
  nu ~ beta(1,1);
  b[1] ~ normal(beta0_mean, beta01_sd);
  b[2] ~ normal(1, beta01_sd);
}
