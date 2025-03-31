data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  real log_AUC[N];        // log AUC
  real ref_log_AUC;       // (average) log AUC of reference dose 3
  real beta2_mean;        // beta2 distribution mean value
  real beta2_sd;          // beta2 distribution standard deviation value
  real log_beta3_mean;    // log_beta3 distribution mean value
  real log_beta3_sd;      // log_beta3 distribution standard deviation value
}

parameters {
  real beta2;
  real log_beta3;
}

model {
  for (n in 1:N){
  y[n] ~ bernoulli_logit(-beta2 + exp(log_beta3)*(log_AUC[n]-ref_log_AUC));
  }
  beta2 ~ normal(beta2_mean, beta2_sd);
  log_beta3 ~ normal(log_beta3_mean, log_beta3_sd);
}
