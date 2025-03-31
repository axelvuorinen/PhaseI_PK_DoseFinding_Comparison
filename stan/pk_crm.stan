data {
  int N;                    // previous patients
  int y[N];                 // binary response
  real sd_adj_auc[N];       // standardized adjustment AUC x_ij (person-specific and dose-specific)
  real alpha_prior;         // alpha variance value for weakly informative normal prior
  real beta_prior;          // beta variance value for weakly informative normal prior
}

parameters {
  real alpha; 
  real log_beta;
}

model {
  for (n in 1:N){
  y[n] ~ bernoulli_logit(alpha + exp(log_beta)*sd_adj_auc[n]);
  }
  alpha ~ normal(0, alpha_prior);
  log_beta ~ normal(0, beta_prior);
}
