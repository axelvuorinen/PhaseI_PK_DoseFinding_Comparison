data {
  int N;                    // previous patients
  int y[N];                 // binary response
  real doses[N];            // standardized adjustment AUC x_ij (person-specific and dose-specific)
  real ref_dose;            // reference dose for standardisation
  vector[2] alpha_prior;    // alpha mean and variance value for weakly informative normal prior
  vector[2] beta_prior;     // beta mean and variance value for weakly informative normal prior
}

parameters {
  real log_alpha; 
  real log_beta;
}

model {
  for (n in 1:N){
  y[n] ~ bernoulli_logit(log_alpha + exp(log_beta)*log(doses[n]/ref_dose));
  }
  log_alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  log_beta ~ normal(beta_prior[1], beta_prior[2]);
}
