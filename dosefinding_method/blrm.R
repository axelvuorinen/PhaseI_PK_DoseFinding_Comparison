#######################################################################################################################
#### Bayesian Logistic Regression Method (BLRM) - Estimation of the toxicity (Dose-exposure-toxicity relationship) #### 
#######################################################################################################################

# Core function
blrm_model <- function(num, doses, ref_dose, x, y, targeted_tox,  proba_safety = 0.9, options, seed) {
  
  ### Setting up the appropriate parameters to instantiate the Bayesian Logistic Regression Method (BLRM)
  model_BLRM <- stan_model("stan/blrm.stan", "BLRM model")
  alpha_prior <- c(f_logit(targeted_tox), 2)
  beta_prior <- c(0, 1)
  
  ### Model estimation using Bayesian inference
  
  data_stan <- list(N = num, y = y, doses = doses[x], ref_dose = ref_dose, alpha_prior = alpha_prior, beta_prior = beta_prior)
  reg <- sampling(model_BLRM, data = data_stan, iter = options$n_iter, chains = options$n_chains, control = list(adapt_delta = options$n_adapt), cores = options$n_cores)
  a <- get_posterior_mean(reg)
  sampl <- rstan::extract(reg)
  params_output <- a[1:2, options$n_chains + 1]
  p_estim_dose_mean <- f_inv_logit(params_output[1] + exp(params_output[2])*log(doses/ref_dose))
  p_estim_sum <- matrix(sapply(1:length(doses), 
                               function(l) f_inv_logit(sampl$log_alpha + exp(sampl$log_beta)*log(doses[l]/ref_dose))), 
                        nrow = options$n_chains*options$n_iter/2, 
                        ncol = length(doses))
  
  ### Check the safety rule(s) to see if the trial need to be stopped before allocating the recommended dose to the next cohort
  proba_stop <- checking_stop_rule(p_estim_sum[,1], target = targeted_tox, error = 0)
  stop_tox <- (proba_stop >= proba_safety)
  stop_trial <- stop_tox
  
  ### Dose recommendation for next cohort or MTD selection if final estimation with all cohorts included in the trial
  # Check if trial has to be stopped or no
  if (stop_trial) {
    new_dose <- NA
    message("The trial stopped based on the stopping rule.\n \n")
  } else { # if no stopping
    new_dose <- order((abs(p_estim_dose_mean - targeted_tox)))[1]
  }
  
  dose_names <- as.vector(sapply(1:length(doses), function (i) paste("Dose", i)))
  ### Compute probability of toxicity 95% credible intervals for each of the dose based on the 2.5- and 97.5- percentiles for the parameters
  credible_interval_per_dose <- matrix(sapply(1:length(doses), function(i) quantile(p_estim_sum[,i], probs = c(0.025, 0.975))), 
                                       nrow = 2,
                                       ncol = length(doses), 
                                       dimnames = list(c("2.5%", "97.5%"), dose_names),
                                       byrow = FALSE)
  
  params <- c(params_output[1], params_output[2])
  names(params) <- c("log_alpha", "log_beta")
  return(list(new_dose = new_dose, p_estim_dose_mean = p_estim_dose_mean, params = params, credible_interval_per_dose = credible_interval_per_dose))
}