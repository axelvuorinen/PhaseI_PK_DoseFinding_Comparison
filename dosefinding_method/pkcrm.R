#############################################################################################
#### TAKEDA (2017) PK-CRM - Estimation of the toxicity (Dose-exposure-toxicity relationship)
#############################################################################################

# Core function
pkcrm_model <- function(num, doses, x, y, auc_s, targeted_tox, proba_safety = 0.9, options, seed) {
  
  ### Setting up the appropriate parameters to instantiate the PK-CRM method
  alpha_prior <- 10
  beta_prior <- 5
  
  ### Model estimation using Bayesian inference
  
  # Performing linear regression model for auc-dose relationship as stated in the model to estimate the parameters (b0, b1)
  reg_auc_dose <- lm(auc_s ~ log(doses[x]))
  estimated_b0_b1 <- as.vector(reg_auc_dose$coefficients)
  # Use the estimation for the parameters of the frequentist linear regression to estimate the AUC for each dose
  log_estimated_auc_dose <- estimated_b0_b1[1] + estimated_b0_b1[2]*log(doses) # Arbitrary solution to solve the problem of negative AUC
  
  # Computation of the standardized adjustement AUC (x_ij) at each dose level using the previously estimated AUC
  sd_adj_auc <- log_estimated_auc_dose - sum(log_estimated_auc_dose)/length(doses)
  sd_adj_auc_ind <- sd_adj_auc[x]
  
  # Computation of the posterior probability of toxicity at each dose level with standardized adjustment AUC using a Bayesian
  # linear regression model
  data_stan <- list(N = num, y = y, sd_adj_auc = sd_adj_auc_ind, alpha_prior = alpha_prior, beta_prior = beta_prior)
  reg <- sampling(model_PK_CRM, data = data_stan, iter = options$n_iter, chains = options$n_chains, control = list(adapt_delta = options$n_adapt), cores = options$n_cores)
  a <- get_posterior_mean(reg)
  sampl <- rstan::extract(reg)
  mean_alpha_post <- a[1, options$n_chains + 1]
  mean_beta_post <- a[2, options$n_chains + 1]
  p_estim_dose_mean <- f_inv_logit(mean_alpha_post + exp(mean_beta_post)*sd_adj_auc)
  p_estim_sum <- matrix(sapply(1:(options$n_chains*options$n_iter/2), 
                               function(l) f_inv_logit(sampl$alpha[l] + exp(sampl$log_beta[l])*sd_adj_auc)), 
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
  
  params <- c(estimated_b0_b1, mean_alpha_post, mean_beta_post)
  names(params) <- c("beta0", "beta1", "alpha", "log_beta")
  return(list(new_dose = new_dose, p_estim_dose_mean = p_estim_dose_mean, params = params, credible_interval_per_dose = credible_interval_per_dose))
}


