##############################################################################################
#### Ursino (2017) PKLOGIT - Estimation of the toxicity (Dose-exposure-toxicity relationship)
##############################################################################################

# Core function
pklogit_model <- function(num, doses, x, y, auc_s, targeted_tox, proba_safety = 0.9, options, seed) {
  
  ### Setting up the appropriate parameters to instantiate the PKLOGIT method
  beta0_mean <- -log(Cl_pop)
  # "Low" uninformative standard deviation for beta0 and beta1 normal distributions
  beta01_sd <- 10
  # "High" uninformative standard deviation for beta0 and beta1 normal distributions - Original method
  #beta01_sd <- sqrt(1000)
  # "Super high" uninformative standard deviation for beta0 and beta1 normal distributions - Original implementation
  #beta01_sd <- 10000
  beta2_mean <- 6
  beta2_sd <- 5
  log_beta3_mean <- 2
  log_beta3_sd <- 1
  
  ### Model estimation using Bayesian inference
  
  # Stan for AUC computation using normal prior distribution
  doses_AUC <- cbind(rep(1, num), log(doses[x]))
  data_stan <- list(N = num, log_AUC = log(auc_s), doses = doses_AUC, beta0_mean = beta0_mean, 
                    beta01_sd = beta01_sd)
  reg_AUC <- sampling(model_PKLOGIT_AUC,
                      data = data_stan, iter = options$n_iter, chains = options$n_chains,
                      control = list(adapt_delta = options$n_adapt),
                      cores = options$n_cores
  )
  res_AUC <- get_posterior_mean(reg_AUC)
  sampl_AUC <- rstan::extract(reg_AUC)
  beta01 <- res_AUC[1:2, options$n_chains + 1]
  nu <- res_AUC[3, options$n_chains + 1]
  
  # Simulation of data using normal distribution of log(AUC) to generate Monte Carlo sample to compute the integral of the product of the dose-exposure density function and the
  # posterior probability of toxicity
  n_simu <- 10000
  simu_sample_logAUC <- matrix(sapply(1:length(doses), function(i) rnorm(n_simu, res_AUC[1, options$n_chains + 1] + res_AUC[2, options$n_chains + 1]*log(doses[i]), nu)),
                               nrow = n_simu,
                               ncol = length(doses),
                               byrow = FALSE)
  simu_ref_log_AUC <- mean(simu_sample_logAUC[, which(doses == ref_dose)])
  
  # Stan for DLT rate computation using uniform priors and a logit link function
  data_stan <- list(N = num, y = y, log_AUC = log(auc_s), ref_log_AUC = simu_ref_log_AUC, 
                    beta2_mean = beta2_mean, beta2_sd = beta2_sd, 
                    log_beta3_mean = log_beta3_mean, log_beta3_sd = log_beta3_sd)
  reg_logit <- sampling(model_PKLOGIT_logit,
                        data = data_stan, iter = options$n_iter, chains = options$n_chains,
                        control = list(adapt_delta = options$n_adapt),
                        cores = options$n_cores
  )
  res_logit <- get_posterior_mean(reg_logit)
  sampl_logit <- rstan::extract(reg_logit)
  beta2 <- res_logit[1, options$n_chains + 1]
  log_beta3 <- res_logit[2, options$n_chains + 1]
  
  # Computation of predictive probability of toxicity to determine the recommended dose to be given to the next cohort OR the MTD using the
  # final posterior probabilities of toxicity for the trial if all cohorts have been included
  p_estim_dose_tox_MCMC <- matrix(sapply(1:length(doses), function(i) f_inv_logit(-res_logit[1, options$n_chains + 1] + exp(res_logit[2, options$n_chains + 1])*(simu_sample_logAUC[,i] - simu_ref_log_AUC))),
                               nrow = n_simu,
                               ncol = length(doses),
                               byrow = FALSE)
  p_estim_dose_mean <- sapply(1:length(doses), function(i) mean(p_estim_dose_tox_MCMC[, i]))
  
  ### The posterior probabilities of toxicity
  centered_sample_logAUC <- matrix(sapply(simu_sample_logAUC, function(i) i - simu_ref_log_AUC), nrow = n_simu, ncol = length(doses), byrow = FALSE)
  beta2_mc <- sampl_logit$beta2
  log_beta3_mc <- sampl_logit$log_beta3
  p_estim_sum <- matrix(sapply(1:((options$n_chains*options$n_iter)/2), function(o) {
    p_estim_sample_logAUC <- f_inv_logit(-beta2_mc[o] + exp(log_beta3_mc[o])*centered_sample_logAUC)
    return(colMeans(p_estim_sample_logAUC))
  }), nrow = options$n_chains*options$n_iter/2, ncol = length(doses), byrow = TRUE)
  
  ### Check the safety rule(s) to see if the trial need to be stopped before allocating the recommended dose to the next cohort
  proba_stop <- checking_stop_rule(p_estim_sum[, 1], target = targeted_tox, error = 0)
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
  
  
  params <- c(beta01, nu, beta2, log_beta3)
  names(params) <- c("beta0", "beta1", "nu", "beta2", "log_beta3")
  return(list(new_dose = new_dose, p_estim_dose_mean = p_estim_dose_mean, p_estim_sum = p_estim_sum, params = params, credible_interval_per_dose = credible_interval_per_dose, ref_log_AUC = simu_ref_log_AUC))
}