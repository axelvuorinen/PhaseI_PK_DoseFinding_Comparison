##############################################################################################
#### Günhan (2021) TITE-PK - Estimation of the toxicity (Dose-exposure-toxicity relationship)
##############################################################################################

# Core function
titepk_model <- function(num, doses, stan_data_titepk, exposure_tox_threshold, targeted_tox, proba_safety,
                         is_overdose_control = FALSE, options, seed) {
  
  ### Setting up the appropriate parameters to instantiate the PKLOGIT method
  # Stan code for estimation
  # TITE-PK: Prior value for log(beta) parameters c(cloglog(0.25), 1.25)
  # Thresholds for overdose control criterion
  targeted_tox_LB <- 0.15 # Lower bound of targeted toxicity interval
  overdosing_threshold_prob <- 0.35
  overdosing_risk_prob <- 0.25
  
  ### Model estimation using Bayesian inference
  # Calculate prior DLT probabilities
  fit_posteriors_daily_tte <- sampling(model_TITE_PK, 
                                       data = stan_data_titepk,
                                       init = "random",
                                       chains = options$n_chains,
                                       iter = options$n_iter,
                                       control = list(adapt_delta = options$n_adapt),
                                       cores = options$n_cores,
                                       open_progress = FALSE)

  # Extract the necessary values
  sampl <- rstan::extract(fit_posteriors_daily_tte)
  vec_name_doses_posterior <- sapply(1:length(doses), function(i) paste0("P_dose[", i, "]"))
  prefs_posteriors_daily_tte <- matrix(summary(fit_posteriors_daily_tte)$summary[vec_name_doses_posterior, 
                                                                                 c(1, 4, 5, 6, 7, 8)], ncol = 6)
  rownames(prefs_posteriors_daily_tte) <- vec_name_doses_posterior
  colnames(prefs_posteriors_daily_tte) <- c("mean", "2.5%", "25%", "50%", "75%", "97.5%")
  p_estim_dose_mean <- prefs_posteriors_daily_tte[,1]
  log_beta <- summary(fit_posteriors_daily_tte)$summary["logbeta",1]
  p_estim_sum <- sampl$P_dose
  
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
    if (is_overdose_control) {
      # Dose recommendation rule with overdose control to have an "EWOC"-like criterion for MTD identification
      proba_overdose_ctrl <- sapply(seq_along(doses), function(l) mean(p_estim_sum[,l] > overdosing_threshold))
      proba_targeted_tox <- sapply(seq_along(doses), function(l) mean(p_estim_sum[,l] > targeted_tox_LB & p_estim_sum[,l] <= overdosing_threshold))
      feasible_doses <- which(proba_overdose_ctrl < overdosing_tolerance_prob)
      new_dose <- which(proba_targeted_tox[feasible_doses] == max(proba_targeted_tox[feasible_doses])) 
    } else {
      # Classical "CRM"-like dose recommendation rule
      new_dose <- order((abs(p_estim_dose_mean - targeted_tox)))[1]
    }
  }
  
  dose_names <- as.vector(sapply(1:length(doses), function (i) paste("Dose", i)))
  ### Compute probability of toxicity 95% credible intervals for each of the dose based on the 2.5- and 97.5- percentiles for the parameters
  credible_interval_per_dose <- matrix(sapply(1:length(doses), function(i) quantile(p_estim_sum[,i], probs = c(0.025, 0.975))), 
                                       nrow = 2,
                                       ncol = length(doses), 
                                       dimnames = list(c("2.5%", "97.5%"), dose_names),
                                       byrow = FALSE)
  
  params <- log_beta
  names(params) <- "log_beta"
  return(list(new_dose = new_dose, p_estim_dose_mean = p_estim_dose_mean, params = params, credible_interval_per_dose = credible_interval_per_dose))
}