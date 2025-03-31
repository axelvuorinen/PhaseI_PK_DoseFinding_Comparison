################################################################################################
#### Micallef (2018) ED-EWOC - Estimation of the toxicity (Dose-exposure-toxicity relationship)
################################################################################################

# Core function
ed_ewoc_model <- function(num, doses, ref_dose, x, y, conci, auc_s, n_sampling_timepoints, time_sampling,
                          real_sampling, PK_parameters, targeted_tox, proba_safety = 0.9, options, seed) {
  
  ### Setting up the appropriate parameters to instantiate the ED-EWOC method
  # Initilization phase with BLRM 
  alpha_prior <- c(f_logit(targeted_tox), 0.4)
  beta_prior <- c(0.6, 1)
  # Main model using ED-EWOC specific methodology
  beta2_mean <- -6
  beta2_sd <- 5
  log_beta3_mean <- 2
  log_beta3_sd <- 1
  # EWOC bundaries
  targeted_tox_LB <- 0.15
  overdosing_threshold_prob <- 0.35
  overdosing_risk_prob <- 0.25
  
  ### Model estimation using Bayesian inference
  if(num <= N/2) {
    
    ### Initialization phase for ED-EWOC model to conduct the first escalation part of the trial without using the proper model
    ### to wait for more patients inclusion before starting the main phase of the model where the correct methodology will be applied
    ### to estimate a popPK model and the toxicity probabilities
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
    stop_trial <- FALSE
    ### Stop the trial if the posterior probability of overdosing exceeds 90%
    proba_stop <- checking_stop_rule(p_estim_sum[,1], target = targeted_tox, error = 0)
    stop_tox <- (proba_stop >= proba_safety)
    if (stop_tox == TRUE) {
      stop_trial <- stop_tox
    }
    
    ### Dose recommendation for next cohort or MTD selection if final estimation with all cohorts included in the trial
    # Check if trial has to be stopped or no
    if (stop_trial) {
      new_dose <- NA
      message("The trial stopped based on the stopping rule or one of the stopping rule(s).\n \n")
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
    
    # Display for rstan estimated parameters
    params <- c(params_output[1], params_output[2])
    names(params) <- c("log_alpha", "log_beta")
    estimated_PK_parameters <- NULL
    simu_ref_log_AUC <- NULL
  }
  else {
    ### Running the main design for the ED-EWOC method once half of the patients have been included using an initialization phase with a BLRM
    ka <- PK_parameters[1]
    Cl_pop <- PK_parameters[2]
    V_pop <- PK_parameters[3]
    # Updating population PK model based on measured concentrations (simulated concentrations with error)
    id <- c()
    seq_doses <- c()
    time <- c()
    for(l in 1:num) {
      id <- c(id, rep(l, n_sampling_timepoints))
      seq_doses <- c(seq_doses, rep(doses[x[l]], n_sampling_timepoints))
      time <- c(time, time_sampling[real_sampling])
      #time <- c(time, real_sampling)
    }
    input_pkdata <- data.frame(id = id, dose = seq_doses, time = time, concentration = conci)
    saemix.data <- saemixData(name.data = input_pkdata, header = TRUE, sep = " ", na = NA, 
                              name.group = c("id"), name.predictors = c("dose", "time"), 
                              name.response = c("concentration"), name.X = "time", units = list(x = "hr", y = "mg/L"))
    model_1cpt <- function(psi, id, xidep) {
      dose <- xidep[,1]
      tim <- xidep[,2]
      ka <- psi[id,1]
      V <- psi[id,2]
      Cl <- psi[id,3]
      k <- Cl/V
      y_pred <- dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
      return(y_pred)
    }
    random_lowU_initialization <- rnorm(3, 0, 0.2)
    saemix.model <- saemixModel(model = model_1cpt, description = "One-compartment model with first-order absorption",
                                psi0 = matrix(c(ka*(1+random_lowU_initialization[1]), V_pop*(1+random_lowU_initialization[2]), Cl_pop*(1+random_lowU_initialization[3])), ncol = 3, byrow = TRUE, dimnames = list(NULL, c("ka", "V", "Cl"))),
                                covariance.model = matrix(c(0,0,0,0,1,0,0,0,1), ncol = 3, byrow = TRUE),
                                omega.init = matrix(c(0,0,0,0,1,0,0,0,1), ncol = 3, byrow = TRUE),
                                error.model = "proportional",
                                transform.par = c(0,1,1), fixed.estim = c(1,1,1),
                                verbose = FALSE)
    saemix.options <- list(seed = seed, nb.chains = 4, nbiter.saemix = c(400, 200), save = FALSE, save.graphs = FALSE, warnings = TRUE)
    saemix.fit <- saemix(saemix.model, saemix.data, saemix.options)
    
    # Computation of new AUCs using updated model to predict concentrations for each dose
    updated_ka <- saemix.fit@results@fixed.effects[1]
    updated_Cl_pop <- saemix.fit@results@fixed.effects[3]
    updated_V_pop <- saemix.fit@results@fixed.effects[2]
    updated_omega_Cl <- saemix.fit@results@omega[3,3]
    updated_omega_V <- saemix.fit@results@omega[2,2]
    updated_Cl_ind <- saemix.fit@results@map.psi[,3]
    updated_V_ind <- saemix.fit@results@map.psi[,2]
    updated_k <- updated_Cl_ind / updated_V_ind
    predicted_auc_updated <- sapply(1:num, function(l) {
      integrand <- function(t) {
        (doses[x[l]] / updated_V_ind[l]) * (updated_ka / (updated_ka - updated_k[l])) * (exp(-updated_k[l] * t) - exp(-updated_ka * t))
      }
      return(integrate(integrand, lower = 0, upper = Inf)$value)
    })
    auc_s <- predicted_auc_updated
    
    # Simulation of data to generate Monte Carlo sample to compute the integral of the product of the dose-exposure density function and the
    # posterior probability of toxicity
    n_simu <- 10000
    eta_Cl_simu <- exp(rnorm(n_simu, 0, updated_omega_Cl))
    eta_V_simu <- exp(rnorm(n_simu, 0, updated_omega_V))
    Cl_ind_simu <- updated_Cl_pop*eta_Cl_simu
    V_ind_simu <- updated_V_pop*eta_V_simu
    # For scenarios in set E and H to avoid triggering critical integration failures when using integrate function and obtain 
    # stable values for the predicted AUCs from the newly estimated pop-PK model, we use hcubature function
    #simu_sample <- AUC_estim_popPK(n_simu, doses, updated_ka, Cl_ind_simu, V_ind_simu)
    simu_sample <- AUC_estim_popPK_stable(n_simu, doses, updated_ka, Cl_ind_simu, V_ind_simu)
    simu_ref_log_AUC <- log(mean(simu_sample[, which(doses == ref_dose)]))
    
    # Stan for posterior parameters distribution estimation
    data_stan <- list(N = num, y = y, log_AUC = log(auc_s), ref_log_AUC = simu_ref_log_AUC,
                      beta2_mean = beta2_mean, beta2_sd = beta2_sd, 
                      log_beta3_mean = log_beta3_mean, log_beta3_sd = log_beta3_sd)
    reg <- sampling(model_ED_EWOC_ED, data = data_stan, iter = options$n_iter, 
                    chains = options$n_chains, control = list(adapt_delta = options$n_adapt), 
                    cores = options$n_cores)
    a <- get_posterior_mean(reg)
    sampl <- rstan::extract(reg)
    
    # Computation of the posterior mean rates of toxicity using MCMC samples for each AUC
    beta2 <- a[1, options$n_chains + 1]
    log_beta3 <- a[2, options$n_chains + 1]
    
    log_centered_AUC_MC_simu <- matrix(sapply(simu_sample, function(i) (log(i) - simu_ref_log_AUC)), nrow = n_simu, ncol = length(doses), byrow = FALSE)
    p_estim_auc_simu <- f_inv_logit(beta2 + exp(log_beta3)*log_centered_AUC_MC_simu)
  
    # Computation of the posterior estimated posterior probabilities of toxicities using MCMC samples for each dose
    p_estim_dose_mean <- sapply(1:length(doses), function(l) mean(p_estim_auc_simu[,l]))
    beta2_mc <- sampl$beta2
    log_beta3_mc <- sampl$log_beta3
    f_estim_sum <- function(o) {
      p_estim_sum_auc <- f_inv_logit(beta2_mc[o] + exp(log_beta3_mc[o])*log_centered_AUC_MC_simu)
      return(colMeans(p_estim_sum_auc))
    }
    p_estim_sum <- matrix(sapply(1:((options$n_chains*options$n_iter)/2), f_estim_sum), nrow = options$n_chains*options$n_iter/2, ncol = length(doses), byrow = TRUE)
    
    ### Computation of the overdose control probability and targeted toxicity probability to have an "EWOC"-like decision rule for dose recommendation
    candidate_doses <- c()
    proba_overdose_ctrl <- rep(NA, length(doses))
    proba_targeted_tox <- rep(NA, length(doses))
    for (l in 1:length(doses)) {
      if (all(is.na(p_estim_sum[,l])) == FALSE) {
        proba_overdose_ctrl[l] <- sum(p_estim_sum[which(!is.na(p_estim_sum[,l])),l] > overdosing_threshold_prob) / length(p_estim_sum[which(!is.na(p_estim_sum[,l])),l] )
        if (proba_overdose_ctrl[l] < overdosing_risk_prob) {
          candidate_doses <- c(candidate_doses,l)
          proba_targeted_tox[l] <- sum(p_estim_sum[which(!is.na(p_estim_sum[,l])),l] > targeted_tox_LB & p_estim_sum[which(!is.na(p_estim_sum[,l])),l] <= overdosing_threshold_prob) / length(p_estim_sum[which(!is.na(p_estim_sum[,l])),l])
        }
      }
    }
    
    ### Check the safety rule(s) to see if the trial need to be stopped before allocating the recommended dose to the next cohort
    stop_trial <- FALSE
    ### Stop the trial if the posterior probability of overdosing exceeds 90%
    proba_stop <- checking_stop_rule(p_estim_sum[,1], target = targeted_tox, error = 0)
    stop_tox <- (proba_stop >= proba_safety)
    if (stop_tox) {
      stop_trial <- stop_tox
    }
    
    ### Dose recommendation for next cohort or MTD selection if final estimation with all cohorts included in the trial
    # Check if trial has to be stopped or no
    if (stop_trial) {
      new_dose <- NA
      message("The trial stopped based on the stopping rule.\n \n")
    } else { # if no stopping
      new_dose <- which(proba_targeted_tox[candidate_doses] == max(proba_targeted_tox[candidate_doses]))
      if (is.null(candidate_doses)){
        new_dose <- 1
      }
    }
    
    dose_names <- as.vector(sapply(1:length(doses), function (i) paste("Dose", i)))
    ### Compute probability of toxicity 95% credible intervals for each of the dose based on the 2.5- and 97.5- percentiles for the parameters
    credible_interval_per_dose <- matrix(sapply(1:length(doses), function(i) quantile(p_estim_sum[,i], probs = c(0.025, 0.975))), 
                                         nrow = 2,
                                         ncol = length(doses), 
                                         dimnames = list(c("2.5%", "97.5%"), dose_names),
                                         byrow = FALSE)
    
    # Display for rstan estimated parameters and pop-PK estimated parameters using SAEMIX
    params <- c(beta2, log_beta3)
    names(params) <- c("beta2", "log_beta3")
    estimated_PK_parameters <- c(updated_ka, updated_Cl_pop, updated_V_pop, updated_omega_Cl, updated_omega_V)
    
  }
  return(list(new_dose = new_dose, auc_s = auc_s, p_estim_dose_mean = p_estim_dose_mean, params = params, credible_interval_per_dose = credible_interval_per_dose, ref_log_AUC = simu_ref_log_AUC, estimated_PK_parameters = estimated_PK_parameters))
}