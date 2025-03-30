##############################################################################
### Execution of the simulation scenario and display of simulation results ###
##############################################################################

# Core function
run_sim_scenario <- function(N, doses, ref_dose, n_trials, cohort_size, targeted_tox, PK_parameters, omega_IIV, omega_alpha, 
                             n_sampling_timepoints, real_sampling, cv, compartmental_model, elimination_kinetics, exposure_metric, 
                             AUC_method, exposure_tox_threshold, options, proba_safety, seed, scenario_id) {

  ########################################################################################################
  ### Computation of patients' drug concentration using simulated PK data based on scenario parameters ###
  ########################################################################################################
  # Simulated PK data for each considered trial to be used inside the toxicity estimation model
  simulated_PK_data <<- simu_data(
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    doses = doses,
    exposure_tox_threshold = exposure_tox_threshold,
    n_sampling_timepoints = n_sampling_timepoints,
    time_sampling = time_sampling,
    N = N,
    n_trials = n_trials,
    exposure_metric = exposure_metric,
    AUC_method = AUC_method,
    compartmental_model = compartmental_model,
    elimination_kinetics = elimination_kinetics
  )
  
  ##########################################################################
  #### Bayesian Logistic Regression Method (BLRM)
  ##########################################################################
  res_blrm <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "blrm"    
  )
  
  ##########################################################################
  #### Ursino (2017) PKLOGIT
  ##########################################################################
  res_pklogit <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "pklogit"
  )
  
  ##########################################################################
  #### Micallef (2018) ED-EWOC
  ##########################################################################
  res_ed_ewoc <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "ed_ewoc"
  )
  
  ##############################################################################
  #### ED (without Micallef (2018) EWOC decision rule)
  ##############################################################################
  res_ed <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "ed"
  )
  
  ##########################################################################
  #### Takeda (2018) PK-CRM
  ##########################################################################
  res_pk_crm <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "pk_crm"
  )
  
  ##########################################################################
  #### Günhan (2021) NAIVE TITE-PK
  ##########################################################################
  res_naive_tite_pk <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "naive_tite_pk"
  )
  
  ##########################################################################
  #### Günhan (2021) INFORMED TITE-PK
  ##########################################################################
  res_informed_tite_pk <<- run_trial_dosefinding(
    N = N,
    doses = doses,
    ref_dose = ref_dose,
    cohort_size = cohort_size,
    real_sampling = real_sampling,
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    simulated_PK_data = simulated_PK_data,
    exposure_tox_threshold = exposure_tox_threshold,
    n_trials = n_trials,
    targeted_tox = targeted_tox,
    proba_safety = proba_safety,
    options = options,
    method = "informed_tite_pk"
  )
  
  
  ##########################################################################
  #### Display the appropriate summary for the results of the simulations
  ##########################################################################
  if (n_trials > 1) {
    scenario <- scenario_id
    models_name <- c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "PK-CRM", "Naive TITE-PK", "Informed TITE-PK")
    dose_selection <- matrix(c(round(res_blrm@new_dose, digits = 3), round(res_pklogit@new_dose, digits = 3), 
                               round(res_ed_ewoc@new_dose, digits = 3), round(res_ed@new_dose, digits = 3),
                               round(res_pk_crm@new_dose, digits = 3), round(res_naive_tite_pk@new_dose, digits = 3),
                               round(res_informed_tite_pk@new_dose, digits = 3)), 
                             nrow = length(models_name), ncol = length(doses)+1, byrow = TRUE)
    allocation_blrm <- rep(NA, length(doses)+1)
    allocation_pklogit <- rep(NA, length(doses)+1)
    allocation_ed_ewoc <- rep(NA, length(doses)+1)
    allocation_ed <- rep(NA, length(doses)+1)
    allocation_pk_crm <- rep(NA, length(doses)+1)
    allocation_naive_tite_pk <- rep(NA, length(doses)+1)
    allocation_informed_tite_pk <- rep(NA, length(doses)+1)
    for(i in 1:length(doses)) {
      n_levels_blrm <- length(which(res_blrm@dose_levels == i))
      n_levels_pklogit <- length(which(res_pklogit@dose_levels == i))
      n_levels_ed_ewoc <- length(which(res_ed_ewoc@dose_levels == i))
      n_levels_ed <- length(which(res_ed@dose_levels == i))
      n_levels_pk_crm <- length(which(res_pk_crm@dose_levels == i))
      n_levels_naive_tite_pk <- length(which(res_naive_tite_pk@dose_levels == i))
      n_levels_informed_tite_pk <- length(which(res_informed_tite_pk@dose_levels == i))
      allocation_blrm[i+1] <- round(n_levels_blrm/length(which(res_blrm@dose_levels != 0)), digits=3)
      allocation_pklogit[i+1] <- round(n_levels_pklogit/length(which(res_pklogit@dose_levels != 0)), digits=3)
      allocation_ed_ewoc[i+1] <- round(n_levels_ed_ewoc/length(which(res_ed_ewoc@dose_levels != 0)), digits=3)
      allocation_ed[i+1] <- round(n_levels_ed/length(which(res_ed@dose_levels != 0)), digits=3)
      allocation_pk_crm[i+1] <- round(n_levels_pk_crm/length(which(res_pk_crm@dose_levels != 0)), digits=3)
      allocation_naive_tite_pk[i+1] <- round(n_levels_naive_tite_pk/length(which(res_naive_tite_pk@dose_levels != 0)), digits=3)
      allocation_informed_tite_pk[i+1] <- round(n_levels_informed_tite_pk/length(which(res_informed_tite_pk@dose_levels != 0)), digits=3)
    }
    allocation_blrm[1] <- "NA"
    allocation_pklogit[1] <- "NA"
    allocation_ed_ewoc[1]  <- "NA"
    allocation_ed[1]  <- "NA"
    allocation_pk_crm[1]  <- "NA"
    allocation_naive_tite_pk[1]  <- "NA"
    allocation_informed_tite_pk[1]  <- "NA"
    dose_allocation <- matrix(c(allocation_blrm, allocation_pklogit, allocation_ed_ewoc, allocation_ed, allocation_pk_crm, 
                                allocation_naive_tite_pk, allocation_informed_tite_pk), nrow = length(models_name), 
                              ncol = length(doses)+1, byrow = TRUE)
    sum_trial_blrm <- rep(NA, n_trials)
    sum_trial_pklogit <- rep(NA, n_trials)
    sum_trial_ed_ewoc <- rep(NA, n_trials)
    sum_trial_ed <- rep(NA, n_trials)
    sum_trial_pk_crm <- rep(NA, n_trials)
    sum_trial_naive_tite_pk <- rep(NA, n_trials)
    sum_trial_informed_tite_pk <- rep(NA, n_trials)
    for (tr in 1:n_trials) {
      sum_trial_blrm[tr] <- sum(res_blrm@toxicity[tr,], na.rm = TRUE)
      sum_trial_pklogit[tr] <- sum(res_pklogit@toxicity[tr,], na.rm = TRUE)
      sum_trial_ed_ewoc[tr] <- sum(res_ed_ewoc@toxicity[tr,], na.rm = TRUE)
      sum_trial_ed[tr] <- sum(res_ed@toxicity[tr,], na.rm = TRUE)
      sum_trial_pk_crm[tr] <- sum(res_pk_crm@toxicity[tr,], na.rm = TRUE)
      sum_trial_naive_tite_pk[tr] <- sum(res_naive_tite_pk@toxicity[tr,], na.rm = TRUE)
      sum_trial_informed_tite_pk[tr] <- sum(res_informed_tite_pk@toxicity[tr,], na.rm = TRUE)
    }
    median_DLTs <- c(round(median(sum_trial_blrm), 0), round(median(sum_trial_pklogit), 0), round(median(sum_trial_ed_ewoc), 0), 
                     round(median(sum_trial_ed), 0), round(median(sum_trial_pk_crm), 0), 
                     round(median(sum_trial_naive_tite_pk), 0), round(median(sum_trial_informed_tite_pk), 0))
    min_DLTs <- c(min(sum_trial_blrm), min(sum_trial_pklogit), min(sum_trial_ed_ewoc), min(sum_trial_ed), 
                  min(sum_trial_pk_crm), min(sum_trial_naive_tite_pk), min(sum_trial_informed_tite_pk))
    max_DLTs <- c(max(sum_trial_blrm), max(sum_trial_pklogit), max(sum_trial_ed_ewoc), max(sum_trial_ed), 
                  max(sum_trial_pk_crm), max(sum_trial_naive_tite_pk), max(sum_trial_informed_tite_pk))
    simulation_results <<- new("simulation_res", scenario = scenario, PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv,
                               exposure_metric = exposure_metric, prior_tox_real = prior_tox_real, N = N, cohort_size = cohort_size, n_trials = n_trials, 
                               n_sampling_timepoints = n_sampling_timepoints, real_time_sampling = time_sampling[real_sampling], doses = doses, 
                               exposure_tox_threshold = exposure_tox_threshold, n_chains = options$n_chains, n_iter = options$n_iter, n_adapt = options$n_adapt,
                               n_cores = options$n_cores, dose_allocation = dose_allocation, dose_selection = dose_selection,
                               targeted_tox = targeted_tox, median_DLTs = median_DLTs, min_DLTs = min_DLTs, max_DLTs = max_DLTs, 
                               models = models_name, seed = res_pklogit@seed)
    
    
    ##########################################################################
    #### Display the appropriate plots for the results of the simulations
    ##########################################################################
    
    # Confirmation plot regarding estimated probabilities of toxicity for each trial and each dose-finding design
    models_name <- c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "Naive TITE-PK", "Informed TITE-PK")
    #models_name <- c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "PK-CRM", "Naive TITE-PK", "Informed TITE-PK")
    MTD_in_scenario <- ifelse(prior_tox_real[1] < 0.40, which(abs(prior_tox_real - targeted_tox) == min(abs(prior_tox_real - targeted_tox))), 0)
    if(1 <= scenario_id & scenario_id <= 5) {
      scenario_name <- paste0("A", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (6 <= scenario_id & scenario_id <= 10) {
      scenario_name <- paste0("B", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (11 <= scenario_id & scenario_id <= 15) {
      scenario_name <- paste0("C", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (16 <= scenario_id & scenario_id <= 20) {
      scenario_name <- paste0("D", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (21 <= scenario_id & scenario_id <= 25) {
      scenario_name <- paste0("E", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (scenario_id == 26) {
      scenario_name <- paste0("F", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (27 <= scenario_id & scenario_id <= 31) {
      scenario_name <- paste0("G", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (31 <= scenario_id & scenario_id <= 36) {
      scenario_name <- paste0("H", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else if (37 <= scenario_id & scenario_id <= 42) {
      scenario_name <- paste0("I", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    } else {
      scenario_name <- paste0("X", as.character(ifelse(MTD_in_scenario == 0, length(doses) + 1, MTD_in_scenario)))
    }
    dose_vector <- sapply(1:length(doses), function(i) paste("Dose", i))
    names(dose_vector) <- paste0("dose", 1:length(doses))
    
    post_proba_tox_blrm <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("BLRM", length(doses)*n_trials), 
      dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
      as.data.frame(matrix(unlist(res_blrm@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    )
    post_proba_tox_pklogit <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("PKLOGIT", length(doses)*n_trials), 
      dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
      as.data.frame(matrix(unlist(res_pklogit@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    )
    post_proba_tox_ed_ewoc <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("ED-EWOC", length(doses)*n_trials),
      dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
      as.data.frame(matrix(unlist(res_ed_ewoc@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    )
    post_proba_tox_ed <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("ED", length(doses)*n_trials),
      dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
      as.data.frame(matrix(unlist(res_ed@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    )
    #post_proba_tox_pk_crm <- cbind(
    #trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
    #method_name = rep("PK-CRM", length(doses)*n_trials),
    #dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
    #as.data.frame(matrix(unlist(res_pk_crm@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    #)
    post_proba_tox_naive_tite_pk <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("Naive TITE-PK", length(doses)*n_trials),
      dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
      as.data.frame(matrix(unlist(res_naive_tite_pk@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    )
    post_proba_tox_informed_tite_pk <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("Informed TITE-PK", length(doses)*n_trials),
      dose_num = rep(sapply(1:length(doses), function(i) paste0("dose", i)), n_trials),
      as.data.frame(matrix(unlist(res_informed_tite_pk@p_estim), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox")))
    )
    aggregate_all_df <- list(post_proba_tox_blrm, post_proba_tox_pklogit, post_proba_tox_ed_ewoc, post_proba_tox_ed, post_proba_tox_naive_tite_pk, post_proba_tox_informed_tite_pk)
    #aggregate_all_df <- list(post_proba_tox_blrm, post_proba_tox_pklogit, post_proba_tox_ed_ewoc, post_proba_tox_ed, post_proba_tox_pk_crm, post_proba_tox_naive_tite_pk, post_proba_tox_informed_tite_pk)
    post_proba_tox_all <- Reduce(function(x, y) merge(x, y, all=TRUE), aggregate_all_df)  
    post_proba_tox_all$method_name <- factor(post_proba_tox_all$method_name, levels = models_name)
    post_proba_tox_all$dose_num <- factor(post_proba_tox_all$dose_num)
    post_mean_tox_scen <- data.frame(dose_num = sapply(1:length(doses), function(i) paste0("dose", i)), prior_proba_tox = res_blrm@prior_tox_real)
    post_proba_tox_plot <<- ggplot() +
      geom_boxplot(data = post_proba_tox_all, aes(x = dose_num, y = proba_tox, fill = method_name)) +
      geom_point(data = post_mean_tox_scen, aes(y = prior_proba_tox, x = dose_num, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Targeted toxicity"), color = "red", linewidth = 1) +
      scale_x_discrete(labels = dose_vector) +
      labs(title = paste("Scenario", scenario_name, "- Estimated probabilities of toxicity from simulated trials for each dose-finding design"), y = "Probability of toxicity", x = "Panel of doses", fill = "Dose-finding designs") +
      scale_fill_brewer(palette="Dark2") +
      scale_linetype_manual(name = "", values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(name = "", values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ggsave(paste0("plots/scenario_", scenario_name, "_PPT_plot.png"), plot = post_proba_tox_plot, dpi = "retina", units = "in", width = 14, height = 8)
    
    # Confirmation plot using 95% credibility intervals of estimated probabilities of toxicity for each trial and each dose-finding design
    BLRM_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_blrm@credible_intervals[[i]], mean = as.vector(res_blrm@p_estim[[i]])))
    post_proba_tox_CI_blrm <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("BLRM", length(doses)*n_trials), 
      dose_num = factor(rep(dose_vector, n_trials)),
      Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
      as.data.frame(matrix(unlist(BLRM_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(BLRM_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(BLRM_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    )
    PKLOGIT_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_pklogit@credible_intervals[[i]], mean = as.vector(res_pklogit@p_estim[[i]])))
    post_proba_tox_CI_pklogit <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("PKLOGIT", length(doses)*n_trials), 
      dose_num = factor(rep(dose_vector, n_trials)),
      Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
      as.data.frame(matrix(unlist(PKLOGIT_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(PKLOGIT_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(PKLOGIT_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    )
    ED_EWOC_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_ed_ewoc@credible_intervals[[i]], mean = as.vector(res_ed_ewoc@p_estim[[i]])))
    post_proba_tox_CI_ed_ewoc <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("ED-EWOC", length(doses)*n_trials),
      dose_num = factor(rep(dose_vector, n_trials)),
      Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
      as.data.frame(matrix(unlist(ED_EWOC_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(ED_EWOC_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(ED_EWOC_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    )
    ED_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_ed@credible_intervals[[i]], mean = as.vector(res_ed@p_estim[[i]])))
    post_proba_tox_CI_ed <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("ED", length(doses)*n_trials),
      dose_num = factor(rep(dose_vector, n_trials)),
      Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
      as.data.frame(matrix(unlist(ED_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(ED_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(ED_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    )
    #PK_CRM_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_pk_crm@credible_intervals[[i]], mean = as.vector(res_pk_crm@p_estim[[i]])))
    #post_proba_tox_CI_pk_crm <- cbind(
    #trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
    #method_name = rep("PK-CRM", length(doses)*n_trials), 
    #dose_num = factor(rep(dose_vector, n_trials)),
    #Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
    #as.data.frame(matrix(unlist(PK_CRM_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(PK_CRM_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(PK_CRM_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    #)
    Naive_TITE_PK_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_naive_tite_pk@credible_intervals[[i]], mean = as.vector(res_naive_tite_pk@p_estim[[i]])))
    post_proba_tox_CI_naive_tite_pk <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("Naive TITE-PK", length(doses)*n_trials),
      dose_num = factor(rep(dose_vector, n_trials)),
      Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
      as.data.frame(matrix(unlist(Naive_TITE_PK_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(Naive_TITE_PK_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(Naive_TITE_PK_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    )
    Informed_TITE_PK_proba_tox_CI_and_mean <- lapply(1:n_trials, function(i) rbind(res_informed_tite_pk@credible_intervals[[i]], mean = as.vector(res_informed_tite_pk@p_estim[[i]])))
    post_proba_tox_CI_informed_tite_pk <- cbind(
      trial_num = as.vector(sapply(1:n_trials, function(i) rep(i, length(doses)))), 
      method_name = rep("Informed TITE-PK", length(doses)*n_trials),
      dose_num = factor(rep(dose_vector, n_trials)),
      Xth_pct = factor(rep(c("2.5th percentile", "97.5th percentile", "mean"), n_trials*length(doses)), levels = c("2.5th percentile", "mean", "97.5th percentile")),
      as.data.frame(matrix(unlist(Informed_TITE_PK_proba_tox_CI_and_mean), nrow = n_trials*length(doses)*nrow(Informed_TITE_PK_proba_tox_CI_and_mean[[1]]), ncol = 1, byrow = TRUE, dimnames = list(c(1:(nrow(Informed_TITE_PK_proba_tox_CI_and_mean[[1]])*length(doses)*n_trials)), "CI_proba_tox")))
    )
    
    BLRM_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_blrm, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      geom_density_ridges() +
      scale_y_discrete(labels = dose_vector) +
      labs(title = "BLRM", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      scale_fill_brewer(palette="Set2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    PKLOGIT_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_pklogit, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      geom_density_ridges() +
      scale_y_discrete(labels = dose_vector) +
      labs(title = "PKLOGIT", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      scale_fill_brewer(palette="Set2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_EWOC_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_ed_ewoc, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      geom_density_ridges() +
      scale_y_discrete(labels = dose_vector) +
      labs(title = "ED-EWOC", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      scale_fill_brewer(palette="Set2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_ed, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      geom_density_ridges() +
      scale_y_discrete(labels = dose_vector) +
      labs(title = "ED", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      scale_fill_brewer(palette="Set2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    #PK_CRM_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_pk_crm, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      #geom_density_ridges() +
      #scale_y_discrete(labels = dose_vector) +
      #labs(title = "PK-CRM", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      #scale_fill_brewer(palette="Set2") +
      #theme_bw() +
      #theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    Naive_TITE_PK_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_naive_tite_pk, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      geom_density_ridges() +
      scale_y_discrete(labels = dose_vector) +
      labs(title = "Naive TITE-PK", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      scale_fill_brewer(palette="Set2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    Informed_TITE_PK_post_proba_tox_CI_plot <- ggplot(data = post_proba_tox_CI_informed_tite_pk, mapping = aes(y = dose_num, x = CI_proba_tox, fill = Xth_pct)) +
      geom_density_ridges() +
      scale_y_discrete(labels = dose_vector) +
      labs(title = "Informed TITE-PK", x = "Probability of toxicity", y = "Panel of doses", fill = "Distributions for estimated probabilities of toxicity") +
      scale_fill_brewer(palette="Set2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    post_proba_tox_CI_plot <<- ggarrange(BLRM_post_proba_tox_CI_plot + rremove("xlab"), PKLOGIT_post_proba_tox_CI_plot + rremove("xlab") + rremove("ylab"), ED_EWOC_post_proba_tox_CI_plot + rremove("xlab"), ED_post_proba_tox_CI_plot + rremove("xlab") + rremove("ylab"), Naive_TITE_PK_post_proba_tox_CI_plot, Informed_TITE_PK_post_proba_tox_CI_plot + rremove("ylab"), 
                                         common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
    post_proba_tox_CI_plot <<- annotate_figure(post_proba_tox_CI_plot, top = text_grob(paste("Scenario", scenario_name, "- Density plots for estimated toxicity probabilities and for lower and upper bounds of 95% credible interval"), color = "black", face = "bold", size = 14))
    ggsave(paste0("plots/scenario_", scenario_name, "_PPT_CI_plot.png"), plot = wrap_elements(post_proba_tox_CI_plot), dpi = "retina", units = "in", width = 14, height = 8)
  
    # Visualization bar plot of the percentage of MTD selection for each possible dose and each dose-finding design
    color_palette <- c("#A3A500", "#00BF7D", "#00B0F6", "#E76BF3", "#D89000", "#F8766D")
    dose_selection_scenario <- data.frame(design = c(rep("BLRM", length(doses)+1),
                                                           rep("PKLOGIT", length(doses)+1), 
                                                           rep("ED-EWOC", length(doses)+1),
                                                           rep("ED", length(doses)+1),
                                                           rep("Naive TITE-PK", length(doses)+1),
                                                           rep("Informed TITE-PK", length(doses)+1)
    ),
    dose = rep(c("STOP", dose_vector), length(models_name)),
    selection = c(t(simulation_results@dose_selection[-5,]))
    )
    #dose_selection_scenario <- data.frame(design = c(rep("BLRM", length(doses)+1),
                                                     #rep("PKLOGIT", length(doses)+1), 
                                                     #rep("ED-EWOC", length(doses)+1),
                                                     #rep("ED", length(doses)+1),
                                                     #rep("PK-CRM", length(doses)+1),
                                                     #rep("Naive TITE-PK", length(doses)+1),
                                                     #rep("Informed TITE-PK", length(doses)+1)
    #),
    #dose = rep(c("STOP", dose_vector), length(models_name)),
    #selection = c(t(simulation_results@dose_selection))
    #)
    
    dose_selection_scenario$design <- factor(dose_selection_scenario$design, levels = c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "Naive TITE-PK", "Informed TITE-PK"))
    #dose_selection_scenario$design <- factor(dose_selection_scenario$design, levels = c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "PK-CRM", "Naive TITE-PK", "Informed TITE-PK"))
    dose_selection_scenario$dose <- factor(dose_selection_scenario$dose, levels = c("STOP", dose_vector))
    if (1 <= MTD_in_scenario & length(doses) >= MTD_in_scenario) {
      real_tox_prob_values <- sapply(1:(length(doses)-1), function(i) ifelse(i == MTD_in_scenario, paste0("\\textbf{", round(prior_tox_real[i], 3), "}", " |"), paste(round(prior_tox_real[i], 3), " |")))
      real_tox_prob_values <- c(real_tox_prob_values, ifelse(length(doses) == MTD_in_scenario, paste0("\\textbf{", round(prior_tox_real[length(doses)], 3), "}"), paste(round(prior_tox_real[i], 3))))
      real_tox_prob_string <- paste(real_tox_prob_values, collapse = " ")
      title_selection_plot <- latex2exp::TeX(paste("\\textbf{Scenario", scenario_name, "- Real probabilities of toxicity (MTD = Dose", MTD_in_scenario, "):}", real_tox_prob_string))
    } else {
      real_tox_prob_values <- sapply(1:(length(doses)-1), function(i) paste(round(prior_tox_real[i], 3), " |"))
      real_tox_prob_values <- c(real_tox_prob_values, paste(round(prior_tox_real[i], 3)))
      real_tox_prob_string <- paste(real_tox_prob_values, collapse = " ")
      title_selection_plot <- latex2exp::TeX(paste("\\textbf{Scenario", scenario_name, "- Real probabilities of toxicity (MTD = None):}", real_tox_prob_string))
    }
    dose_selection_scenario_barplot <<- ggplot(dose_selection_scenario, aes(fill = dose, y = selection*100, x = design)) +
      geom_bar(position="dodge", stat="identity") +
      scale_fill_manual(values = c(tail(color_palette, n = 1), color_palette[1:length(doses)])) +
      labs(x = "Dose-finding design", y = "Pourcentage of dose selection (%)", fill = "Selected MTD", title = title_selection_plot) +
      theme_bw() +
      theme(legend.title = element_text(face = "bold"))
    ggsave(paste0("plots/scenario_", scenario_name, "_dose_selection_plot.png"), plot = dose_selection_scenario_barplot, dpi = "retina", units = "in", width = 14, height = 8)
    
    # Visualization bar plot of the percentage of dose allocation for each dose-finding design
    dose_allocation_scenario <- data.frame(design = c(rep("BLRM", length(doses)),
                                                     rep("PKLOGIT", length(doses)), 
                                                     rep("ED-EWOC", length(doses)),
                                                     rep("ED", length(doses)),
                                                     rep("Naive TITE-PK", length(doses)),
                                                     rep("Informed TITE-PK", length(doses))
    ),
    dose = rep(dose_vector, length(models_name)),
    allocation = as.numeric(c(t(simulation_results@dose_allocation[-5,-1])))
    )
    #dose_allocation_scenario <- data.frame(design = c(rep("BLRM", length(doses)),
                                                      #rep("PKLOGIT", length(doses)), 
                                                      #rep("ED-EWOC", length(doses)),
                                                      #rep("ED", length(doses)),
                                                      #rep("PK-CRM", length(doses)),
                                                      #rep("Naive TITE-PK", length(doses)),
                                                      #rep("Informed TITE-PK", length(doses))
    #),
    #dose = rep(dose_vector, length(models_name)),
    #allocation = as.numeric(c(t(simulation_results@dose_allocation[-5,-1])))
    #)
    
    dose_allocation_scenario$design <- factor(dose_allocation_scenario$design, levels = c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "Naive TITE-PK", "Informed TITE-PK"))
    #dose_allocation_scenario$design <- factor(dose_allocation_scenario$design, levels = c("BLRM", "PKLOGIT", "ED-EWOC", "ED", "PK-CRM", "Naive TITE-PK", "Informed TITE-PK"))
    dose_allocation_scenario$dose <- factor(dose_allocation_scenario$dose, levels = dose_vector)
    dose_allocation_scenario_barplot <<- ggplot(dose_allocation_scenario, aes(fill = dose, y = allocation*100, x = design)) +
      geom_bar(position="dodge", stat="identity") +
      scale_fill_manual(values = color_palette[1:length(doses)]) +
      labs(x = "Dose-finding design", y = "Pourcentage of dose allocation (%)", fill = "Allocated dose", title = paste("Scenario", scenario_name)) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ggsave(paste0("plots/scenario_", scenario_name, "_dose_allocation_plot.png"), plot = dose_allocation_scenario_barplot, dpi = "retina", units = "in", width = 14, height = 8)
    
    # Confirmation plot for estimated PK population parameters and random effects variances by ED-EWOC pop-PK model
    ED_EWOC_popPK_estimated_parameters_names <- c("ka", "Cl_pop", "V_pop", "omega_Cl", "omega_V")
    #ED_EWOC_popPK_estimated_parameters_names <- c("Cl_pop", "V_pop", "omega_Cl", "omega_V")
    ED_EWOC_popPK_estimated_parameters_aggregated <- res_ed_ewoc@popPK_estimated_parameters
    ED_EWOC_popPK_estimated_parameters_aggregated <- lapply(1:length(ED_EWOC_popPK_estimated_parameters_aggregated), function (i) ED_EWOC_popPK_estimated_parameters_aggregated[[i]][complete.cases(ED_EWOC_popPK_estimated_parameters_aggregated[[i]]),])
    ED_EWOC_num_notdf <- which(sapply(1:length(ED_EWOC_popPK_estimated_parameters_aggregated), function(i) !is.matrix(ED_EWOC_popPK_estimated_parameters_aggregated[[i]])))
    for(i in ED_EWOC_num_notdf) {
      ED_EWOC_popPK_estimated_parameters_aggregated[[i]] <- matrix(as.vector(ED_EWOC_popPK_estimated_parameters_aggregated[[i]]), nrow = 1, ncol = length(ED_EWOC_popPK_estimated_parameters_aggregated[[i]]), byrow = TRUE, dimnames = list(1, ED_EWOC_popPK_estimated_parameters_names))
    }
    ED_EWOC_num_trials_notempty <- which(lapply(ED_EWOC_popPK_estimated_parameters_aggregated, nrow) != 0)
    ED_EWOC_popPK_estimated_parameters_aggregated_clear <- ED_EWOC_popPK_estimated_parameters_aggregated[ED_EWOC_num_trials_notempty]
    ED_EWOC_popPK_estimated_parameters_all_raw <- Reduce(function(x, y) merge(x, y, all=TRUE), ED_EWOC_popPK_estimated_parameters_aggregated_clear)
    ED_EWOC_popPK_estimated_parameters_all <- data.frame(
      trial_num = rep(unlist(sapply(ED_EWOC_num_trials_notempty, function(i) rep(i, nrow(ED_EWOC_popPK_estimated_parameters_aggregated[[i]])))), length(ED_EWOC_popPK_estimated_parameters_names)),
      PK_parameters = as.vector(sapply(ED_EWOC_popPK_estimated_parameters_names, function(i) {
        total_length <- sum(sapply(1:length(ED_EWOC_popPK_estimated_parameters_aggregated_clear), function(j) nrow(ED_EWOC_popPK_estimated_parameters_aggregated_clear[[j]])))
        return(rep(i, total_length))
      })),
      estimated_value = as.vector(sapply(ED_EWOC_popPK_estimated_parameters_names, function(i) ED_EWOC_popPK_estimated_parameters_all_raw[, i]))
    )
    ED_EWOC_popPK_estimated_parameters_all$PK_parameters <- factor(ED_EWOC_popPK_estimated_parameters_all$PK_parameters, levels = ED_EWOC_popPK_estimated_parameters_names)
    
    ED_EWOC_popPK_ka_plot <- ggplot() +
      geom_violin(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "ka", y = ka), fill = "deepskyblue") +
      geom_boxplot(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "ka", y = ka), fill = "deepskyblue", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = PK_parameters[1], lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated (population) $k_{a}$}"), y = "Estimated value", x = "", lty = NULL) +
      scale_x_discrete(labels = unname(TeX("$k_{a}$"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_EWOC_popPK_Cl_pop_plot <- ggplot() +
      geom_violin(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "Cl_pop", y = Cl_pop), fill = "brown1") +
      geom_boxplot(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "Cl_pop", y = Cl_pop), fill = "brown1", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = PK_parameters[2], lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated $Cl_{pop}$}"), y = "Estimated value", x = "", lty = NULL) +
      scale_x_discrete(labels = unname(TeX("$Cl_{pop}$"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_EWOC_popPK_V_pop_plot <- ggplot() +
      geom_violin(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "V_pop", y = V_pop), fill = "darkseagreen1") +
      geom_boxplot(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "V_pop", y = V_pop), fill = "darkseagreen1", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = PK_parameters[3], lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated $V_{pop}$}"), y = "Estimated value", x = "",  lty = NULL) +
      scale_x_discrete(labels = unname(TeX("$V_{pop}$"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_EWOC_popPK_omega_pop_plot <- ggplot() +
      geom_violin(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_Cl", y = omega_Cl), fill = "brown1") +
      geom_boxplot(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_Cl", y = omega_Cl), fill = "brown1", width = 0.1, col = "lightgrey") +
      geom_violin(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_V", y = omega_V), fill = "darkseagreen1") +
      geom_boxplot(data = ED_EWOC_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_V", y = omega_V), fill = "darkseagreen1", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = omega_IIV, lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      scale_x_discrete(labels = unname(TeX(c("$omega_{Cl}$", "$omega_{V}$")))) +
      labs(title = TeX("\\textbf{Pop-PK estimated $omega_{Cl}$ and $omega_{V}$}"), y = "Estimated value", x = "", lty = NULL) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_EWOC_popPK_plot <<- ggarrange(ED_EWOC_popPK_ka_plot, ED_EWOC_popPK_Cl_pop_plot + rremove("ylab"), ED_EWOC_popPK_V_pop_plot, ED_EWOC_popPK_omega_pop_plot + rremove("ylab"), common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 2)
    ED_EWOC_popPK_plot <<- annotate_figure(ED_EWOC_popPK_plot, top = text_grob(paste("Scenario", scenario_name, "- Estimated parameters of ED-EWOC population PK model for all trials"), color = "black", face = "bold", size = 14))
    ggsave(paste0("plots/scenario_", scenario_name, "_ED_EWOC_popPK_plot.png"), plot = wrap_elements(ED_EWOC_popPK_plot), dpi = "retina", units = "in", width = 14, height = 8, limitsize = FALSE)
    
    # Confirmation plot for estimated PK population parameters and random effects variances by ED pop-PK model
    ED_popPK_estimated_parameters_names <- c("ka", "Cl_pop", "V_pop", "omega_Cl", "omega_V")
    #ED_popPK_estimated_parameters_names <- c("Cl_pop", "V_pop", "omega_Cl", "omega_V")
    ED_popPK_estimated_parameters_aggregated <- res_ed@popPK_estimated_parameters
    ED_popPK_estimated_parameters_aggregated <- lapply(1:length(ED_popPK_estimated_parameters_aggregated), function (i) ED_popPK_estimated_parameters_aggregated[[i]][complete.cases(ED_popPK_estimated_parameters_aggregated[[i]]),])
    ED_num_notdf <- which(sapply(1:length(ED_popPK_estimated_parameters_aggregated), function(i) !is.matrix(ED_popPK_estimated_parameters_aggregated[[i]])))
    for(i in ED_num_notdf) {
      ED_popPK_estimated_parameters_aggregated[[i]] <- matrix(as.vector(ED_popPK_estimated_parameters_aggregated[[i]]), nrow = 1, ncol = length(ED_popPK_estimated_parameters_aggregated[[i]]), byrow = TRUE, dimnames = list(1, ED_popPK_estimated_parameters_names))
    }
    ED_num_trials_notempty <- which(lapply(ED_popPK_estimated_parameters_aggregated, nrow) != 0)
    ED_popPK_estimated_parameters_aggregated_clear <- ED_popPK_estimated_parameters_aggregated[ED_num_trials_notempty]
    ED_popPK_estimated_parameters_all_raw <- Reduce(function(x, y) merge(x, y, all=TRUE), ED_popPK_estimated_parameters_aggregated_clear)
    ED_popPK_estimated_parameters_all <- data.frame(
      trial_num = rep(unlist(sapply(ED_num_trials_notempty, function(i) rep(i, nrow(ED_popPK_estimated_parameters_aggregated[[i]])))), length(ED_popPK_estimated_parameters_names)),
      PK_parameters = as.vector(sapply(ED_popPK_estimated_parameters_names, function(i) {
        total_length <- sum(sapply(1:length(ED_popPK_estimated_parameters_aggregated_clear), function(j) nrow(ED_popPK_estimated_parameters_aggregated_clear[[j]])))
        return(rep(i, total_length))
      })),
      estimated_value = as.vector(sapply(ED_popPK_estimated_parameters_names, function(i) ED_popPK_estimated_parameters_all_raw[, i]))
    )
    ED_popPK_estimated_parameters_all$PK_parameters <- factor(ED_popPK_estimated_parameters_all$PK_parameters, levels = ED_popPK_estimated_parameters_names)
    
    ED_popPK_ka_plot <- ggplot() +
      geom_violin(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "ka", y = ka), fill = "deepskyblue") +
      geom_boxplot(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "ka", y = ka), fill = "deepskyblue", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = PK_parameters[1], lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated (population) $k_{a}$}"), y = "Estimated value", x = "", lty = NULL) +
      scale_x_discrete(labels = unname(TeX("$k_{a}$"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_popPK_Cl_pop_plot <- ggplot() +
      geom_violin(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "Cl_pop", y = Cl_pop), fill = "brown1") +
      geom_boxplot(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "Cl_pop", y = Cl_pop), fill = "brown1", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = PK_parameters[2], lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated $Cl_{pop}$}"), y = "Estimated value", x = "", lty = NULL) +
      scale_x_discrete(labels = unname(TeX("$Cl_{pop}$"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_popPK_V_pop_plot <- ggplot() +
      geom_violin(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "V_pop", y = V_pop), fill = "darkseagreen1") +
      geom_boxplot(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "V_pop", y = V_pop), fill = "darkseagreen1", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = PK_parameters[3], lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated $V_{pop}$}"), y = "Estimated value", x = "",  lty = NULL) +
      scale_x_discrete(labels = unname(TeX("$V_{pop}$"))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_popPK_omega_pop_plot <- ggplot() +
      geom_violin(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_Cl", y = omega_Cl), fill = "brown1") +
      geom_boxplot(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_Cl", y = omega_Cl), fill = "brown1", width = 0.1, col = "lightgrey") +
      geom_violin(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_V", y = omega_V), fill = "darkseagreen1") +
      geom_boxplot(data = ED_popPK_estimated_parameters_all_raw, mapping = aes(x = "omega_V", y = omega_V), fill = "darkseagreen1", width = 0.1, col = "lightgrey") +
      geom_hline(aes(yintercept = omega_IIV, lty = "True parameter value"), col = "red") +
      scale_linetype_manual(values = 1) +
      labs(title = TeX("\\textbf{Pop-PK estimated $omega_{Cl}$ and $omega_{V}$}"), y = "Estimated value", x = "", lty = NULL) +
      scale_x_discrete(labels = unname(TeX(c("$omega_{Cl}$", "$omega_{V}$")))) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
    ED_popPK_plot <<- ggarrange(ED_popPK_ka_plot, ED_popPK_Cl_pop_plot + rremove("ylab"), ED_popPK_V_pop_plot, ED_popPK_omega_pop_plot + rremove("ylab"), common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 2)
    ED_popPK_plot <<- annotate_figure(ED_popPK_plot, top = text_grob(paste("Scenario", scenario_name, "- Estimated parameters of ED population PK model for all trials"), color = "black", face = "bold", size = 14))
    ggsave(paste0("plots/scenario_", scenario_name, "_ED_popPK_plot.png"), plot = wrap_elements(ED_popPK_plot), dpi = "retina", units = "in", width = 14, height = 8, limitsize = FALSE)
    
    # Return formatted display after running simulation scenario for multiple trials to highlight important results
    return(simulation_results)
  } else {
    # Return formatted display after running simulation scenario for one clinical trial to highlight important results
    return(c(res_blrm, res_pklogit, res_ed_ewoc, res_ed, res_pk_crm, res_naive_tite_pk, res_informed_tite_pk))
  }
}
