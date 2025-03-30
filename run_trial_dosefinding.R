##################################################################################################################
#### Configuring the architecture of each trial and running trials according to a selected dose-finding model ####
##################################################################################################################

# Core function
run_trial_dosefinding <- function(N, doses, ref_dose, cohort_size, real_sampling, PK_parameters, omega_IIV, omega_alpha, cv, simulated_PK_data, 
                                  exposure_tox_threshold, n_trials, method, targeted_tox, options, seed = 12345, proba_safety = 0.9) {
  
  ### Setting up the appropriate parameters to instantiate the method for correct trials running procedure
  AUCs <- NULL
  MTD <- NULL
  dose_levels <- NULL
  toxicity <- NULL
  p_estim_mean <- list()
  credible_intervals <- list()
  ref_log_AUC <- list()
  popPK_estimated_parameters <- list()
  seed_trial <- matrix(NA, nrow = 1, ncol = n_trials + 1)
  rownames(seed_trial) <- "Seed"
  colnames(seed_trial) <- c("Initial", paste("Trial ", 1:n_trials, sep = ""))
  seed_trial[1, 1] <- seed
  sel <- rep(0, length(doses) + 1)
  ka <- PK_parameters[1]
  Cl_pop <- PK_parameters[2]
  V_pop <- PK_parameters[3]
  
  if (method == "blrm") {
    model_BLRM <<- stan_model("stan/blrm.stan", "BLRM model")
  } else if(method == "pklogit") {
    # Hierarchical dose-AUC model by integrating parameter nu in the standard deviation of the normal distributions for beta0 and beta1
    model_PKLOGIT_AUC <<- stan_model("stan/pklogit_AUC_hierarchical.stan", "PK logit model - Hierarchical AUC")
    # Non hierarchical dose-AUC model by keeping the standard deviation of the normal distributions for beta0 and beta1 equal to 10
    model_PKLOGIT_logit <<- stan_model("stan/pklogit_logit.stan", "PK logit model - Toxicity")
  } else if(method == "ed_ewoc" | method == "ed") {
    model_BLRM <<- stan_model("stan/blrm.stan", "BLRM model")
    model_ED_EWOC_ED <<- stan_model("stan/ed_ewoc_ed_logit.stan", "ED-EWOC/ED logit model")
    popPK_estimated_parameters <- lapply(1:n_trials, matrix, data = NA, nrow = floor((N/2)/cohort_size), ncol = 5, dimnames = list(1:floor((N/2)/cohort_size), c("ka", "Cl_pop", "V_pop", "omega_Cl", "omega_V")))
  } else if(method == "pk_crm") {
    model_PK_CRM <<- stan_model("stan/pk_crm.stan", "PK-CRM model")
  } else if(method == "naive_tite_pk" | method == "informed_tite_pk") {
    ############################################
    ## Load the R files for TITE-PK method
    ###########################################
    source("stan/TITEPK_sequential-master/src/utils.R")
    source("stan/TITEPK_sequential-master/src/utils_tite_pk.R")
    source("stan/TITEPK_sequential-master/lib/tools.R")
    is_overdose_control <- TRUE # FALSE if no overdose control shall be used for the dose recommendation
    ############################################
    ## Compiling Stan files for TITE-PK method
    ###########################################
    stan_model_code_titepk <- stanc_builder("stan/TITEPK_sequential-master/src/tite_pk.stan")
    cat(stan_model_code_titepk$model_code, file="stan/TITEPK_sequential-master/src/tite_pk_run.stan")
    cat("\n", file="stan/TITEPK_sequential-master/src/tite_pk_run.stan", append = TRUE)
    ## compile Stan model functions via stanc_builder which avoids model
    ## name obfuscation and then admits caching in case nothing changes
    model_exposed_cached <- stanc_builder("stan/TITEPK_sequential-master/src/tite_pk_run.stan")
    cat("Exposing Stan functions...")
    expose_stan_functions(model_exposed_cached)
    cat("done.\n")
    ## Stan code for estimation
    model_TITE_PK <<- stan_model("stan/TITEPK_sequential-master/src/tite_pk_run.stan", "TITE-PK One parameter model")
    ############################################
    if (method == "naive_tite_pk") {
      ## NAIVE CHOICE of fixed pseudo-PK model parameters
      ## k_e: elimination rate constant / T_e: elimination half life
      ## Calculated from PK analysis (mean of the estimated T_e values)
      k_e <<- runif(1, 0, 50)/runif(1, 50, 200)
      T_e <<- log(2) / k_e
      ## k: kinetic constant which govern delay between concentration in central compart. and effect site
      ## calculated from PK analysis
      k <<- runif(1, 0, 10)
      naive_TITEPK_pseudoparameters <<- c("k" = k, "k_e" = k_e)
    } else {
      ## INFORMED CHOICE of fixed pseudo-PK model parameters
      ## k_e: elimination rate constant / T_e: elimination half life
      ## Calculated from PK analysis (mean of the estimated T_e values)
      k_e <<- Cl_pop/V_pop
      T_e <<- log(2) / k_e
      ## k: kinetic constant which govern delay between concentration in central compart. and effect site
      ## calculated from PK analysis
      k <<- ka
    }
    #########################################################
    ## Setting up the TITE-PK Regimens to be considered
    #########################################################
    ## This is similar to a NONMEM data structure
    ## Choosing the reference regimen
    ref_tau <<- 0
    ## Since cycle length is given as 1 day
    tref_month <<- 1/30.4167
    tau_drug <<- 0
    addl <<- 0
    # Calculating reference scale
    source("stan/TITEPK_sequential-master/src/pk_param.R")
    # Creating provisional regimen (Doses and frequencies)
    source("stan/TITEPK_sequential-master/src/regimen_ref.R")
    dosings <- regimen_ref_packed
    dosings <- dosings[rep(seq_len(nrow(dosings)), each=length(doses)),]
    dosings$amt <- doses
    dosings$lamt <- log(doses)
    regimens_daily <- dosings
    pk_calc <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
    # Calculating C(t_ref), E(t_ref), AUC_E(t_ref) for reference regimen
    pksys <<- ddply(dosings, .(amt), pk_subject, time = 1:tref_h, pk_fun = pk_calc)
    ## TITE-PK: Prior value for log(beta) parameters
    prior_titepk <- c(f_cloglog(0.25), 1.25)
  }
  
  # Executing the method for n_trials
  for (tr in 1:n_trials) {
    set.seed(seed + tr)
    conc <- simulated_PK_data[[tr]]$simu_exact_conc
    conc_with_error <- simulated_PK_data[[tr]]$simu_conc_with_error
    exact_auc <- simulated_PK_data[[tr]]$simu_exact_exposure
    sensitivity_auc <- simulated_PK_data[[tr]]$sensitivity_auc
    tox <- simulated_PK_data[[tr]]$toxicity
    Cl_ind <- simulated_PK_data[[tr]]$individual_parameters[2:(N+1)]
    V_ind <- simulated_PK_data[[tr]]$individual_parameters[(N+2):((2*N)+1)]
    alpha <- simulated_PK_data[[tr]]$alpha
    
    if (method == "naive_tite_pk" | method == "informed_tite_pk") {
      PK_data_titepk <- data.frame(
        cmt = numeric(),
        time = numeric(),
        amt = numeric(),
        week = numeric(),
        dv = numeric(),
        mdv = numeric(),
        evid = numeric(),
        lamt = numeric(),
        event_grade = numeric(),
        tau = numeric(),
        addl = numeric(),
        id = numeric(),
        cohort = numeric()
      )
    }
    
    # Initialization step - Stage 1 of the trial
    x <- rep(1, cohort_size)
    y <- tox[cbind(1:length(x), x)]
    M <- N / cohort_size
    j <- 1:cohort_size
    auc_s <- c()
    conci_all <- c()
    tr_ref_log_AUC <- c()
    for(k in j) {
      ### Required initialization step - Stage 1 of the trial - BLRM model
      if (method == "blrm") {
        conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
        conci_all <- c(conci_all, conci)
        auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
      } ### Required initialization step - Stage 1 of the trial - PKLOGIT model
      else if (method == "pklogit") {
        conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
        conci_all <- c(conci_all, conci)
        auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
      } ### Required initialization step - Stage 1 of the trial - ED-EWOC model or ED model
      else if (method == "ed_ewoc" | method == "ed"){
        conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
        conci_all <- c(conci_all, conci)
        auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
      } ### Required initialization step - Stage 1 of the trial - PK-CRM model
      else if (method == "pk_crm") {
        conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
        conci_all <- c(conci_all, conci)
        auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
      } ### Required initialization step - Stage 1 of the trial - TITE-PK model
      else if (method == "naive_tite_pk" | method == "informed_tite_pk") {
        if(y[k] == 1) { # patient has a DLT for the administered dose
          nb_lines <- 3 # number of lines per patient id in PK dataframe (with one line as the toxicity event)
          cmt <- c(1,10,11)
          time <- c(1,48,48)
          amt <- c(doses[x[k]],0,0)
          week <- rep(0,nb_lines)
          dv <- c(0,y[k],0)
          mdv <- c(1,0,0)
          evid <- c(1,0,0)
          lamt <- c(log(doses[x[k]]),-15,-15)
          event_grade <- rep(0, nb_lines)
          tau <- rep(tau_drug, nb_lines)
          addl_df <- rep(addl, nb_lines)
          id <- rep(k, nb_lines)
          cohort <- rep(1, nb_lines)
        } else { # patient has no DLT
          nb_lines <- 2 # number of lines per patient id in PK dataframe
          cmt <- c(1,11)
          time <- c(1,48)
          amt <- c(doses[x[k]],0)
          week <- rep(0,nb_lines)
          dv <- c(0,y[k])
          mdv <- c(1,0)
          evid <- c(1,0)
          lamt <- c(log(doses[x[k]]),-15)
          event_grade <- rep(0, nb_lines)
          tau <- rep(tau_drug, nb_lines)
          addl_df <- rep(addl, nb_lines)
          id <- rep(k, nb_lines)
          cohort <- rep(1, nb_lines)
        }
        ## Data needs to be in a list format for Stan
        PK_data_titepk <- rbind(PK_data_titepk,
                                data.frame(
                                  cmt = cmt,
                                  time = time,
                                  amt = amt,
                                  week = week,
                                  dv = dv,
                                  mdv = mdv,
                                  evid = evid,
                                  lamt = lamt,
                                  event_grade = event_grade,
                                  tau = tau,
                                  addl = addl_df,
                                  id = id,
                                  cohort = cohort
                                )
        )
        conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
        conci_all <- c(conci_all, conci)
        auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
      }
    }
    stage1 <- TRUE
    for (i in 2:M) {
      j <- (cohort_size * (i - 1) + 1):(cohort_size * i) # position
      ### Starting dose until first toxicity appears
      if (stage1) {
        x <- c(x, rep(min((max(x) + 1), length(doses)), cohort_size))
        y <- c(y, tox[cbind(j, x[j])])
        for (k in j) {
          if (method == "blrm") {
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if (method == "pklogit") {
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if (method == "ed_ewoc" | method == "ed"){
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if (method == "pk_crm"){
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if(method == "naive_tite_pk" | method == "informed_tite_pk") {
            if(y[k] == 1) { # Patient has a DLT for the administered dose
              nb_lines <- 3 # Number of lines per patient id in PK dataframe (with one line as the toxicity event)
              cmt <- c(1,10,11)
              time <- c(1,48,48)
              amt <- c(doses[x[k]],0,0)
              week <- rep(0,nb_lines)
              dv <- c(0,y[k],0)
              mdv <- c(1,0,0)
              evid <- c(1,0,0)
              lamt <- c(log(doses[x[k]]),-15,-15)
              event_grade <- rep(0, nb_lines)
              tau <- rep(tau_drug, nb_lines)
              addl_df <- rep(addl, nb_lines)
              id <- rep(k, nb_lines)
              cohort <- rep(i, nb_lines)
            } else { # Patient has no DLT
              nb_lines <- 2 # Number of lines per patient id in PK dataframe
              cmt <- c(1,11)
              time <- c(1,48)
              amt <- c(doses[x[k]],0)
              week <- rep(0,nb_lines)
              dv <- c(0,y[k])
              mdv <- c(1,0)
              evid <- c(1,0)
              lamt <- c(log(doses[x[k]]),-15)
              event_grade <- rep(0, nb_lines)
              tau <- rep(tau_drug, nb_lines)
              addl_df <- rep(addl, nb_lines)
              id <- rep(k, nb_lines)
              cohort <- rep(i, nb_lines)
            }
            ## Data needs to be in a list format for Stan
            PK_data_titepk <- rbind(PK_data_titepk,
                                    data.frame(
                                      cmt = cmt,
                                      time = time,
                                      amt = amt,
                                      week = week,
                                      dv = dv,
                                      mdv = mdv,
                                      evid = evid,
                                      lamt = lamt,
                                      event_grade = event_grade,
                                      tau = tau,
                                      addl = addl_df,
                                      id = id,
                                      cohort = cohort
                                    )
            )
            stan_data_titepk <- c(PK_data_titepk, 
                                  list(N = nrow(PK_data_titepk), 
                                       theta = log(c(T_e, k_e)),
                                       ref_dose = ref_dose,
                                       ref_tau = ref_tau,
                                       tref_month = tref_month,
                                       Nregimens = nrow(regimens_daily),
                                       regimens_lamt = regimens_daily$lamt,
                                       regimens_tau = regimens_daily$tau,
                                       addls = regimens_daily[, 7],
                                       params_prior = prior_titepk))
            
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          }
        }
        if (any(y == "1")) {
          stage1 <- FALSE
        }
      } ### After a first toxicity is declared, end stage 1 and proceed to model estimation
      else {
        # For Stan
        num <- length(x) # how many patients included so far
        if (method == "blrm") {
          res_model <- blrm_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, targeted_tox = targeted_tox,
                                  proba_safety = proba_safety, options = options, seed = seed + tr)
          new_dose <- res_model$new_dose
          p_estim_dose_mean <- res_model$p_estim_dose_mean
          credible_interval_per_dose <- res_model$credible_interval_per_dose
          
        } else if (method == "pklogit") {
          res_model <- pklogit_model(num = num, doses = doses, x = x, y = y, auc_s = auc_s, targeted_tox = targeted_tox,
                                     proba_safety = proba_safety, options = options, seed = seed + tr)
          new_dose <- res_model$new_dose
          p_estim_dose_mean <- res_model$p_estim_dose_mean
          credible_interval_per_dose <- res_model$credible_interval_per_dose
          tr_ref_log_AUC <- append(tr_ref_log_AUC, res_model$ref_log_AUC)
          
        } else if (method == "ed_ewoc") {
          res_model <- ed_ewoc_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, conci = conci_all, auc_s = auc_s, n_sampling_timepoints = n_sampling_timepoints, time_sampling = time_sampling, 
                                     real_sampling = real_sampling, PK_parameters = PK_parameters, targeted_tox = targeted_tox, proba_safety = proba_safety, options = options, seed = seed + tr)
          new_dose <- res_model$new_dose
          auc_s <- res_model$auc_s
          p_estim_dose_mean <- res_model$p_estim_dose_mean
          credible_interval_per_dose <- res_model$credible_interval_per_dose
          tr_ref_log_AUC <- append(tr_ref_log_AUC, res_model$ref_log_AUC)
          if (!is.null(res_model$estimated_PK_parameters)) {
            popPK_estimated_parameters[[tr]][(i-ceiling((N/2)/cohort_size)),] <- res_model$estimated_PK_parameters
          }
          
        } else if (method == "ed") {
          res_model <- ed_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, conci = conci_all, auc_s = auc_s, n_sampling_timepoints = n_sampling_timepoints, time_sampling = time_sampling, 
                                real_sampling = real_sampling, PK_parameters = PK_parameters, targeted_tox = targeted_tox, proba_safety = proba_safety, options = options, seed = seed + tr)
          new_dose <- res_model$new_dose
          auc_s <- res_model$auc_s
          p_estim_dose_mean <- res_model$p_estim_dose_mean
          credible_interval_per_dose <- res_model$credible_interval_per_dose
          tr_ref_log_AUC <- append(tr_ref_log_AUC, res_model$ref_log_AUC)
          if (!is.null(res_model$estimated_PK_parameters)) {
            popPK_estimated_parameters[[tr]][(i-ceiling((N/2)/cohort_size)),] <- res_model$estimated_PK_parameters
          }
          
        } else if (method == "pk_crm"){
          res_model <- pkcrm_model(num = num, doses = doses, x = x, y = y, auc_s = auc_s, targeted_tox = targeted_tox,
                                   proba_safety = proba_safety, options = options, seed = seed + tr)
          new_dose <- res_model$new_dose
          p_estim_dose_mean <- res_model$p_estim_dose_mean
          credible_interval_per_dose <- res_model$credible_interval_per_dose
          
        } else if (method == "naive_tite_pk" | method == "informed_tite_pk"){
          res_model <- titepk_model(num = num, doses = doses, stan_data_titepk = stan_data_titepk, targeted_tox = targeted_tox, proba_safety = proba_safety,
                                    is_overdose_control = is_overdose_control, options = options, seed = seed + tr)
          new_dose <- res_model$new_dose
          p_estim_dose_mean <- res_model$p_estim_dose_mean
          credible_interval_per_dose <- res_model$credible_interval_per_dose
        } 
        
        # Stop the trial if necessary
        if (any(is.na(new_dose)) == "TRUE") break
        
        # Assign the recommended dose to the next cohort if no dose has been missed in the escalation process;
        # otherwise, assign the higher dose to the previous one (dose skipping rule) within the limit of the dose pannel
        new_dose <- min(new_dose, max(x) + 1)
        if(new_dose > length(doses)) {
          new_dose <- new_dose - 1
        }
        
        # Add crucial information on AUC for the new cohort enrolled at the recommended dose for new iteration of the
        # dose escalation process
        x <- c(x, rep(new_dose, cohort_size))
        y <- c(y, tox[cbind(j, x[j])])
        for (k in j) {
          if (method == "blrm") {
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if (method == "pklogit") {
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if (method == "ed_ewoc" | method == "ed") {
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
          } else if (method == "pk_crm") {
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          } else if (method == "naive_tite_pk" | method == "informed_tite_pk") {
            if(y[k] == 1) { # patient has a DLT for the administered dose
              nb_lines <- 3 # number of lines per patient id in PK dataframe (with one line as the toxicity event)
              cmt <- c(1,10,11)
              time <- c(1,48,48)
              amt <- c(doses[x[k]],0,0)
              week <- rep(0,nb_lines)
              dv <- c(0,y[k],0)
              mdv <- c(1,0,0)
              evid <- c(1,0,0)
              lamt <- c(log(doses[x[k]]),-15,-15)
              event_grade <- rep(0, nb_lines)
              tau <- rep(tau_drug, nb_lines)
              addl_df <- rep(addl, nb_lines)
              id <- rep(k, nb_lines)
              cohort <- rep(i, nb_lines)
            } else { # patient has no DLT
              nb_lines <- 2 # number of lines per patient id in PK dataframe
              cmt <- c(1,11)
              time <- c(1,48)
              amt <- c(doses[x[k]],0)
              week <- rep(0,nb_lines)
              dv <- c(0,y[k])
              mdv <- c(1,0)
              evid <- c(1,0)
              lamt <- c(log(doses[x[k]]),-15)
              event_grade <- rep(0, nb_lines)
              tau <- rep(tau_drug, nb_lines)
              addl_df <- rep(addl, nb_lines)
              id <- rep(k, nb_lines)
              cohort <- rep(i, nb_lines)
            }
            ## Data needs to be in a list format for Stan
            PK_data_titepk <- rbind(PK_data_titepk,
                                    data.frame(
                                      cmt = cmt,
                                      time = time,
                                      amt = amt,
                                      week = week,
                                      dv = dv,
                                      mdv = mdv,
                                      evid = evid,
                                      lamt = lamt,
                                      event_grade = event_grade,
                                      tau = tau,
                                      addl = addl_df,
                                      id = id,
                                      cohort = cohort
                                    )
            )
            stan_data_titepk <- c(PK_data_titepk, 
                                  list(N = nrow(PK_data_titepk), 
                                       theta = log(c(T_e, k_e)),
                                       ref_dose = ref_dose,
                                       ref_tau = ref_tau,
                                       tref_month = tref_month,
                                       Nregimens = nrow(regimens_daily),
                                       regimens_lamt = regimens_daily$lamt,
                                       regimens_tau = regimens_daily$tau,
                                       addls = regimens_daily[, 7],
                                       params_prior = prior_titepk))
            conci <- as.vector(conc_with_error[[x[k]]][k, real_sampling])
            conci_all <- c(conci_all, conci)
            auc_s <- c(auc_s, PK::auc(conc = conci, time = time_sampling[real_sampling], design = "complete")$est)
          }
        }
      }
    }
    # Save important data at the end of the entire trial run in predefined lists
    trial <- paste0("Trial:", tr, sep = "")
    ## Check if the trial has been stopped before including all patients
    if (length(x) < N) {
      ### The trial has been stopped at some point
      n_stop <- N - length(x)
      MtD <- 0 # Final MTD after the trial
      mtd <- new_dose # Candidate MTD
      if (is.na(mtd) == "TRUE") {
        mtd <- 0
      }
      MTD <- c(MTD, mtd) # Vector of MTD for all trials
      sel[MtD + 1] <- sel[MtD + 1] + 1
      dose_levels <- rbind(dose_levels, c(x, rep(0, n_stop)))
      toxicity <- rbind(toxicity, c(y, rep(NA, n_stop)))
      AUCs <- rbind(AUCs, c(auc_s, rep(NA, n_stop)))
      ref_log_AUC[[trial]] <- tr_ref_log_AUC
      p_estim_mean[[trial]] <- p_estim_dose_mean
      credible_intervals[[trial]] <- credible_interval_per_dose
      seed_trial[1, tr + 1] <- seed + tr
    } else {
      ### The trial has been conducted until the end without any stopping
      # For Stan
      num <- length(x) # all patients included at the end of the trial
      if (method == "blrm") {
        res_model_final <- blrm_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, targeted_tox = targeted_tox, 
                                      proba_safety = proba_safety, options = options, seed = seed + tr)
        MtD <- res_model_final$new_dose
        p_estim_dose_mean <- res_model_final$p_estim_dose_mean
        credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        
      } else if (method == "pklogit") {
        res_model_final <- pklogit_model(num = num, doses = doses, x = x, y = y, auc_s = auc_s, targeted_tox = targeted_tox,
                                         proba_safety = proba_safety, options = options, seed = seed + tr)
        MtD <- res_model_final$new_dose
        p_estim_dose_mean <- res_model_final$p_estim_dose_mean
        credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        tr_ref_log_AUC <- append(tr_ref_log_AUC, res_model_final$ref_log_AUC)
        
      } else if (method == "ed_ewoc") {
        res_model_final <- ed_ewoc_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, conci = conci_all, auc_s = auc_s, n_sampling_timepoints = n_sampling_timepoints, time_sampling = time_sampling, 
                                         real_sampling = real_sampling, PK_parameters = PK_parameters, targeted_tox = targeted_tox, proba_safety = proba_safety, options = options, seed = seed + tr)
        MtD <- res_model_final$new_dose
        auc_s <- res_model_final$auc_s
        p_estim_dose_mean <- res_model_final$p_estim_dose_mean
        credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        tr_ref_log_AUC <- append(tr_ref_log_AUC, res_model_final$ref_log_AUC)
        popPK_estimated_parameters[[tr]][(i-ceiling((N/2)/cohort_size)),] <- res_model_final$estimated_PK_parameters
        
      } else if (method == "ed") {
        res_model_final <- ed_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, conci = conci_all, auc_s = auc_s, n_sampling_timepoints = n_sampling_timepoints, time_sampling = time_sampling, 
                                    real_sampling = real_sampling, PK_parameters = PK_parameters, targeted_tox = targeted_tox, proba_safety = proba_safety, options = options, seed = seed + tr)
        MtD <- res_model_final$new_dose
        auc_s <- res_model_final$auc_s
        p_estim_dose_mean <- res_model_final$p_estim_dose_mean
        credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        tr_ref_log_AUC <- append(tr_ref_log_AUC, res_model_final$ref_log_AUC)
        popPK_estimated_parameters[[tr]][(i-ceiling((N/2)/cohort_size)),] <- res_model_final$estimated_PK_parameters
        
      } else if (method == "pk_crm") {
        res_model_final <- pkcrm_model(num = num, doses = doses, x = x, y = y, auc_s = auc_s, targeted_tox = targeted_tox, 
                                       proba_safety = proba_safety, options = options, seed = seed + tr)
        MtD <- res_model_final$new_dose
        p_estim_dose_mean <- res_model_final$p_estim_dose_mean
        credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        
      } else if (method == "naive_tite_pk" | method == "informed_tite_pk"){
        if (all(y == "0")) {
          # If no DLTs are observed during a trial and the maximum sample size has been reached, a BLRM is used to estimate the probabilities of toxicity at the end.
          # This approach is necessary because TITE-PK encounters critical issues in the absence of observed DLTs!
          res_model_final <-  blrm_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, targeted_tox = targeted_tox, 
                                         proba_safety = proba_safety, options = options, seed = seed + tr)
          MtD <- res_model_final$new_dose
          p_estim_dose_mean <- res_model_final$p_estim_dose_mean
          credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        } else {
          res_model_final <- titepk_model(num = num, doses = doses, stan_data_titepk = stan_data_titepk, targeted_tox = targeted_tox, proba_safety = proba_safety,
                                          is_overdose_control = is_overdose_control, options = options, seed = seed + tr)
          MtD <- res_model_final$new_dose
          p_estim_dose_mean <- res_model_final$p_estim_dose_mean
          credible_interval_per_dose <- res_model_final$credible_interval_per_dose
        }
      } 
      
      # Declare that all doses are too toxic at the end of the trial after all cohorts are included
      if (any(is.na(MtD)) == "TRUE") {
        MtD <- 0 # Final MTD after the trial
        MTD <- c(MTD, MtD) # Vector of MTD for all trials
        sel[MtD + 1] <- sel[MtD + 1] + 1
        dose_levels <- rbind(dose_levels, x)
        toxicity <- rbind(toxicity, y)
        AUCs <- rbind(AUCs, auc_s)
        ref_log_AUC[[trial]] <- tr_ref_log_AUC
        p_estim_mean[[trial]] <- p_estim_dose_mean
        credible_intervals[[trial]] <- credible_interval_per_dose
        seed_trial[1, tr + 1] <- seed + tr
      } else {
        MTD <- c(MTD, MtD)
        sel[MtD + 1] <- sel[MtD + 1] + 1
        dose_levels <- rbind(dose_levels, x)
        toxicity <- rbind(toxicity, y)
        AUCs <- rbind(AUCs, auc_s)
        ref_log_AUC[[trial]] <- tr_ref_log_AUC
        p_estim_mean[[trial]] <- p_estim_dose_mean
        credible_intervals[[trial]] <- credible_interval_per_dose
        seed_trial[1, tr + 1] <- seed + tr
      }
    }
  }
  # Define patient id for "dose_finding_res" S4 class that will prompt the results of the simulations
  patient_id <- c(1:N)
  # Rename rows of matrices for dose levels, toxicity and AUCs obtained for each trial with the number of the trial
  rownames(dose_levels) <- sapply(1:n_trials, function(tr) paste0("Trial:", tr, sep = ""))
  rownames(toxicity) <- sapply(1:n_trials, function(tr) paste0("Trial:", tr, sep = ""))
  rownames(AUCs) <- sapply(1:n_trials, function(tr) paste0("Trial:", tr, sep = ""))
  # MTD if TR=1 otherwise the percentage of MTD selection for each dose level after all trials starting from dose 0;
  # equals to 0 if the trial has stopped before the end, according to the stopping rules.
  if (n_trials == 1) {
    new_dose <- MtD
  } else {
    new_dose <- sel / n_trials
  }
  
  
  # Output
  res <- new("dose_finding_res", patient_id = patient_id, N = N, 
             n_sampling_timepoints = n_sampling_timepoints, real_time_sampling = time_sampling[real_sampling],
             doses = doses, conc = conci, exposure_tox_threshold = exposure_tox_threshold, 
             n_chains = options$n_chains, n_iter = options$n_iter, n_adapt = options$n_adapt, 
             new_dose = new_dose, MTD_each_trial = MTD, MTD_final = MtD, 
             targeted_tox = targeted_tox, dose_levels = dose_levels, toxicity = toxicity, 
             AUCs = AUCs, ref_log_AUC = ref_log_AUC, n_trials = n_trials,
             prior_tox_real = prior_tox_real, p_estim = p_estim_mean,
             credible_intervals = credible_intervals,
             popPK_estimated_parameters = popPK_estimated_parameters,
             model = method, seed = seed_trial
  )
  return(res)
}