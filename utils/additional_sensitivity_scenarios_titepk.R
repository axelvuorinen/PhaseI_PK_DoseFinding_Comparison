# Evaluate the performance of the TITE-PK method using simulated toxicity time from a Weibull distribution
# instead of choosing to implement the toxicity event as a binary endpoint that can only occur
# at the end of the observation cycle as we did for the previous simulations.
library(ReIns)
seed <- 12345
set.seed(seed)

N <- 30
n_trials <- 1000
cohort_size <- 2
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3]
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
time_sampling <- seq(0, 24, length.out = 48)
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
exposure_metric <- "AUC"
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
targeted_tox <- 0.25
proba_safety <- 0.9
options <- list(n_chains = 4, n_iter = 8000, n_adapt = 0.8, n_cores = getOption("mc.cores", 1L))

# Simulated uniform/weibull latent toxicity times
lambda <- 0.8
rho <- 5
tox_times <- lapply(1:n_trials, function(i) 1+ReIns::rtweibull(N, shape = lambda, scale = rho, endpoint = 23))

# Scenario A1 (MTD = Dose 1)
scenario_id <- 1
exposure_tox_threshold <- 6.6
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
TITEPK_scenario_simulated_tox_time <- simu_data(
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
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)

# Scenario A2 (MTD = Dose 2)
scenario_id <- 2
exposure_tox_threshold <- 10.96
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
TITEPK_scenario_simulated_tox_time <- simu_data(
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
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)

# Scenario A3 (MTD = Dose 3)
scenario_id <- 3
exposure_tox_threshold <- 18.2
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
TITEPK_scenario_simulated_tox_time <- simu_data(
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
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)

# Scenario A4 (MTD = Dose 4)
scenario_id <- 4
exposure_tox_threshold <- 29.5
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
TITEPK_scenario_simulated_tox_time <- simu_data(
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
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)

# Scenario A5 (MTD = None: Stopping rule)
scenario_id <- 5
exposure_tox_threshold <- 3.9
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
TITEPK_scenario_simulated_tox_time <- simu_data(
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
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)



# Run the selected scenario among set A to compare the performance of the informed version of TITE-PK 
# having toxicity times to the one with only binary toxicity data
############################################
## Load the R files for TITE-PK method
###########################################
source("stan/TITEPK_sequential-master/src/utils.R")
source("stan/TITEPK_sequential-master/src/utils_tite_pk.R")
source("stan/TITEPK_sequential-master/lib/tools.R")
is_overdose_control <- FALSE
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
## INFORMED CHOICE of fixed pseudo-PK model parameters
## T_e: half elimination rate constant
## Calculated from PK analysis (mean of the estimated T_e values)
k_e <<- Cl_pop/V_pop
## k_eff: kinetic constant which govern delay betwenn concentration in central compart. and effect site
## calculated from PK analysis
k <<- ka
T_e <<- log(2) / k
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
prior_titepk <- c(f_cloglog(0.25), 1.25) #cloglog(0.25)

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

for (tr in 1:n_trials) {
  set.seed(seed + tr)
  conc <- TITEPK_scenario_simulated_tox_time[[tr]]$simu_exact_conc
  conc_with_error <- TITEPK_scenario_simulated_tox_time[[tr]]$simu_conc_with_error
  exact_auc <- TITEPK_scenario_simulated_tox_time[[tr]]$simu_exact_exposure
  sensitivity_auc <- TITEPK_scenario_simulated_tox_time[[tr]]$sensitivity_auc
  tox <- TITEPK_scenario_simulated_tox_time[[tr]]$toxicity
  Cl_ind <- TITEPK_scenario_simulated_tox_time[[tr]]$individual_parameters[2:(N+1)]
  V_ind <- TITEPK_scenario_simulated_tox_time[[tr]]$individual_parameters[(N+2):((2*N)+1)]
  alpha <- TITEPK_scenario_simulated_tox_time[[tr]]$alpha
  toxicity_time_per_patient <- rep(NA, N)
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
  
  # Initialization step - Stage 1 of the trial
  x <- rep(1, cohort_size)
  y <- tox[cbind(1:length(x), x)]
  M <- N / cohort_size
  j <- 1:cohort_size
  toxicity_time_per_patient[j] <- sapply(j, function(z) ifelse(y[z] == 1, tox_times[[tr]][z], NA))
  auc_s <- c()
  conci_all <- c()
  tr_ref_log_AUC <- c()
  
  for(k in j) {
    ### Required initialization step - Stage 1 of the trial - TITE-PK model
    if(y[k] == 1) { # patient has a DLT for the administered dose
      nb_lines <- 3 # number of lines per patient id in PK dataframe (with one line as the toxicity event)
      cmt <- c(1,10,11)
      time <- c(1,toxicity_time_per_patient[k],24)
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
  stage1 <- TRUE
  for (i in 2:M) {
    j <- (cohort_size * (i - 1) + 1):(cohort_size * i) # position
    ### Starting dose until first toxicity appears
    if (stage1) {
      x <- c(x, rep(min((max(x) + 1), length(doses)), cohort_size))
      y <- c(y, tox[cbind(j, x[j])])
      toxicity_time_per_patient[j] <- sapply(j, function(z) ifelse(y[z] == 1, tox_times[[tr]][z], NA))
      for (k in j) {
        if(y[k] == 1) { # Patient has a DLT for the administered dose
          nb_lines <- 3 # Number of lines per patient id in PK dataframe (with one line as the toxicity event)
          cmt <- c(1,10,11)
          time <- c(1,toxicity_time_per_patient[k],24)
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
          time <- c(1,24)
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
      if (any(y == "1")) {
        stage1 <- FALSE
      }
    } ### After a first toxicity is declared, end stage 1 and proceed to model estimation
    else {
      # For Stan
      num <- length(x) # how many patients included so far
      res_model <- titepk_model(num = num, doses = doses, stan_data_titepk = stan_data_titepk, targeted_tox = targeted_tox, proba_safety = proba_safety,
                                is_overdose_control = is_overdose_control, options = options, seed = seed + tr)
      new_dose <- res_model$new_dose
      p_estim_dose_mean <- res_model$p_estim_dose_mean
      credible_interval_per_dose <- res_model$credible_interval_per_dose
      
      # Stop the trial if necessary
      if (any(is.na(new_dose)) == "TRUE") break
      
      # Assign the recommended dose to the next cohort if no dose has been missed in the escalation process;
      # otherwise, assign the higher dose to the previous one (dose skipping rule) within the limit of the dose panel
      new_dose <- min(new_dose, max(x) + 1)
      if(new_dose > length(doses)) {
        new_dose <- new_dose - 1
      }
      
      # Add crucial information on AUC for the new cohort enrolled at the recommended dose for new iteration of the
      # dose escalation process
      x <- c(x, rep(new_dose, cohort_size))
      y <- c(y, tox[cbind(j, x[j])])
      toxicity_time_per_patient[j] <- sapply(j, function(z) ifelse(y[z] == 1, tox_times[[tr]][z], NA))
      for (k in j) {
        if(y[k] == 1) { # patient has a DLT for the administered dose
          nb_lines <- 3 # number of lines per patient id in PK dataframe (with one line as the toxicity event)
          cmt <- c(1,10,11)
          time <- c(1,toxicity_time_per_patient[k],24)
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
          time <- c(1,24)
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
    if (all(y == "0")) {
      # If no DLTs are observed during a trial and the maximum sample size has been reached, a BLRM is used to estimate the probabilities of toxicity at the end.
      # This approach is necessary because TITE-PK encounters critical issues in the absence of observed DLTs!
      res_model_final <<-  blrm_model(num = num, doses = doses, ref_dose = ref_dose, x = x, y = y, targeted_tox = targeted_tox, 
                                      proba_safety = proba_safety, options = options, seed = seed + tr)
      MtD <- res_model_final$new_dose
      p_estim_dose_mean <- res_model_final$p_estim_dose_mean
      credible_interval_per_dose <- res_model_final$credible_interval_per_dose
    } else {
      res_model_final <<- titepk_model(num = num, doses = doses, stan_data_titepk = stan_data_titepk, targeted_tox = targeted_tox, proba_safety = proba_safety,
                                       is_overdose_control = is_overdose_control, options = options, seed = seed + tr)
      MtD <- res_model_final$new_dose
      p_estim_dose_mean <- res_model_final$p_estim_dose_mean
      credible_interval_per_dose <- res_model_final$credible_interval_per_dose
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
           model = "weibull_informed_tite_pk", seed = seed_trial
)