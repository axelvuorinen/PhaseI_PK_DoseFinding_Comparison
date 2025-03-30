##############################################################################
################### Creating simulation scenarios ############################
##############################################################################

# In case, there is any bug in the execution of the simulations that can be attributed to rstan
# and the asssociated dependencies used by rstan, execute the line under
# install.packages(c("BH", "StanHeaders", "Rcpp", "RcppEigen", "RcppParallel", "inline", "loo", "pkgbuild", "rstan"))

### File settings
rm(list = ls())
folder_path <- "..." # Folder path 
setwd(folder_path)
# Initialization seed for each scenario and each method
seed <- 12345
set.seed(seed)



### Loading libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(ggridges)
library(patchwork)
library(latex2exp)
library(PK)
library(cubature)
library(deSolve)
library(rstan)
library(saemix)
library(dplyr)
library(functional)
library(loo)
library(assertthat)





### Manually import all necessary R files with their functions into the R environment to run the simulations
source(paste0(folder_path, "run_sim_scenario.R"))
source(paste0(folder_path, "run_trial_dosefinding.R"))
source(paste0(folder_path, "simu_data.R"))
source(paste0(folder_path, "dosefinding_method/blrm.R"))
source(paste0(folder_path, "dosefinding_method/pklogit.R"))
source(paste0(folder_path, "dosefinding_method/ed_ewoc.R"))
source(paste0(folder_path, "dosefinding_method/ed.R"))
source(paste0(folder_path, "dosefinding_method/pkcrm.R"))
source(paste0(folder_path, "dosefinding_method/titepk.R"))
source(paste0(folder_path, "utils/basic_functions.R"))
source(paste0(folder_path, "utils/S4_classes.R"))



### Setting up the fundamental simulation parameters that will remained unchanged across all scenarios
# Number of patients per trial
N <- 30
# Number of trials to be simulated
n_trials <- 1000
# Toxicity threshold at 25% to determine the recommended dose or/and the MTD at each trial
targeted_tox <- 0.25
# Timing PK sampling assessment for PK data simulation
time_sampling <- seq(0, 24, length.out = 48)
# Rstan parameters for Bayesian inference and maximum tolerated probability of toxicity for safety monitoring
options <- list(n_chains = 4, n_iter = 8000, n_adapt = 0.8, n_cores = getOption("mc.cores", 1L))
proba_safety <- 0.9



########################################################################################################
### Scenario A1 - Scenario with a one-compartment PK model and AUC-estimated toxicity - MTD = Dose 1 ###
########################################################################################################

scenario_id <- 1

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC"
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 6.6 #New = 6.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_A1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_A1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_A1.Rdata"
)



########################################################################################################
### Scenario A2 - Scenario with a one-compartment PK model and AUC-estimated toxicity - MTD = Dose 2 ###
########################################################################################################

scenario_id <- 2

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 10.96 #New = 10.3
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_A2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_A2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_A2.Rdata"
)



####################################################################################################################
### Scenario A3 (REFERENCE) - Scenario with a one-compartment PK model and AUC-estimated toxicity - MTD = Dose 3 ###
####################################################################################################################

scenario_id <- 3

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 18.2 #New = 18.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_A3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_A3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_A3.Rdata"
)



########################################################################################################
### Scenario A4 - Scenario with a one-compartment PK model and AUC-estimated toxicity - MTD = Dose 4 ###
########################################################################################################

scenario_id <- 4

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 29.5 #New = 30.4
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_A4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_A4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_A4.Rdata"
)



#####################################################################################################################
### Scenario A5 - Scenario with a one-compartment PK model and AUC-estimated toxicity - MTD = None: Stopping rule ###
#####################################################################################################################

scenario_id <- 5

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 3.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_A5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_A5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_A5.Rdata"
)


########################################################################################################
### Scenario B1 - Scenario with a two-compartment PK model and AUC-estimated toxicity - MTD = Dose 1 ###
########################################################################################################

scenario_id <- 6

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 6.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_2cmt(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_B1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_B1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_B1.Rdata"
)



########################################################################################################
### Scenario B2 - Scenario with a two-compartment PK model and AUC-estimated toxicity - MTD = Dose 2 ###
########################################################################################################

scenario_id <- 7

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC"
AUC_method <- "integration"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 10.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_2cmt(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_B2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_B2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_B2.Rdata"
)



########################################################################################################
### Scenario B3 - Scenario with a two-compartment PK model and AUC-estimated toxicity - MTD = Dose 3 ###
########################################################################################################

scenario_id <- 8

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 18.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_2cmt(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_B3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_B3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_B3.Rdata"
)



########################################################################################################
### Scenario B4 - Scenario with a two-compartment PK model and AUC-estimated toxicity - MTD = Dose 4 ###
########################################################################################################

scenario_id <- 9

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 30.3
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_2cmt(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_B4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_B4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_B4.Rdata"
)



#####################################################################################################################
### Scenario B5 - Scenario with a two-compartment PK model and AUC-estimated toxicity - MTD = None: Stopping rule ###
#####################################################################################################################

scenario_id <- 10

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 3.8
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_2cmt(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_B5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_B5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_B5.Rdata"
)



#########################################################################################################
### Scenario C1 - Scenario with a one-compartment PK model and Cmax-estimated toxicity - MTD = Dose 1 ###
#########################################################################################################

scenario_id <- 11

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 0.49
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                       exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                       time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_C1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_C1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_C1.Rdata"
)



#########################################################################################################
### Scenario C2 - Scenario with a one-compartment PK model and Cmax-estimated toxicity - MTD = Dose 2 ###
#########################################################################################################

scenario_id <- 12

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 0.82
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_C2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_C2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_C2.Rdata"
)


#########################################################################################################
### Scenario C3 - Scenario with a one-compartment PK model and Cmax-estimated toxicity - MTD = Dose 3 ###
#########################################################################################################

scenario_id <- 13

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 1.51
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_C3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_C3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_C3.Rdata"
)


#########################################################################################################
### Scenario C4 - Scenario with a one-compartment PK model and Cmax-estimated toxicity - MTD = Dose 4 ###
#########################################################################################################

scenario_id <- 14

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 2.43
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_C4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_C4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_C4.Rdata"
)



######################################################################################################################
### Scenario C5 - Scenario with a one-compartment PK model and Cmax-estimated toxicity - MTD = None: Stopping rule ###
######################################################################################################################

scenario_id <- 15

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 0.31
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_C5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_C5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_C5.Rdata"
)



#########################################################################################################
### Scenario D1 - Scenario with a two-compartment PK model and Cmax-estimated toxicity - MTD = Dose 1 ###
#########################################################################################################

scenario_id <- 16

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 0.45
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_D1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_D1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_D1.Rdata"
)



##################################################################################################################
### Scenario D2 - Standard scenario with a two-compartment PK model and Cmax-estimated toxicity - MTD = Dose 2 ###
##################################################################################################################

scenario_id <- 17

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 0.75
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                                          exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                                          time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_D2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_D2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_D2.Rdata"
)



#########################################################################################################
### Scenario D3 - Scenario with a two-compartment PK model and Cmax-estimated toxicity - MTD = Dose 3 ###
#########################################################################################################

scenario_id <- 18

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 1.38
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_D3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_D3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_D3.Rdata"
)


#########################################################################################################
### Scenario D4 - Scenario with a two-compartment PK model and Cmax-estimated toxicity - MTD = Dose 4 ###
#########################################################################################################

scenario_id <- 19

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 2.22
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_D4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_D4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_D4.Rdata"
)


######################################################################################################################
### Scenario D5 - Scenario with a two-compartment PK model and Cmax-estimated toxicity - MTD = None: Stopping rule ###
######################################################################################################################

scenario_id <- 20

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's Cmax for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold Cmax - if Cmax_alpha is above this threshold, the DLT is declared.
exposure_metric <- "Cmax"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 0.29
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_Cmax(PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                        exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                        time_sampling = time_sampling, compartmental_model = compartmental_model)


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_D5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_D5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_D5.Rdata"
)



####################################################################################################################################################
### Scenario E1 - Scenario with a one-compartment PK model using Michaelis-Menten elimination kinetics and AUC-estimated toxicity - MTD = Dose 1 ###
####################################################################################################################################################

scenario_id <- 21

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
Vmax_pop <- 10
Km_pop <- 2
PK_parameters <- c(ka, Cl_pop, V_pop, Vmax_pop, Km_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "Michaelis-Menten"
exposure_tox_threshold <- 11.5
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_MM(PK_parameters = PK_parameters, omega_IIV = , omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                     exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                     time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_E1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_E1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_E1.Rdata"
)



####################################################################################################################################################
### Scenario E2 - Scenario with a one-compartment PK model using Michaelis-Menten elimination kinetics and AUC-estimated toxicity - MTD = Dose 2 ###
####################################################################################################################################################


scenario_id <- 22

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
Vmax_pop <- 10
Km_pop <- 2
PK_parameters <- c(ka, Cl_pop, V_pop, Vmax_pop, Km_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "Michaelis-Menten"
exposure_tox_threshold <- 20
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_MM(PK_parameters = PK_parameters, omega_IIV = , omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                     exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                     time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_E2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_E2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_E2.Rdata"
)



####################################################################################################################################################
### Scenario E3 - Scenario with a one-compartment PK model using Michaelis-Menten elimination kinetics and AUC-estimated toxicity - MTD = Dose 3 ###
####################################################################################################################################################

scenario_id <- 23

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
Vmax_pop <- 10
Km_pop <- 2
PK_parameters <- c(ka, Cl_pop, V_pop, Vmax_pop, Km_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "Michaelis-Menten"
exposure_tox_threshold <- 41
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_MM(PK_parameters = PK_parameters, omega_IIV = , omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                     exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                     time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_E3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_E3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_E3.Rdata"
)



####################################################################################################################################################
### Scenario E4 - Scenario with a one-compartment PK model using Michaelis-Menten elimination kinetics and AUC-estimated toxicity - MTD = Dose 4 ###
####################################################################################################################################################

scenario_id <- 24

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
Vmax_pop <- 10
Km_pop <- 2
PK_parameters <- c(ka, Cl_pop, V_pop, Vmax_pop, Km_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "Michaelis-Menten"
exposure_tox_threshold <- 75
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_MM(PK_parameters = PK_parameters, omega_IIV = , omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                     exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                     time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_E4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_E4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_E4.Rdata"
)



#################################################################################################################################################################
### Scenario E5 - Scenario with a one-compartment PK model using Michaelis-Menten elimination kinetics and AUC-estimated toxicity - MTD = None: Stopping rule ###
#################################################################################################################################################################

scenario_id <- 25

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
Vmax_pop <- 10
Km_pop <- 2
PK_parameters <- c(ka, Cl_pop, V_pop, Vmax_pop, Km_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "Michaelis-Menten"
exposure_tox_threshold <- 8
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- f_prior_tox_AUC_MM(PK_parameters = PK_parameters, omega_IIV = , omega_IIV, omega_alpha = omega_alpha, cv = cv, doses = doses,
                                     exposure_tox_threshold = exposure_tox_threshold, n_sampling_timepoints = n_sampling_timepoints,
                                     time_sampling = time_sampling, compartmental_model = compartmental_model)

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_E5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_E5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_E5.Rdata"
)



###############################################################################################################################################################
### Scenario F - Scenario with a one-compartment PK model using AUC-estimated toxicity and with an optimal reduced number of sampling times - MTD = Dose 3 ###
###############################################################################################################################################################

scenario_id <- 26

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 4
real_sampling <- c(2, 3, 17, 48) # Sampling at 30 minutes, 1 hour, 8 hours and 24 hours.
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC"
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 18.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))


### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_F <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_F", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_F.Rdata"
)


##############################################################################################
### Scenario G1 - Standard scenario but with a lower IIV around 20% (< 70%) - MTD = Dose 1 ###
##############################################################################################

scenario_id <- 27

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.2 # Lower IIV (0.2 <<< 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 5.3
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_G1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_G1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_G1.Rdata"
)


##############################################################################################
### Scenario G2 - Standard scenario but with a lower IIV around 20% (< 70%) - MTD = Dose 2 ###
##############################################################################################

scenario_id <- 28

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.2 # Lower IIV (0.2 <<< 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 8.8
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_G2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_G2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_G2.Rdata"
)



##############################################################################################
### Scenario G3 - Standard scenario but with a lower IIV around 20% (< 70%) - MTD = Dose 3 ###
##############################################################################################

scenario_id <- 29

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.2 # Lower IIV (0.2 <<< 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 16.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_G3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                           PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                           real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                           exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                           options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_G3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_G3.Rdata"
)



##############################################################################################
### Scenario G4 - Standard scenario but with a lower IIV around 20% (< 70%) - MTD = Dose 4 ###
##############################################################################################

scenario_id <- 30

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.2 # Lower IIV (0.2 <<< 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 25.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_G4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_G4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_G4.Rdata"
)



###########################################################################################################
### Scenario G5 - Standard scenario but with a lower IIV around 20% (< 70%) - MTD = None: Stopping rule ###
###########################################################################################################

scenario_id <- 31

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.2 # Lower IIV (0.2 <<< 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 3.6
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_G5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_G5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_G5.Rdata"
)



#################################################################################################
### Scenario H1 - Standard scenario but with a higher IIV around 120% (> 100%) - MTD = Dose 1 ###
#################################################################################################

scenario_id <- 32

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 1.2 # Higher IIV (1.2 >>> 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 7.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_H1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_H1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_H1.Rdata"
)



#################################################################################################
### Scenario H2 - Standard scenario but with a higher IIV around 120% (> 100%) - MTD = Dose 2 ###
#################################################################################################

scenario_id <- 33

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 1.2 # Higher IIV (1.2 >>> 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 13.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_H2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_H2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_H2.Rdata"
)



#################################################################################################
### Scenario H3 - Standard scenario but with a higher IIV around 120% (> 100%) - MTD = Dose 3 ###
#################################################################################################

scenario_id <- 34

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 1.2 # Higher IIV (1.2 >>> 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 24.4
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_H3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_H3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_H3.Rdata"
)



#################################################################################################
### Scenario H4 - Standard scenario but with a higher IIV around 120% (> 100%) - MTD = Dose 4 ###
#################################################################################################

scenario_id <- 35

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 1.2 # Higher IIV (1.2 >>> 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 39.1
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_H4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_H4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_H4.Rdata"
)


##############################################################################################################
### Scenario H5 - Standard scenario but with a higher IIV around 120% (> 100%) - MTD = None: Stopping rule ###
##############################################################################################################

scenario_id <- 36

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 1.2 # Higher IIV (1.2 >>> 0.7)
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 4.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_H5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_H5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_H5.Rdata"
)


####################################################################################################################################################
### Scenario I1 - Add an additional high dose to the explored dose panel to improve the modeling range of the dose-response curve - MTD = Dose 1 ###
####################################################################################################################################################

scenario_id <- 37

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37, 220.65)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 6.2
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_I1 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_I1", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_I1.Rdata"
)



####################################################################################################################################################
### Scenario I2 - Add an additional high dose to the explored dose panel to improve the modeling range of the dose-response curve - MTD = Dose 2 ###
####################################################################################################################################################

scenario_id <- 38

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37, 220.65)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 10.3
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_I2 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_I2", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_I2.Rdata"
)



####################################################################################################################################################
### Scenario I3 - Add an additional high dose to the explored dose panel to improve the modeling range of the dose-response curve - MTD = Dose 3 ###
####################################################################################################################################################

scenario_id <- 39

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37, 220.65)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 18.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_I3 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_I3", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_I3.Rdata"
)



####################################################################################################################################################
### Scenario I4 - Add an additional high dose to the explored dose panel to improve the modeling range of the dose-response curve - MTD = Dose 4 ###
####################################################################################################################################################

scenario_id <- 40

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37, 220.65)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 30.4
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_I4 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_I4", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_I4.Rdata"
)



####################################################################################################################################################
### Scenario I5 - Add an additional high dose to the explored dose panel to improve the modeling range of the dose-response curve - MTD = Dose 5 ###
####################################################################################################################################################

scenario_id <- 41

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37, 220.65)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 44.6
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_I5 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_I5", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_I5.Rdata"
)



#################################################################################################################################################################
### Scenario I6 - Add an additional high dose to the explored dose panel to improve the modeling range of the dose-response curve - MTD = None: Stopping rule ###
#################################################################################################################################################################

scenario_id <- 42

### Setting up the fundamental parameters of all simulated trials (doses panel, reference dose, etc.)
doses <- c(30.6, 50.69, 93.69, 150.37, 220.65)
ref_dose <- doses[3] # Reference dose chosen for all methods and tests in the current scenario
# Cohort size for patient inclusion
cohort_size <- 2
### Setting up the PK parameters of the scenario for all simulated trial
# First-order absorption linear one compartment model parameters' (ka, Cl and V) with standard deviation
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7 # 0.7 for high IIV
# Standard deviation of an alpha parameter related to the sensitivity of the subject's AUC for determining toxicity data
omega_alpha <- 0.8 # 1.17 for more discrepancies between subject-specific threshold
# Timing and number of PK sampling assessment
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
# Proportional error model of 20% for computed concentration with a measure error
cv <- 0.20
# Biological toxicity threshold AUC - if AUC_alpha is above this threshold, the DLT is declared.
exposure_metric <- "AUC" 
AUC_method <- "compartmental"
elimination_kinetics <- "linear"
exposure_tox_threshold <- 3.9
# Prior probability of toxicity as a function of dose and fixed pharmacokinetic parameters
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))

### Plugging all the scenario parameters into the scenario execution function
system.time(scenario_I6 <- run_sim_scenario(N = N, doses = doses, ref_dose = ref_dose, n_trials = n_trials, cohort_size = cohort_size, targeted_tox = targeted_tox,
                                            PK_parameters = PK_parameters, omega_IIV = omega_IIV, omega_alpha = omega_alpha, n_sampling_timepoints = n_sampling_timepoints,
                                            real_sampling = real_sampling, cv = cv, compartmental_model = compartmental_model, elimination_kinetics = elimination_kinetics,
                                            exposure_metric = exposure_metric, AUC_method = AUC_method, exposure_tox_threshold = exposure_tox_threshold,
                                            options = options, seed = seed, proba_safety = proba_safety, scenario_id = scenario_id))

### Saving the most important simulation output for later confirmatory analysis in the case of 
### inconsistent results 
save(
  list = c("scenario_I6", "simulated_PK_data", "res_blrm", "res_pklogit", "res_ed_ewoc", "res_ed", 
           "res_pk_crm", "res_naive_tite_pk", "res_informed_tite_pk", "post_proba_tox_plot", 
           "post_proba_tox_CI_plot", "dose_selection_scenario_barplot", "dose_allocation_scenario_barplot",
           "ED_EWOC_popPK_plot", "ED_popPK_plot", "show"), 
  file = "run_scen/scenario_I6.Rdata"
)