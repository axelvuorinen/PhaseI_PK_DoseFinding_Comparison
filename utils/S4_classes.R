setClassUnion("class_new_dose", c("numeric", "logical", "NULL"))

setClass("dose_finding_res", slots = list(patient_id = "numeric", N = "numeric", n_sampling_timepoints = "numeric", real_time_sampling = "numeric", doses = "numeric", conc = "numeric", 
                                          exposure_tox_threshold = "numeric",  n_chains = "numeric", n_iter = "numeric", n_adapt = "numeric", new_dose = "class_new_dose", 
                                          MTD_each_trial = "class_new_dose", MTD_final = "class_new_dose", targeted_tox = "numeric", dose_levels = "matrix", 
                                          toxicity = "matrix", AUCs = "matrix", ref_log_AUC = "list", n_trials = "numeric", prior_tox_real = "vector", p_estim = "list", 
                                          credible_intervals = "list", popPK_estimated_parameters = "list", model = "character", seed = "matrix"))

#' An S4 class to represent a dose finding results.
#'
#' @slot patient_id The patient's ID provided in the study.
#' @slot N  The total sample size per trial.
#' @slot n_sampling_timepoints  The number of concentration sampling time points obtained during each trial.
#' @slot real_time_sampling The sampling time points.
#' @slot doses A vector with the doses panel.
#' @slot conc The estimated concentration values for each patient at each dose.
#' @slot exposure_tox_threshold  The AUC threshold or Cmax threshold to be set before starting the trial to determine the toxicity outcome.
#' @slot n_chains The number of chains for the Stan model.
#' @slot n_iter The number of iterations for the Stan model.
#' @slot n_adapt  The number of warmup iterations for the Stan model.
#' @slot new_dose  The next maximum tolerated dose (MTD) if TR=1 otherwise the percentage of MTD selection for each dose level after all trials starting from dose 0; equals to 0 if the trial has stopped before the end, according to the stopping rules.
#' @slot MTD_each_trial A vector containing the next maximum tolerated doses (MTD) of each trial (TR); equals to 0 if the trial has stopped before the end, according to the stopping rules.
#' @slot MTD_final The final next maximum tolerated (MTD) dose after all the trials.
#' @slot targeted_tox  The toxicity target.
#' @slot dose_levels A vector of dose levels assigned to patients in the trial.
#' @slot toxicity The estimated toxicity outcome.
#' @slot AUCs A vector with the computed (estimated) AUC values of each patient based on concentrations with error.
#' @slot ref_log_AUC The estimated reference log(AUC) used for normalization in the computation of AUC-toxicity model for some PK dose-finding designs (PKLOGIT, ED-EWOC/ED, PK-CRM)
#' @slot n_trials The total number of trials to be simulated.
#' @slot prior_tox_real The prior toxicity for each dose computed using equation (13) of Ursino (2017) based on the parameters of the simulated scenario.
#' @slot p_estim The estimated mean probabilities of toxicity.
#' @slot credible_intervals 95% credible intervals for the mean posterior probability of toxicity for each dose.
#' @slot popPK_estimated_parameters Only for ED-EWOC and ED methods - Estimated PK parameters by pop-PK model for each trial and each cohort of patients included.
#' @slot model A character string to specify the selected dose-finding model.
#' @slot seed The seed of the random number generator that is used during the initialization phase and at the beginning of each trial.
#@slot saved_warnings = vector saved_warnings All warnings obtained from iterative rstan (Bayesian) adaptive model to include all cohort of patients for all trials.
#' @import methods
#' @export


#' An S4 class to get a better printing format for the above S4 class on dose finding results.
setGeneric("show")
#' @export 
setMethod(f = "show", signature = "dose_finding_res", definition = function(object)
{
  cat("Today: ", date(), "\n") 
  cat("\n","A. Data Summary (", object@model,"model)", "\n")
  cat("Number of simulations:", object@n_trials, "\n")
  cat("Total number of patients in the trial:", object@N, "\n")
  cat("Number of concentration sampling time points:", object@n_sampling_timepoints, "\n")
  cat("Sampling time after the drug administration (hours):", round(object@real_time_sampling, digits = 3), "\n")
  n <- object@N
  cat("Levels of Doses:", round(object@doses, digits=3), "\n")
  cat("Concentration of the drug:", round(object@conc, digits = 3), "\n")
  cat("Initialization seed:",object@seed[1], "\n")
  cat("\n","B. STAN Model's Options \n")
  cat("The Stan model runs with", object@n_chains, "MCMC chains ")
  cat("which each chain has", object@n_iter, "iterations ")
  cat("and", object@n_adapt, "warmup iterations \n")
  if(object@n_trials == "1"){
    cat("\n","C. Dose-Finding Results: \n")
    cat("PID", "\t", "Level", "\t", "Toxicity", "\t", "AUCs", "\n")
    for (i in 1:n){
      cat(i,"\t", object@dose_levels[i],"\t", object@toxicity[i],"\t", "\t", round(object@AUCs[i], digits=3) ,"\n")
    }
    cat("\nThe prior probabilities of toxicity:", round(object@prior_tox_real, digits=4), "\n")
    cat("Next recommended dose level:", object@new_dose, "\n")
    cat("Recommendation is based on a target toxicity probability of:",object@targeted_tox, "\n")
    cat("Initial seed and seed of the trial:",object@seed, "\n")
  }else{
    cat("\n","C. Dose-Finding Results: \n")
    doselevels <- as.vector(object@dose_levels)
    t <- matrix(NA, nrow=4, ncol=length(object@doses)+1)
    rownames(t) <- c("Dose", "Truth Probabilities", "Dose-Allocation (%)", "Selected % MTD")
    colnames(t) <- rep("", length(object@doses)+1)
    stop <- as.character("STOP")
    dd <- paste("", 1:length(object@doses), sep = "")
    t[1, ] <- c(stop, dd)
    t[2, ] <- c("NA", round(object@prior_tox_real, digits=3))
    t[3,1] <- "NA"
    for(i in 1:length(object@doses)){
      doselevels_withoutSTOP <- doselevels
      doselevels_withoutSTOP[doselevels_withoutSTOP == 0] <- NA
      n_levels <- length(which(doselevels_withoutSTOP == i))
      t[3, i+1] <- round(n_levels/length(na.omit(doselevels_withoutSTOP)), digits=3)
    }
    t[4, ] <- round(object@new_dose, digits=3)
    print(t, quote = FALSE)
    cat("Recommendation is based on a target toxicity probability of:",object@targeted_tox, "\n")
    cat("Seed of each trial:",object@seed[-1], "\n")
  }
}
)




setClass("simulation_res", slots = list(scenario = "numeric", PK_parameters = "vector", omega_IIV = "numeric", omega_alpha = "numeric",
                                        cv = "numeric", exposure_metric = "character", prior_tox_real = "vector", N = "numeric", 
                                        cohort_size = "numeric", n_trials = "numeric", n_sampling_timepoints = "numeric", 
                                        real_time_sampling = "numeric", doses = "numeric", exposure_tox_threshold = "numeric",  n_chains = "numeric", 
                                        n_iter = "numeric", n_adapt = "numeric", n_cores = "numeric", dose_allocation = "matrix", 
                                        dose_selection = "matrix", targeted_tox = "numeric", median_DLTs = "vector", min_DLTs = "vector", 
                                        max_DLTs = "vector", models = "vector", seed = "matrix"))

#' An S4 class to represent the simulation results of all methods for a specific scenario.
#'
#' @slot scenario The simulated scenario.
#' @slot PK_parameters The value for population parameters of the PK model from which the data of each data is simulated.
#' @slot omega_IIV The Inter-Individual Variation (IIV) value used as standard deviation in the log-normal distribution to compute individual PK parameters.
#' @slot omega_alpha The standard deviation for the measure of subject sensitivity (alpha is a sensitivity parameter of the patient assumed to come from a log-normal distribution) to compute subject-specific AUCs related to toxicity tolerance.
#' @slot cv The standard deviation of the normal distribution for the concentration proportional error model (residual error). 
#' @slot exposure_metric The reference exposure measurement to be used in the dose-exposure-toxicity relationship for simulating toxicity data (AUC or Cmax).
#' @slot prior_tox_real The prior toxicity for each dose computed using equation (13) of Ursino (2017) based on the parameters of the simulated scenario.
#' @slot N The total sample size per trial.
#' @slot n_trials The total number of trials to be simulated.
#' @slot n_sampling_timepoints  The number of concentration sampling time points obtained during each trial.
#' @slot real_time_sampling The real sampling time points.
#' @slot doses A vector with the doses panel.
#' @slot exposure_tox_threshold  The AUC threshold or Cmax threshold to be set before starting the trial to determine the toxicity outcome.
#' @slot n_chains The number of chains for the Stan model.
#' @slot n_iter The number of iterations for the Stan model.
#' @slot n_adapt  The number of warmup iterations for the Stan model.
#' @slot n_cores  The number of CPUs used to compute the H-MC chains for the Stan model.
#' @slot dose_allocation The percentage of MTD selection for each dose level after all trials starting from dose 0; equals to 0 if the trial has stopped before the end, according to the stopping rules.
#' @slot dose_selection A vector containing the next maximum tolerated doses (MTD) of each trial (TR); equals to 0 if the trial has stopped before the end, according to the stopping rules.
#' @slot targeted_tox  The toxicity target.
#' @slot median_DLTs The median number of DLTs observed in all trials.
#' @slot min_DLTs The minimum number of DLTs observed in all trials.
#' @slot max_DLTs The maximum number of DLTs observed in all trials.
#' @slot models A vector of character string to specify the selected dose-finding models to run the simulation.
#' @slot seed The seed of the random number generator that is used during the initialization phase and at the beginning of each trial.
#' @import methods
#' @export


#' An S4 class to get a better printing format for the above S4 class on simulation results of all methods for a specific scenario.
setGeneric("show")
#' @export 
setMethod(f = "show", signature = "simulation_res", definition = function(object)
{
  cat("Today: ", date(), "\n") 
  cat("\n","A. Scenario Settings (Scenario ", paste0(object@scenario,")"), "\n")
  cat("Number of simulations:", object@n_trials, "\n")
  cat("Total number of patients in each trial:", object@N, "\n")
  if (compartmental_model == "1-cmt"){
    cat("Population PK parameters (ka, Cl_pop, V_pop):", object@PK_parameters, "\n")
  } else {
    cat("Population PK parameters (ka, Cl_pop, V_1_pop, Q_pop, V_2_pop):", object@PK_parameters, "\n")
  }
  cat("Pourcentage of Inter-Individual Variability (IIV):", object@omega_IIV, "\n")
  cat("Standard deviation of subject-specific AUC with a sensitivity parameter to toxicity:", object@omega_alpha, "\n")
  cat("Standard deviation of the concentration proportional error model:", object@cv, "\n")
  cat("Reference exposure measurement used to simulate toxicity data:", object@exposure_metric, "\n")
  cat("Prior toxicity for all doses:", object@prior_tox_real, "\n")
  cat("\n","B. Data Summary (",object@models,"models)", "\n")
  cat("Cohort size for patients inclusion:", object@cohort_size, "\n")
  cat("Number of concentration sampling time points:", object@n_sampling_timepoints, "\n")
  cat("Sampling time after the drug administration (hours):", round(object@real_time_sampling, digits = 3), "\n")
  cat("Levels of Doses:", round(object@doses, digits=3), "\n")
  cat("Initialization seed for each method:",object@seed[1], "\n")
  cat("\n", "Rstan computation specs: \n")
  cat("The Stan models run with", object@n_chains, "MCMC chains ")
  cat("which each chain has", object@n_iter, "iterations ")
  cat("and", object@n_adapt, "warmup iterations ")
  cat("using", object@n_cores, "CPUs. \n")
  n <- object@N
  cat("\n","C. Dose-Finding Simulation Results: \n")
  models_name <- object@models
  matrix_selection <- matrix(NA, nrow=length(object@models)+2, ncol=length(object@doses)+1)
  matrix_allocation <- matrix(NA, nrow=length(object@models)+1, ncol=length(object@doses)+1)
  matrix_additional <- matrix(NA, nrow=length(object@models), ncol=3)
  rownames(matrix_selection) <- c("Dose", "Prior toxicity probability", models_name)
  rownames(matrix_allocation) <- c("Dose", models_name)
  rownames(matrix_additional) <- models_name
  colnames(matrix_allocation) <- rep("", length(object@doses)+1)
  colnames(matrix_selection) <- rep("", length(object@doses)+1)
  colnames(matrix_additional) <- c("Median (DLTs)", "Min (DLTs)", "Max DLTs")
  stop <- as.character("STOP")
  dd <- paste("", 1:length(object@doses), sep = "")
  matrix_selection[1, ] <- c(stop, dd)
  matrix_allocation[1, ] <- c(stop, dd)
  matrix_selection[2, ] <- c("NA", round(object@prior_tox_real, digits=3))
  for(i in 1:length(object@models)){
    matrix_selection[2+i, ] <- object@dose_selection[i, ]
    matrix_allocation[1+i, ] <- object@dose_allocation[i, ]
  }
  matrix_additional[, 1] <- object@median_DLTs
  matrix_additional[, 2] <- object@min_DLTs
  matrix_additional[, 3] <- object@max_DLTs
  cat("\n","I. Pourcentage (%) of dose selection \n")
  print(matrix_selection, quote = FALSE)
  cat("\n","II. Pourcentage (%) of dose allocation \n")
  print(matrix_allocation, quote = FALSE)
  cat("\n","III. Number of DLTs \n")
  print(matrix_additional, quote = FALSE)
  cat("Recommendation is based on a target toxicity probability of:",object@targeted_tox, "\n")
  cat("Seed of each trial for all methods:",object@seed[-1], "\n")
}
)