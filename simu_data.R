####################################################################################
### Simulate appropriate PK data depending on the number of trials and the model ###
####################################################################################

# Core function
simu_data <- function(PK_parameters, omega_IIV, omega_alpha, cv, doses, exposure_tox_threshold, n_sampling_timepoints, time_sampling, N, 
                      n_trials, exposure_metric = "AUC", AUC_method = "compartmental", compartmental_model = "1-cmt", 
                      elimination_kinetics = "linear", seed=12345) {

  simulated_conc_data <- list()
  
  for(tr in 1:n_trials){
    
    set.seed(seed + tr)
    # Computation of a list containing matrices on the concentrations of drug of each individual for each PK assessment time point for 
    # different potentials doses with and without a measure error used in the output
    # AND
    # Computation of a list containing toxicity vectors for each individual depending on the dose and the number of the trial
    # Assessment of toxicity outcome for simulated PK data based on the predicted exposure metric (AUC, TAUC or Cmax)
    conc <- list()
    conc_with_error <- list()
    auc <- matrix(data = NA, nrow = N, ncol = length(doses))
    tauc <- matrix(data = NA, nrow = N, ncol = length(doses))
    cmax <- matrix(data = NA, nrow = N, ncol = length(doses))
    s_auc <- matrix(data = NA, nrow = N, ncol = length(doses))
    s_tauc <- matrix(data = NA, nrow = N, ncol = length(doses))
    s_cmax <- matrix(data = NA, nrow = N, ncol = length(doses))
    tox <- matrix(data = NA, nrow = N, ncol = length(doses))
    
    if (compartmental_model == "1-cmt") {
      if (elimination_kinetics == "linear") {
        # Essential population and individual parameters of a one-compartment PK model with linear elimination and first order absorption
        ka <- PK_parameters[1]
        Cl_pop <- PK_parameters[2]
        V_pop <- PK_parameters[3]
        eta_Cl <- exp(rnorm(N, 0, omega_IIV))
        eta_V <- exp(rnorm(N, 0, omega_IIV))
        Cl_ind <- Cl_pop * eta_Cl
        V_ind <- V_pop * eta_V
        s_alpha <- exp(rnorm(N, 0, omega_alpha))
        # Computation of a matrix containing the value of individual PK parameters for each individual used to simulate PK data
        parameters_ind <- c(ka, Cl_ind, V_ind)
        # Computation of concentration values for a oral one-compartment PK model with linear elimination kinetics
        for (i in 1:length(doses)) {
          conc[[i]] <- matrix(data = NA, nrow = N, ncol = length(time_sampling))
          conc_with_error[[i]] <- matrix(data = NA, nrow = N, ncol = length(time_sampling))
          for (j in 1:N) {
            for (k in 1:length(time_sampling)) {
              t <- time_sampling[k]
              # Computation of exact concentrations
              conc[[i]][j, k] <- (doses[i] / V_ind[j]) * (ka / (ka - (Cl_ind[j] / V_ind[j]))) * (exp(-(Cl_ind[j] / V_ind[j]) * t) - exp(-ka * t))
              # Computation of concentrations with proportional error model
              conc_with_error[[i]][j, k] <- conc[[i]][j, k] * (1 + rnorm(1, 0, cv))
            }
            # Computation of AUC when selected as the exposure metric for the scenario
            if (exposure_metric == "AUC") {
              if (AUC_method == "integration") {
                # Computation of AUC using integration method
                integrand <- function(t) {
                  (doses[i] / V_ind[j]) * (ka / (ka - (Cl_ind[j] / V_ind[j]))) * (exp(-(Cl_ind[j] / V_ind[j]) * t) - exp(-ka * t))
                }
                auc[j, i] <- integrate(integrand, lower = 0, upper = Inf)$value
              } else if (AUC_method == "compartmental") {
                # Computation of AUC using proportional formula between AUC, doses and clearance (CL)
                auc[j, i] <- doses[i] / Cl_ind[j]
              } else {
                stop("Error in AUC_method: Invalid method for the AUC estimation. Please choose between 'integration' (computation based on integral method) or 'compartmental' (computation based on assumed porportionality between AUC, dose and CL).")
              }
              # Computation of a subject-specific toxicity threshold based on AUC using a linear function of AUC (s(AUC))
              s_auc[j, i] <- s_alpha[j]*auc[j, i]
              # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
              if (s_auc[j, i] >= exposure_tox_threshold) {
                tox[j, i] <- 1
              } else {
                tox[j, i] <- 0
              }
            } else if (exposure_metric == "TAUC") {
              # Computation of TAUC using integration method when selected as the exposure metric for the scenario
              integrand <- function(t) {
                (doses[i] / V_ind[j]) * (ka / (ka - (Cl_ind[j] / V_ind[j]))) * (exp(-(Cl_ind[j] / V_ind[j]) * t) - exp(-ka * t))
              }
              tauc[j, i] <- integrate(integrand, lower = 0, upper = tail(time_sampling, n = 1))$value
              # Computation of a subject-specific toxicity threshold based on truncated AUC using a linear function of TAUC (s(TAUC))
              s_tauc[j, i] <- s_alpha[j]*tauc[j, i]
              # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
              if (s_tauc[j, i] >= exposure_tox_threshold) {
                tox[j, i] <- 1
              } else {
                tox[j, i] <- 0
              } 
            } else if (exposure_metric == "Cmax") {
              # Extraction of the maximal concentration (Cmax) for each subject at each potential dose
              cmax[j, i] <- max(conc[[i]][j,])
              # Computation of a subject-specific toxicity threshold based on Cmax using a linear function of Cmax (s(Cmax))
              s_cmax[j, i] <- s_alpha[j]*cmax[j, i]
              # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
              if (s_cmax[j, i] >= exposure_tox_threshold) {
                tox[j, i] <- 1
              } else {
                tox[j, i] <- 0
              }
            } else {
              stop("Error in exposure_metric: Invalid metric of exposure choosen for PK data. Please choose between 'AUC', 'TAUC' and 'Cmax' to use one of the two available metrics as measure of exposure.")
            }
          }
        }
      } else if (elimination_kinetics == "Michaelis-Menten") {
        # Computation of concentration values for a oral one-compartment PK model with Michaelis-Menten elimination kinetics
        ka <- PK_parameters[1]
        V_pop <- PK_parameters[3]
        eta_V <- exp(rnorm(N, 0, omega_IIV))
        V_ind <- V_pop * eta_V
        Vmax <- PK_parameters[4] # Maximum rate of metabolism (Michaelis-Menten)
        Km <- PK_parameters[5] # Michaelis constant
        s_alpha <- exp(rnorm(N, 0, omega_alpha))
        # Computation of a matrix containing the value of individual PK parameters for each individual used to simulate PK data
        parameters_ind <- c(ka, V_ind, Vmax, Km)
        state <- c(C = 0, AUC = 0)
        for (i in 1:length(doses)) {
          conc[[i]] <- matrix(data = NA, nrow = N, ncol = length(time_sampling))
          conc_with_error[[i]] <- matrix(data = NA, nrow = N, ncol = length(time_sampling))
          for (j in 1:N) {
            parameters <- c(dose = doses[i], ka = ka, Vmax = Vmax, V = V_ind[j], Km = Km)
            # Computation of exact concentrations
            conc[[i]][j,] <- ode(y = state, times = time_sampling, func = oral_1cmt_Michaelis_Menten_equation, parms = parameters, method = "lsoda")[, "C"]
            # Computation of concentrations with proportional error model
            conc_with_error[[i]][j,] <- conc[[i]][j,] * (1 + rnorm(1, 0, cv))
            # Computation of AUC when selected as the exposure metric for the scenario
            if (exposure_metric == "AUC") {
              # Computation of AUC using integration method
              AUC_Inf_approximation <- 900
              auc[j, i] <- ode(y = state, times = (0:AUC_Inf_approximation), func = oral_1cmt_Michaelis_Menten_equation, parms = parameters, method = "lsoda")[AUC_Inf_approximation, "AUC"]
              # Computation of a subject-specific toxicity threshold based on AUC using a linear function of AUC (s(AUC))
              s_auc[j, i] <- s_alpha[j]*auc[j, i]
              # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
              if (s_auc[j, i] >= exposure_tox_threshold) {
                tox[j, i] <- 1
              } else {
                tox[j, i] <- 0
                }
              } else if (exposure_metric == "TAUC") {
              # Computation of TAUC using integration method when selected as the exposure metric for the scenario
              tauc[j, i] <- ode(y = state, times = time_sampling, func = oral_1cmt_Michaelis_Menten_equation, parms = parameters, method = "lsoda")[tail(time_sampling, n = 1), "AUC"]
              # Computation of a subject-specific toxicity threshold based on truncated AUC using a linear function of TAUC (s(TAUC))
              s_tauc[j, i] <- s_alpha[j]*tauc[j, i]
              # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
              if (s_tauc[j, i] >= exposure_tox_threshold) {
                tox[j, i] <- 1
              } else {
                tox[j, i] <- 0
              } 
            } else if (exposure_metric == "Cmax") {
              # Extraction of the maximal concentration (Cmax) for each subject at each potential dose
              cmax[j, i] <- max(conc[[i]][j,])
              # Computation of a subject-specific toxicity threshold based on Cmax using a linear function of Cmax (s(Cmax))
              s_cmax[j, i] <- s_alpha[j]*cmax[j, i]
              # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
              if (s_cmax[j, i] >= exposure_tox_threshold) {
                tox[j, i] <- 1
              } else {
                tox[j, i] <- 0
              }
            } else {
              stop("Error in exposure_metric: Invalid metric of exposure choosen for PK data. Please choose between 'AUC', 'TAUC' and 'Cmax' to use one of the two available metrics as measure of exposure.")
            }
          }
        }
      } else {
        stop("Error in elimination_kinetics: Invalid mechanism of drug elimination choosen for PK data. Please choose between 'linear' and 'Michaelis-Menten' elimination kinetics!")
      }
    } else if (compartmental_model == "2-cmt") {
      # Essential population and individual parameters of a two-compartment PK model with linear elimination and first order absorption
      ka <- PK_parameters[1]
      Cl_pop <- PK_parameters[2]
      V_1_pop <- PK_parameters[3]
      Q_pop <- PK_parameters[4]
      V_2_pop <- PK_parameters[5]
      eta_Cl <- exp(rnorm(N, 0, omega_IIV))
      eta_V_1 <- exp(rnorm(N, 0, omega_IIV))
      eta_Q <- exp(rnorm(N, 0, omega_IIV))
      eta_V_2 <- exp(rnorm(N, 0, omega_IIV))
      Cl_ind <- Cl_pop * eta_Cl
      V_1_ind <- V_1_pop * eta_V_1
      Q_ind <- Q_pop * eta_Q
      V_2_ind <- V_2_pop * eta_V_2
      s_alpha <- exp(rnorm(N, 0, omega_alpha))
      # Computation of a matrix containing the value of individual PK parameters for each trial used in the output
      parameters_ind <- c(ka, Cl_ind, V_1_ind, Q_ind, V_2_ind)
      for (i in 1:length(doses)) {
        conc[[i]] <- matrix(data = NA, nrow = N, ncol = length(time_sampling))
        conc_with_error[[i]] <- matrix(data = NA, nrow = N, ncol = length(time_sampling))
        for (j in 1:N) {
          # Creation of storage variables to make the final concentration formula easier to read and write
          beta <- (1/2)*(Q_ind[j]/V_1_ind[j] + Q_ind[j]/V_2_ind[j] + Cl_ind[j]/V_1_ind[j] - sqrt((Q_ind[j]/V_1_ind[j] + Q_ind[j]/V_2_ind[j] + Cl_ind[j]/V_1_ind[j])^2 - 4*(Q_ind[j]/V_2_ind[j])*(Cl_ind[j]/V_1_ind[j])))
          alpha <- (Q_ind[j]/V_2_ind[j])*(Cl_ind[j]/V_1_ind[j])/beta
          A <- (ka/V_1_ind[j])*(((Q_ind[j]/V_2_ind[j])-alpha)/((ka-alpha)*(beta-alpha)))
          B <- (ka/V_1_ind[j])*(((Q_ind[j]/V_2_ind[j])-beta)/((ka-beta)*(alpha-beta)))
          for (k in 1:length(time_sampling)) {
            t <- time_sampling[k]
            # Computation of exact concentrations
            conc[[i]][j, k] <- doses[i]*(A*exp(-alpha*t) + B*exp(-beta*t) - (A + B)*exp(-ka*t))
            # Computation of concentrations with proportional error model
            conc_with_error[[i]][j, k] <- conc[[i]][j, k] * (1 + rnorm(1, 0, cv))
          }
          if (exposure_metric == "AUC") {
            # Computation of AUC using integration method when selected as the exposure metric for the scenario
            integrand <- function(t) {
              doses[i]*(A*exp(-alpha*t) + B*exp(-beta*t) - (A + B)*exp(-ka*t))
            }
            auc[j, i] <- integrate(integrand, lower = 0, upper = Inf)$value
            # Computation of a subject-specific toxicity threshold based on AUC or Cmax using a linear function of AUC (s(AUC))
            s_auc[j, i] <- s_alpha[j]*auc[j, i]
            # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
            if (s_auc[j, i] >= exposure_tox_threshold) {
              tox[j, i] <- 1
            } else{
              tox[j, i] <- 0
            }
          } else if (exposure_metric == "TAUC") {
            # Computation of TAUC using integration method when selected as the exposure metric for the scenario
            integrand <- function(t) {
              doses[i]*(A*exp(-alpha*t) + B*exp(-beta*t) - (A + B)*exp(-ka*t))
            }
            tauc[j, i] <- integrate(integrand, lower = 0, upper = tail(time_sampling, n = 1))$value
            # Computation of a subject-specific toxicity threshold based on truncated AUC using a linear function of TAUC (s(TAUC))
            s_tauc[j, i] <- s_alpha[j]*tauc[j, i]
            # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
            if (s_tauc[j, i] >= exposure_tox_threshold) {
              tox[j, i] <- 1
            } else {
              tox[j, i] <- 0
            } 
        } else if (exposure_metric == "Cmax") {
            # Extraction of the maximal concentration (Cmax) for each subject at each potential dose
            cmax[j, i] <- max(conc[[i]][j,])
            # Computation of a subject-specific toxicity threshold based on Cmax using a linear function of Cmax (s(Cmax))
            s_cmax[j, i] <- s_alpha[j]*cmax[j, i]
            # Computation of toxicity matrix to reference DLT event based on subject-specific toxicity threshold using simulated PK data
            if (s_cmax[j, i] >= exposure_tox_threshold) {
              tox[j, i] <- 1
            } else {
              tox[j, i] <- 0
            }
          } else {
            stop("Error in exposure_metric: Invalid metric of exposure choosen for PK data. Please choose between 'AUC', 'TAUC' and 'Cmax' to use one of the two available metrics as measure of exposure.")
          }
        }
      }
    } else {
      stop("Error in PK compartmental model selection: Wrong value for parameter. Please type either '1-cmt' or '2-cmt' to choose between a one-compartment or a two-compartment model with first order absorption and linear elimination for oral admistration.")
    }
    if (exposure_metric == "AUC") {
      # Exposure metric: AUC
      simulated_conc_data[[tr]] <- list("simu_exact_conc" = conc, "simu_conc_with_error" = conc_with_error, "simu_exact_exposure" = auc,
                                        "sensitivity_auc" = s_auc, "toxicity" = tox, "individual_parameters" = parameters_ind, "alpha" = s_alpha)
    } else if (exposure_metric == "TAUC") {
      # Exposure metric: Truncated AUC (0h-48h)
      simulated_conc_data[[tr]] <- list("simu_exact_conc" = conc, "simu_conc_with_error" = conc_with_error, "simu_exact_exposure" = tauc,
                                        "sensitivity_exposure" = s_tauc, "toxicity" = tox, "individual_parameters" = parameters_ind, "alpha" = s_alpha)
    } else if (exposure_metric == "Cmax") {
      # Exposure metric: Cmax
      simulated_conc_data[[tr]] <- list("simu_exact_conc" = conc, "simu_conc_with_error" = conc_with_error, "simu_exact_exposure" = cmax,
                                        "sensitivity_exposure" = s_cmax, "toxicity" = tox, "individual_parameters" = parameters_ind, "alpha" = s_alpha)
    }
  }
  return(simulated_conc_data)
}