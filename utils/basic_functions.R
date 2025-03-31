########################################################
#### Basic (fundamental) functions for main scripts ####
########################################################

# Logit function
f_logit <- function(p){
  return(log(p/(1-p)))
}

# Inverse logit function
f_inv_logit <- function(x){
  return(1/(1+exp(-x)))
}

# Cloglog function
f_cloglog <- function(x){
  return(log(-log(1-x)))
}


# Stopping rule function
checking_stop_rule <- function(x, target, error){
  if(any(is.na(x)) == TRUE){
    x <- x[-which(is.na(x))]
  }
  sum(x>(target+error))/length(x)              
}

# Oral 1-cmt PK model with first order absorption - Equation for Michaelis Menten elimination kinetics
oral_1cmt_Michaelis_Menten_equation <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dC <- -(((Vmax/V)*C)/(Km+C)) + (dose/V)*ka*exp(-ka*t)
    dAUC <- C
    return(list(c(dC, dAUC)))
  })
}

# Prior ("real") probability of toxicity computed for two-compartment PK model with AUC as exposure metric function
f_prior_tox_AUC_2cmt <- function(PK_parameters, omega_IIV, omega_alpha, cv, doses, exposure_tox_threshold, n_sampling_timepoints, time_sampling){
  N <- 10000
  prior_tox_simu_data <- simu_data(
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    doses = doses,
    exposure_tox_threshold = exposure_tox_threshold,
    n_sampling_timepoints = n_sampling_timepoints,
    time_sampling = time_sampling,
    N = N,
    n_trials = 1,
    exposure_metric = "AUC",
    AUC_method = "integration",
    compartmental_model = "2-cmt",
    elimination_kinetics = "linear"
  )
  auc_per_dose <- prior_tox_simu_data[[1]]$simu_exact_exposure
  tox <- prior_tox_simu_data[[1]]$toxicity
  auc_with_DLT_per_dose <- lapply(1:length(doses), function (i) auc_per_dose[which(tox[,i] == 1),i])
  prior_tox_real_AUC <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose[,i]))
  return(prior_tox_real_AUC)
}

# Prior ("real") probability of toxicity computed for Cmax as exposure metric function
f_prior_tox_Cmax <- function(PK_parameters, omega_IIV, omega_alpha, cv, doses, exposure_tox_threshold, n_sampling_timepoints, time_sampling,
                             compartmental_model){
  N <- 10000
  prior_tox_simu_data <- simu_data(
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    doses = doses,
    exposure_tox_threshold = exposure_tox_threshold,
    n_sampling_timepoints = n_sampling_timepoints,
    time_sampling = time_sampling,
    N = N,
    n_trials = 1,
    exposure_metric = "Cmax",
    AUC_method = "integration",
    compartmental_model = compartmental_model,
    elimination_kinetics = "linear"
  )
  conc_per_dose <- lapply(1:length(doses), function(i) prior_tox_simu_data[[1]]$simu_exact_conc[[i]])
  cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N)
  for (i in 1:length(doses)){
    for (j in 1:N){
      cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
    }
  }
  tox <- prior_tox_simu_data[[1]]$toxicity
  cmax_with_DLT_per_dose <- lapply(1:length(doses), function (i) cmax_per_dose[which(tox[,i] == 1),i])
  prior_tox_real_Cmax <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
  
  return(prior_tox_real_Cmax)
}

# Prior ("real") probability of toxicity computed for AUC as exposure metric in the context of Michaelis-Menten elimination kinetics
f_prior_tox_AUC_MM <- function(PK_parameters, omega_IIV, omega_alpha, cv, doses, exposure_tox_threshold, n_sampling_timepoints,
                               time_sampling, compartmental_model){
  N <- 10000
  prior_tox_simu_data <- simu_data(
    PK_parameters = PK_parameters,
    omega_IIV = omega_IIV,
    omega_alpha = omega_alpha,
    cv = cv,
    doses = doses,
    exposure_tox_threshold = exposure_tox_threshold,
    n_sampling_timepoints = n_sampling_timepoints,
    time_sampling = time_sampling,
    N = N,
    n_trials = 1,
    exposure_metric = "AUC",
    AUC_method = "integration",
    compartmental_model = compartmental_model,
    elimination_kinetics = "Michaelis-Menten"
  )
  auc_per_dose_MM <- prior_tox_simu_data[[1]]$simu_exact_exposure
  tox <- prior_tox_simu_data[[1]]$toxicity
  auc_with_DLT_per_dose <- lapply(1:length(doses), function(i) auc_per_dose_MM[which(tox[,i] == 1),i])
  prior_tox_real_AUC_MM <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose_MM[,i]))
  
  return(prior_tox_real_AUC_MM)
}

# Defining several useful functions for the ED-EWOC method and other methods 

### ED-EWOC and ED AUC estimation function
AUC_estim_popPK <- function(N, doses, ka, Cl_ind, V_ind) {
  auc <- matrix(data = NA, nrow = N, ncol = length(doses))
  for (i in 1:length(doses)) {
    auc[,i] <- sapply(1:N, function(j) {
      integrand <- function(t) {
        (doses[i] / V_ind[j]) * (ka / (ka - (Cl_ind[j] / V_ind[j]))) * (exp(-(Cl_ind[j] / V_ind[j]) * t) - exp(-ka * t))
      }
      return(integrate(integrand, lower = 0, upper = Inf)$value)
    })
  }
  return(auc)
}

AUC_estim_popPK_stable <- function(N, doses, ka, Cl_ind, V_ind) {
  auc <- matrix(data = NA, nrow = N, ncol = length(doses))
  for (i in 1:length(doses)) {
    auc[,i] <- sapply(1:N, function(j) {
      integrand <- function(t) {
        (doses[i] / V_ind[j]) * (ka / (ka - (Cl_ind[j] / V_ind[j]))) * (exp(-(Cl_ind[j] / V_ind[j]) * t) - exp(-ka * t))
      }
      return(hcubature(integrand, lower = 0, upper = Inf)$integral)
    })
  }
  return(auc)
}

### Reference AUC computation function
ref_AUC_computation <- function(ref_dose, PK_parameters, cv,  time_sampling, compartmental_model, elimination_kinetics = "linear") {
  if (compartmental_model == "1-cmt") {
    if (elimination_kinetics == "linear") {
      ka <- PK_parameters[1]
      Cl_pop <- PK_parameters[2]
      V_pop <- PK_parameters[3]
      integrand <- function(t) {
        (ref_dose / V_pop) * (ka / (ka - (Cl_pop / V_pop))) * (exp(-(Cl_pop / V_pop) * t) - exp(-ka * t)) * (1 + rnorm(1, 0, cv))
      }
      ref_auc <- integrate(integrand, lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, stop.on.error = TRUE)$value
      return(ref_auc)
  
    } else if (elimination_kinetics == "Michaelis-Menten") {
      ka <- PK_parameters[1]
      V_pop <- PK_parameters[3]
      Vmax <- PK_parameters[4]
      Km <- PK_parameters[5]
      parameters <- c(dose = ref_dose, ka = ka, V = V_pop, Vmax = Vmax, Km = Km)
      state <- c(C = 0, AUC = 0)
      AUC_Inf_approximation <- 1000
      ref_auc <- ode(y = state, times = (0:AUC_Inf_approximation), func = oral_1cmt_Michaelis_Menten_equation, parms = parameters, method = "lsoda")[AUC_Inf_approximation, "AUC"]
      return(ref_auc)
      
    } else {
      stop("Error in elimination_kinetics: Invalid mechanism of drug elimination choosen for PK data. Please choose between 'linear' and 'Michaelis-Menten' elimination kinetics!")
    }
  
  } else if (compartmental_model == "2-cmt") {
    ka <- PK_parameters[1]
    Cl_pop <- PK_parameters[2]
    V_1_pop <- PK_parameters[3]
    Q_pop <- PK_parameters[4]
    V_2_pop <- PK_parameters[5]
    beta <- (1/2)*(Q_pop/V_1_pop + Q_pop/V_2_pop + Cl_pop/V_1_pop - sqrt((Q_pop/V_1_pop + Q_pop/V_2_pop + Cl_pop/V_1_pop)^2 - 4*(Q_pop/V_2_pop)*(Cl_pop/V_1_pop)))
    alpha <- (Q_pop/V_2_pop)*(Cl_pop/V_1_pop)/beta
    A <- (ka/V_1_pop)*(((Q_pop/V_2_pop)-alpha)/((ka-alpha)*(beta-alpha)))
    B <- (ka/V_1_pop)*(((Q_pop/V_2_pop)-beta)/((ka-beta)*(alpha-beta)))
    integrand <- function(t) {
      ref_dose*(A*exp(-alpha*t) + B*exp(-beta*t) - (A + B)*exp(-ka*t))
    }
    ref_auc <- integrate(integrand, lower = 0, upper = Inf)$value
    return(ref_auc)
    
  } else {
    stop("Error in PK compartmental model selection: Wrong value for parameter. Please type either '1-cmt' or '2-cmt' to choose between a one-compartment or a two-compartment model with first order absorption and linear elimination for oral admistration.")
  }
}

# Defining several useful functions for the PKLOGIT method 

### Logit function for PKLOGIT posterior probability of toxicity estimation
f_logit_estim <- function(v, ref_auc, lambda, parmt){
  f_inv_logit(-lambda[1]+exp(lambda[2])*(v-ref_auc))*dnorm(v,parmt[1],parmt[2])
}
### Logit function for PKLOGIT probabilities of toxicity estimation using parameters' values from MCMC chains
f2_logit_estim <- function(v, ref_auc, lambda1, lambda2, parmt1, parmt2){
  f_inv_logit(-lambda1+exp(lambda2)*(v-ref_auc))*dnorm(v,parmt1,parmt2)
}

### AUC estimation function for dfpk package (alternative solution since installation is failing)

pk.estim <-
  function(par,t,dose,conc){
    ka <- par[1]
    CL <- par[2]
    V <- par[3]
    s<-dose*ka/V/(ka-CL/V)* (exp(-CL/V*t)-exp(-ka*t))
    if (any(is.nan(s)) ) cat("ka=", ka, "\n V=",V,"\n CL=",CL, "\n")
    s[s==0] = rep(2^(-1074),length(s[s==0] ))
    sum( (log(s)-log(conc))^2 )
  }

AUC.estim <-
  function(t,conc,dose, method = 2){
    if(method == "1"){
      out_pk = optim(c(1,5,50), pk.estim, t=t, dose=dose, conc=conc, method = "L-BFGS-B", lower = c(0.1,0.2,1))
      dose/out_pk$par[2]
    }else if(method == "2"){
      a <- auc(conc, t, design="complete")
      a$est
    }else{
      stop("Error in AUCmethod: Invalid method for the AUC estimation. See the R documentation of the function for the possible methods")
    }
  }
