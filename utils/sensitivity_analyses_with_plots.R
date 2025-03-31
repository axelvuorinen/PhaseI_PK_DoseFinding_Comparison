# SENSITIVITY ANALYSIS FOR CHOICE OF 2-CMT PK MODEL AUC THRESHOLD
N_simu <- 10000
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
time_sampling <- seq(0, 24, length.out = 48)
cv <- 0.20
threshold_AUC_dose_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 5, # Random threshold, doesn't matter for AUC histograms
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
auc_per_dose <- threshold_AUC_dose_2cmt_PKdata[[1]]$simu_exact_exposure

hist(auc_per_dose[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "AUC", main = "Histogram of AUC of all simulated patients for dose 1")
lines(density(auc_per_dose[,1]), lwd = 2, col = "chocolate")
hist(auc_per_dose[,2], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "AUC", main = "Histogram of AUC of all simulated patients for dose 2")
lines(density(auc_per_dose[,2]), lwd = 2, col = "chocolate")
hist(auc_per_dose[,3], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "AUC", main = "Histogram of AUC of all simulated patients for dose 3")
lines(density(auc_per_dose[,3]), lwd = 2, col = "chocolate")
hist(auc_per_dose[,4], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "AUC", main = "Histogram of AUC of all simulated patients for dose 4")
lines(density(auc_per_dose[,4]), lwd = 2, col = "chocolate")

# Threshold at first glance (dose-dependent) : 5, 10, 20, 40

# MTD = Dose 1 / Threshold = 6.2 -> p_tox = (0.2508, 0.4253, 0.6563, 0.7996)

threshold_AUC_DLT_dose1_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 6.2,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
auc_per_dose <- threshold_AUC_DLT_dose1_2cmt_PKdata[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose1_2cmt_PKdata[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function (i) auc_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose1 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose[,i]))
prior_tox_real_dose1


# MTD = Dose 2 / Threshold = 10.2 -> p_tox = (0.1252, 0.2535, 0.4701, 0.6476) 

threshold_AUC_DLT_dose2_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 10.2,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
auc_per_dose <- threshold_AUC_DLT_dose2_2cmt_PKdata[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose2_2cmt_PKdata[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function (i) auc_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose2 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose[,i]))
prior_tox_real_dose2


# MTD = Dose 3 / Threshold = 18.9 -> p_tox = (0.0414, 0.1040, 0.2521, 0.4142)

threshold_AUC_DLT_dose3_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 18.9,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
auc_per_dose <- threshold_AUC_DLT_dose3_2cmt_PKdata[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose3_2cmt_PKdata[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function (i) auc_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose3 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose[,i]))
prior_tox_real_dose3


# MTD = Dose 4 / Threshold = 30.3 -> p_tox = (0.0143, 0.0435, 0.1307, 0.2526)

threshold_AUC_DLT_dose4_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 30.3,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
auc_per_dose <- threshold_AUC_DLT_dose4_2cmt_PKdata[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose4_2cmt_PKdata[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function (i) auc_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose4 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose[,i]))
prior_tox_real_dose4


# MTD = Dose 0 / Threshold = 3.8 -> p_tox = (0.4195 0.6096, 0.8031, 0.8996)

threshold_AUC_DLT_dose0_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 3.8,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
auc_per_dose <- threshold_AUC_DLT_dose0_2cmt_PKdata[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose0_2cmt_PKdata[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function (i) auc_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose0 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose[,i]))
prior_tox_real_dose0




# SENSITIVITY ANALYSIS FOR CHOICE OF CMAX THRESHOLD FOR 1-CMT PK MODEL WITH LINEAR ELIMINATION KINETICS (HISTOGRAM)

N_simu <- 10000
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
time_sampling <- seq(0, 24, length.out = 48)
cv <- 0.20
threshold_Cmax_dose_1cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 5, # Random threshold, doesn't matter for Cmax histograms
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "1-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_dose_1cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}

hist(cmax_per_dose[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 1")
lines(density(cmax_per_dose[,1]), lwd = 2, col = "chocolate")
hist(cmax_per_dose[,2], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 2")
lines(density(cmax_per_dose[,2]), lwd = 2, col = "chocolate")
hist(cmax_per_dose[,3], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 3")
lines(density(cmax_per_dose[,3]), lwd = 2, col = "chocolate")
hist(cmax_per_dose[,4], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 4")
lines(density(cmax_per_dose[,4]), lwd = 2, col = "chocolate")

# Threshold at first glance (dose-dependent) : 0.5, 1.5, 2.5, 4

# MTD = Dose 1 / Threshold = 0.49 -> p_tox = (0.2558, 0.4348, 0.6690, 0.8179)

threshold_Cmax_DLT_dose1_1cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 0.49,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "1-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose1_1cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose1_1cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose1 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose1


# MTD = Dose 2 / Threshold = 0.82 -> p_tox = (0.12, 0.2520, 0.4720, 0.6529) 

threshold_Cmax_DLT_dose2_1cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 0.82,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "1-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose2_1cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose2_1cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose2 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose2


# MTD = Dose 3 / Threshold = 1.51 -> p_tox = (0.0391, 0.1035, 0.2529, 0.4200)

threshold_Cmax_DLT_dose3_1cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 1.51,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "1-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose3_1cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose3_1cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose3 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose3


# MTD = Dose 4 / Threshold = 2.43 -> p_tox = (0.0139, 0.0422, 0.1292, 0.2524)

threshold_Cmax_DLT_dose4_1cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 2.43,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "1-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose4_1cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose4_1cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose4 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose4


# MTD = Dose 0 / Threshold = 0.31 -> p_tox = (0.4174, 0.6102, 0.8143, 0.9138)

threshold_Cmax_DLT_dose0_1cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 0.31,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "1-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose0_1cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose0_1cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose0 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose0



# SENSITIVITY ANALYSIS FOR CHOICE OF CMAX THRESHOLD FOR 2-CMT PK MODEL WITH LINEAR ELIMINATION KINETICS (HISTOGRAM)

N_simu <- 10000
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
time_sampling <- seq(0, 24, length.out = 48)
cv <- 0.20
threshold_Cmax_dose_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 5, # Random threshold, doesn't matter for Cmax histograms
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_dose_2cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}

hist(cmax_per_dose[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 1")
lines(density(cmax_per_dose[,1]), lwd = 2, col = "chocolate")
hist(cmax_per_dose[,2], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 2")
lines(density(cmax_per_dose[,2]), lwd = 2, col = "chocolate")
hist(cmax_per_dose[,3], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 3")
lines(density(cmax_per_dose[,3]), lwd = 2, col = "chocolate")
hist(cmax_per_dose[,4], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "Cmax", main = "Histogram of Cmax of all simulated patients for dose 4")
lines(density(cmax_per_dose[,4]), lwd = 2, col = "chocolate")

# Threshold at first glance (dose-dependent) : 0.5, 1.5, 2.5, 4

# MTD = Dose 1 / Threshold = 0.45 -> p_tox = (0.2530, 0.4391, 0.6865, 0.8312)

threshold_Cmax_DLT_dose1_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 0.45,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose1_2cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose1_2cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose1 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose1


# MTD = Dose 2 / Threshold = 0.75 -> p_tox = (0.1157, 0.2511, 0.4794, 0.6737) 

threshold_Cmax_DLT_dose2_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 0.75,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose2_2cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose2_2cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose2 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose2


# MTD = Dose 3 / Threshold = 1.38 -> p_tox = (0.0329, 0.0977, 0.2525, 0.4259)

threshold_Cmax_DLT_dose3_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 1.38,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose3_2cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose3_2cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose3 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose3


# MTD = Dose 4 / Threshold = 2.22 -> p_tox = (0.0110, 0.0358, 0.1227, 0.2516)

threshold_Cmax_DLT_dose4_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 2.22,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose4_2cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose4_2cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose4 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose4


# MTD = Dose 0 / Threshold = 0.29 -> p_tox = (0.4139, 0.6199, 0.8226, 0.9182)

threshold_Cmax_DLT_dose0_2cmt_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 0.29,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "Cmax",
  AUC_method = "integration",
  compartmental_model = "2-cmt",
  elimination_kinetics = "linear"
)
conc_per_dose <- lapply(1:length(doses), function(i) threshold_Cmax_DLT_dose0_2cmt_PKdata[[1]]$simu_exact_conc[[i]])
cmax_per_dose <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)){
  for (j in 1:N_simu){
    cmax_per_dose[j,i] <- max(conc_per_dose[[i]][j,])
  }
}
tox <- threshold_Cmax_DLT_dose0_2cmt_PKdata[[1]]$toxicity
cmax_with_DLT_per_dose <- lapply(1:length(doses), function(i) cmax_per_dose[which(tox[,i] == 1),i])
prior_tox_real_dose0 <- sapply(1:length(doses), function(i) length(cmax_with_DLT_per_dose[[i]])/length(cmax_per_dose[,i]))
prior_tox_real_dose0



#####
#####
#####


# ADAPTATION OF THE EWOC DECISION CRITERION FOR DOSE ALLOCATION: MICALLEF ED-EWOC METHOD
N <- 10000
n_trials <- 1
omega_IIV <- 0.2
test_ED_EWOC_PKdata <- simu_data(
  PK_parameters = test_PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N,
  n_trials = n_trials
)

model_BLRM <- stan_model("stan/blrm.stan", "BLRM model")
model_ED_EWOC <- stan_model("stan/ed_ewoc_logit.stan", "ED-EWOC logit model")

set.seed(12345)
conc <- test_ED_EWOC_PKdata[[1]]$simu_exact_conc
conc_with_error <- test_ED_EWOC_PKdata[[1]]$simu_conc_with_error
exact_auc <- test_ED_EWOC_PKdata[[1]]$simu_exact_exposure
sensitivity_auc <- test_ED_EWOC_PKdata[[1]]$sensitivity_exposure
tox <- test_ED_EWOC_PKdata[[1]]$toxicity
Cl_ind <- test_ED_EWOC_PKdata[[1]]$individual_parameters[2:(N+1)]
V_ind <- test_ED_EWOC_PKdata[[1]]$individual_parameters[(N+2):((2*N)+1)]
alpha <- test_ED_EWOC_PKdata[[1]]$alpha

x <- c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4))
y <- tox[cbind(1:length(x), x)]
conci <- c()
M <- N / cohort_size
j <- 1:cohort_size
auc_s <- c()
for(k in 1:N) {
  temp_conc <- as.vector(conc_with_error[[x[k]]][k, ])
  conci <- c(conci, temp_conc[real_sampling])
}

### The trial has been conducted until the end without any stopping
# For Stan
final_model <- ed_ewoc_model(num = N, doses = doses, ref_dose = ref_dose, x = x, y = y, conci = conci, 
                             n_sampling_timepoints = n_sampling_timepoints, real_sampling = real_sampling, PK_parameters = PK_parameters, Cl_ind = Cl_ind, V_ind = V_ind,
                             exposure_tox_threshold = exposure_tox_threshold, targeted_tox = targeted_tox, options = options, seed = seed, 
                             proba = 0.9)
MtD <- final_model$new_dose

final_model$p_estim_dose_mean
targeted_tox_lower_bound <- 0.15
overdosing_threshold_prob <- 0.35
overdosing_risk_prob <- 0.25
candidate_doses <- c()
proba_overdose_ctrl <- rep(NA, length(doses))
proba_targeted_tox <- rep(NA, length(doses))
for (l in 1:length(doses)) {
  if (all(is.na(p_estim_sum[,l])) == FALSE) {
    proba_overdose_ctrl[l] <- sum(p_estim_sum[which(!is.na(p_estim_sum[,l])),l] > overdosing_threshold) / length(p_estim_sum[which(!is.na(p_estim_sum[,l])),l] )
    if (proba_overdose_ctrl[l] < overdosing_risk) {
      candidate_doses <- c(candidate_doses,l)
      proba_targeted_tox[l] <- sum(p_estim_sum[which(!is.na(p_estim_sum[,l])),l] > targeted_tox_inf & p_estim_sum[which(!is.na(p_estim_sum[,l])),l] <= overdosing_threshold) / length(p_estim_sum[which(!is.na(p_estim_sum[,l])),l])
    }
  }
}
print(proba_overdose_ctrl)
print(proba_targeted_tox)



#####
#####
#####

# COMPARING THE DISTRIBUTION OF THE SIMULATED CONCENTRATIONS USING A 2-CMT PK MODEL TO THE ONE OF 
# THE SIMULATED CONCENTRATIONS USING 1-CMT PK MODEL FOR DOSE 1, 2, 3 AND 4

### Settings (1-CMT PK MODEL)
N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
targeted_tox <- 0.25
time_sampling <- seq(0, 24, by = 0.1)
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
exposure_metric <- "AUC"
AUC_method <- "compartmental"
exposure_tox_threshold <- 6.6
# PK data simulation for each dose
simulated_PK_data_1cmt <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model
)

### Settings - Old parameters' choice (2-CMT PK MODEL)
N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
targeted_tox <- 0.25
time_sampling <- seq(0, 24, by = 0.1)
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 1
V_2_pop <- 10
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
exposure_metric <- "AUC" 
AUC_method <- "integration"
exposure_tox_threshold <- 4.4
# PK data simulation for each dose
simulated_PK_data_2cmt_naive <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model
)

### Settings - New parameters' choice (2-CMT PK MODEL)
N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
targeted_tox <- 0.25
time_sampling <- seq(0, 24, by = 0.1)
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
exposure_metric <- "AUC" 
AUC_method <- "integration"
exposure_tox_threshold <- 4.4
# PK data simulation for each dose
simulated_PK_data_2cmt_informed <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model
)

### CONCENTRATION PLOTS TO COMPARE THE CURVE SHAPE/DIFFERENCE BETWEEN THE TWO WAY OF SIMULATING THE DATA USING COMPARTMENTAL MODEL WITH OLD CHOICE OF PARAMETERS
pop_mean_concentrations_1cmt_dose1 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[1]][,t]))
pop_mean_concentrations_2cmt_dose1 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_2cmt_naive[[1]]$simu_exact_conc[[1]][,t]))
#plot(x = time_sampling, y = pop_mean_concentrations_1cmt_dose1, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curves for dose 1 obtained by linking the samplings comparing 1-cmt and old 2-cmt PK models", col = "red", lwd=2)
#lines(x = time_sampling, y = pop_mean_concentrations_2cmt_dose1, type = "l", col = "chocolate", lwd=2)
#legend(18, 0.3, legend = c("One-compartment curve", "Two-compartment curve"), col = c("red", "chocolate"), lty=1, lwd=3, cex=0.8)
pop_mean_concentrations_dose1_df <- data.frame(time_sampling = time_sampling, pk_model = c(rep("1cmt", length(time_sampling)), rep("2cmt", length(time_sampling))), pop_mean_concentrations_dose1 = c(pop_mean_concentrations_1cmt_dose1, pop_mean_concentrations_2cmt_dose1))
pop_mean_concentrations_dose1_df$pk_model <- factor(pop_mean_concentrations_dose1_df$pk_model, levels = c("1cmt", "2cmt"), labels = c("One-compartment", "Two-compartment"))
ggplot(data = pop_mean_concentrations_dose1_df, aes(x = time_sampling, y = pop_mean_concentrations_dose1, color = pk_model)) +
  geom_line(lwd = 1, lty = 1) +
  labs(y = "Concentration (mg/L)", x = "Time (hours)", title = "Population concentration curves for dose 1 obtained by linking the samplings comparing 1-cmt and old 2-cmt PK models", color = "Population concentraton curve for the PK model") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.position = "bottom")

pop_mean_concentrations_1cmt_dose4 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[4]][,t]))
pop_mean_concentrations_2cmt_dose4 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_2cmt_naive[[1]]$simu_exact_conc[[4]][,t]))
#plot(x = time_sampling, y = pop_mean_concentrations_1cmt_dose4, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curves for dose 4 obtained by linking the samplings comparing 1-cmt and old 2-cmt PK models", col = "red", lwd=2)
#lines(x = time_sampling, y = pop_mean_concentrations_2cmt_dose4, type = "l", col = "chocolate", lwd=2)
#legend(18, 1.4, legend = c("One-compartment curve", "Two-compartment curve"), col = c("red", "chocolate"), lty=1, lwd=3, cex=0.8)
pop_mean_concentrations_dose4_df <- data.frame(time_sampling = time_sampling, pk_model = c(rep("1cmt", length(time_sampling)), rep("2cmt", length(time_sampling))), pop_mean_concentrations_dose4 = c(pop_mean_concentrations_1cmt_dose4, pop_mean_concentrations_2cmt_dose4))
pop_mean_concentrations_dose4_df$pk_model <- factor(pop_mean_concentrations_dose4_df$pk_model, levels = c("1cmt", "2cmt"), labels = c("One-compartment", "Two-compartment"))
ggplot(data = pop_mean_concentrations_dose4_df, aes(x = time_sampling, y = pop_mean_concentrations_dose4, color = pk_model)) +
  geom_line(lwd = 1, lty = 1) +
  labs(y = "Concentration (mg/L)", x = "Time (hours)", title = "Population concentration curves for dose 4 obtained by linking the samplings comparing 1-cmt and new 2-cmt PK models", color = "Population concentraton curve for the PK model") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.position = "bottom")

### CONCENTRATION PLOTS TO COMPARE THE CURVE SHAPE/DIFFERENCE BETWEEN THE TWO WAY OF SIMULATING THE DATA USING COMPARTMENTAL MODEL WITH NEW CHOICE OF PARAMETERS
pop_mean_concentrations_1cmt_dose1 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[1]][,t]))
pop_mean_concentrations_2cmt_dose1 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[1]][,t]))
#plot(x = time_sampling, y = pop_mean_concentrations_1cmt_dose1, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curves for dose 1 obtained by linking the samplings comparing 1-cmt and new 2-cmt simulated concentrations", col = "red", lwd=2)
#lines(x = time_sampling, y = pop_mean_concentrations_2cmt_dose1, type = "l", col = "chocolate", lwd=2)
#legend(18, 0.3, legend = c("One-compartment curve", "Two-compartment curve"), col = c("red", "chocolate"), lty=1, lwd=3, cex=0.8)
pop_mean_concentrations_dose1_df <- data.frame(time_sampling = time_sampling, pk_model = c(rep("1cmt", length(time_sampling)), rep("2cmt", length(time_sampling))), pop_mean_concentrations_dose1 = c(pop_mean_concentrations_1cmt_dose1, pop_mean_concentrations_2cmt_dose1))
pop_mean_concentrations_dose1_df$pk_model <- factor(pop_mean_concentrations_dose1_df$pk_model, levels = c("1cmt", "2cmt"), labels = c("One-compartment", "Two-compartment"))
ggplot(data = pop_mean_concentrations_dose1_df, aes(x = time_sampling, y = pop_mean_concentrations_dose1, color = pk_model)) +
  geom_line(lwd = 1, lty = 1) +
  labs(y = "Concentration (mg/L)", x = "Time (hours)", title = "Population concentration curves for dose 1 obtained by linking the samplings comparing 1-cmt and new 2-cmt PK models", color = "Population concentraton curve for the PK model") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.position = "bottom")

pop_mean_concentrations_1cmt_dose4 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[4]][,t]))
pop_mean_concentrations_2cmt_dose4 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[4]][,t]))
#plot(x = time_sampling, y = pop_mean_concentrations_1cmt_dose4, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curves for dose 4 obtained by linking the samplings comparing comparing 1-cmt and new 2-cmt PK models", col = "red", lwd=2)
#lines(x = time_sampling, y = pop_mean_concentrations_2cmt_dose4, type = "l", col = "chocolate", lwd=2)
#legend(18, 1.4, legend = c("One-compartment curve", "Two-compartment curve"), col = c("red", "chocolate"), lty=1, lwd=3, cex=0.8)
pop_mean_concentrations_dose4_df <- data.frame(time_sampling = time_sampling, pk_model = c(rep("1cmt", length(time_sampling)), rep("2cmt", length(time_sampling))), pop_mean_concentrations_dose4 = c(pop_mean_concentrations_1cmt_dose4, pop_mean_concentrations_2cmt_dose4))
pop_mean_concentrations_dose4_df$pk_model <- factor(pop_mean_concentrations_dose4_df$pk_model, levels = c("1cmt", "2cmt"), labels = c("One-compartment", "Two-compartment"))
pop_mean_concentrations_dose4_plot <- ggplot(data = pop_mean_concentrations_dose4_df, aes(x = time_sampling, y = pop_mean_concentrations_dose4, color = pk_model)) +
  geom_line(lwd = 1, lty = 1) +
  labs(y = "Concentration (mg/L)", x = "Time (hours)", title = "", color = "Population concentraton curve for the PK model") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold", size = 12), legend.position = "bottom", axis.title = element_text(size=12), legend.text = element_text(size=12))
#ggsave(paste0("plots/1cmt_vs_informed2cmt_concentration_curve_dose4.png"), plot = wrap_elements(pop_mean_concentrations_dose4_plot), dpi = "retina", units = "in", width = 14, height = 8, limitsize = FALSE)

                                                                                                                                                                                                                               
### CONCENTRATIONS CURVES WITH MEDIAN, 5-PERCENTILE AND 95-PERCENTILE FOR DOSE 1 AND DOSE 4 USING SIMULATED CONCENTRATIONS WITH 1-CMT VS INFORMED 2-CMT PK MODEL

pop_median_concentrations_1cmt_dose1 <- sapply(1:length(time_sampling), function(t) median(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[1]][,t]))
pop_median_concentrations_2cmt_dose1 <- sapply(1:length(time_sampling), function(t) median(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[1]][,t]))
pop_P5_concentrations_1cmt_dose1 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[1]][,t], probs = 0.05))
pop_P5_concentrations_2cmt_dose1 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[1]][,t], probs = 0.05))
pop_P95_concentrations_1cmt_dose1 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[1]][,t], probs = 0.95))
pop_P95_concentrations_2cmt_dose1 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[1]][,t], probs = 0.95))
plot(x = time_sampling, y = pop_median_concentrations_1cmt_dose1, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curve from the one-compartment PK model for dose 1 obtained by linking the samplings (5-percentile - median - 95-percentile)", col = "red", ylim = c(0,0.65))
lines(x = time_sampling, y = pop_P5_concentrations_1cmt_dose1, type = "l", lty = 2, col = "red")
lines(x = time_sampling, y = pop_P95_concentrations_1cmt_dose1, type = "l", lty = 2, col = "red")
plot(x = time_sampling, y = pop_median_concentrations_2cmt_dose1, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curve from the two-compartment PK model for dose 1 obtained by linking the samplings (5-percentile - median - 95-percentile)", col = "chocolate", ylim = c(0,0.65))
lines(x = time_sampling, y = pop_P5_concentrations_2cmt_dose1, type = "l", lty = 2, col = "chocolate")
lines(x = time_sampling, y = pop_P95_concentrations_2cmt_dose1, type = "l", lty = 2, col = "chocolate")

pop_median_concentrations_1cmt_dose4 <- sapply(1:length(time_sampling), function(t) median(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[4]][,t]))
pop_median_concentrations_2cmt_dose4 <- sapply(1:length(time_sampling), function(t) median(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[4]][,t]))
pop_P5_concentrations_1cmt_dose4 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[4]][,t], probs = 0.05))
pop_P5_concentrations_2cmt_dose4 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[4]][,t], probs = 0.05))
pop_P95_concentrations_1cmt_dose4 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[4]][,t], probs = 0.95))
pop_P95_concentrations_2cmt_dose4 <- sapply(1:length(time_sampling), function(t) quantile(simulated_PK_data_2cmt_informed[[1]]$simu_exact_conc[[4]][,t], probs = 0.95))
plot(x = time_sampling, y = pop_median_concentrations_1cmt_dose4, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curve from the one-compartment PK model for dose 4 obtained by linking the samplings (5-percentile - median - 95-percentile)", col = "red", ylim = c(0,3.2))
lines(x = time_sampling, y = pop_P5_concentrations_1cmt_dose4, type = "l", lty = 2, col = "red")
lines(x = time_sampling, y = pop_P95_concentrations_1cmt_dose4, type = "l", lty = 2, col = "red")
plot(x = time_sampling, y = pop_median_concentrations_2cmt_dose4, type = "l", ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curve from the two-compartment PK model for dose 4 obtained by linking the samplings (5-percentile - median - 95-percentile)", col = "chocolate", ylim = c(0,3.2))
lines(x = time_sampling, y = pop_P5_concentrations_2cmt_dose4, type = "l", lty = 2, col = "chocolate")
lines(x = time_sampling, y = pop_P95_concentrations_2cmt_dose4, type = "l", lty = 2, col = "chocolate")



#####
#####
#####

### Sensitivity analysis for the choice of parameters (Vm, Km) for Michaelis-Menten elimination kinetics implementation in a 1-cmt oral 
### PK model using concentration plots to compare with standard linear elimination kinetics

# First order absorption 1-cmt PK model with Michaelis Menten elimination
# Vm = maximum eliminaton rate (in amount per time unit)
# Km = Michaelis-Menten constant (in concentration unit)
Michaelis_Menten_equation <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dC <- -(((Vm/V)*C)/(Km+C)) + (dose/V)*ka*exp(-ka*t)
    dAUC <- C
    return(list(c(dC, dAUC)))
  })
}
res_MM_conc_per_dose <- matrix(sapply(1:length(doses), function(i) {
  parameters <- c(dose = doses[i], ka = PK_parameters[1], Vm = PK_parameters[2], V = PK_parameters[3], Km = 2)
  state <- c(C = 0, AUC = 0)
  return(ode(y = state, times = time_sampling, func = Michaelis_Menten_equation, parms = parameters, method = "rk4")[, "C"])
}),
ncol = length(doses), nrow = length(time_sampling), byrow = FALSE)
res_MM_AUC_per_dose <- matrix(sapply(1:length(doses), function(i) {
  parameters <- c(dose = doses[i], ka = PK_parameters[1], Vm = PK_parameters[2], V = PK_parameters[3], Km = 2)
  state <- c(C = 0, AUC = 0)
  return(ode(y = state, times = (0:1000), func = Michaelis_Menten_equation, parms = parameters, method = "rk4")[, "AUC"])
}),
ncol = length(doses), nrow = length(0:1000), byrow = FALSE)

pop_mean_concentrations_1cmt_dose1 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[1]][,t]))
pop_mean_concentrations_1cmt_dose2 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[2]][,t]))
pop_mean_concentrations_1cmt_dose3 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[3]][,t]))
pop_mean_concentrations_1cmt_dose4 <- sapply(1:length(time_sampling), function(t) mean(simulated_PK_data_1cmt[[1]]$simu_exact_conc[[4]][,t]))
plot(x = time_sampling, y = pop_mean_concentrations_1cmt_dose1, type = "l", lwd = 3, ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curve from a one-compartment PK model with linear elimination kinetics for all doses", col = "purple", ylim = c(0, 1.4))
lines(x = time_sampling, y = pop_mean_concentrations_1cmt_dose2, type = "l", lwd = 3, col = "brown")
lines(x = time_sampling, y = pop_mean_concentrations_1cmt_dose3, type = "l", lwd = 3, col = "orange")
lines(x = time_sampling, y = pop_mean_concentrations_1cmt_dose4, type = "l", lwd = 3, col = "lightblue")
legend(20, 1.2, legend = c("Dose 1", "Dose 2", "Dose 3", "Dose 4"), col = c("purple", "brown", "orange", "lightblue"), lty = 1, lwd = 3, cex = 0.8)

plot(x = time_sampling, y = res_MM_conc_per_dose[, 1], type = "l", lwd = 3, ylab = "Concentration (mg/L)", xlab = "Time (hours)", main = "Concentration curve from a one-compartment PK model with Michaelis Menten elimination kinetics for all doses", col = "purple", ylim = c(0, 1.4))
lines(x = time_sampling, y = res_MM_conc_per_dose[, 2], type = "l", lwd = 3, col = "brown")
lines(x = time_sampling, y = res_MM_conc_per_dose[, 3], type = "l", lwd = 3, col = "orange")
lines(x = time_sampling, y = res_MM_conc_per_dose[, 4], type = "l", lwd = 3, col = "lightblue")
legend(20, 1.2, legend = c("Dose 1", "Dose 2", "Dose 3", "Dose 4"), col = c("purple", "brown", "orange", "lightblue"), lty = 1, lwd = 3, cex = 0.8)

# # SENSITIVITY ANALYSIS FOR THE CHOICE OF AUC THRESHOLD USING 1-CMT ORAL PK MODEL WITH
# MICHAEL MENTIS ELIMINATION ACCORDING TO THE POSITION OF THE MTD

N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
compartmental_model <- "1-cmt"
ka <- 2
V_pop <- 100
Vm_pop <- 10
Km_pop <- 2
PK_parameters <- c(ka, V_pop, Vm_pop, Km_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
exposure_metric <- "AUC" 
AUC_method <- "integration"
elimination_kinetics <- "Michaelis-Menten"

threshold_AUC_dose_PKdata_MM <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 10, # Random threshold, doesn't matter for AUC histograms
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model,
  elimination_kinetics = elimination_kinetics
)

auc_per_dose_MM <- threshold_AUC_dose_PKdata_MM[[1]]$simu_exact_exposure
hist(auc_per_dose_MM[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "AUC", main = "Michaelis-Menten elimination - Histogram of AUC of all simulated patients for dose 1")
lines(density(auc_per_dose_MM[,1]), lwd = 2, col = "chocolate")
hist(auc_per_dose_MM[,2], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "AUC", main = "Michaelis-Menten elimination - Histogram of AUC of all simulated patients for dose 2")
lines(density(auc_per_dose_MM[,2]), lwd = 2, col = "chocolate")
hist(auc_per_dose_MM[,3], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "AUC", main = "Michaelis-Menten elimination - Histogram of AUC of all simulated patients for dose 3")
lines(density(auc_per_dose_MM[,3]), lwd = 2, col = "chocolate")
hist(auc_per_dose_MM[,4], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "AUC", main = "Michaelis-Menten elimination - Histogram of AUC of all simulated patients for dose 4")
lines(density(auc_per_dose_MM[,4]), lwd = 2, col = "chocolate")
# Threshold at first glance (dose-dependent): (8-11), (15-25), (40-70), (80-160)

# SENSITIVITY ANALYSIS OF AUC THRESHOLD ACCORDING TO DLT PROBABILITY FOR ORAL 1-CMT PK MODEL WITH MICHAELIS-MENTEN ELIMINATION

# MTD = Dose 1 / Threshold = 11.5 -> p_tox = (0.2507, 0.5145, 0.8200, 0.9476)

threshold_AUC_DLT_dose1_PKdata_MM <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 11.5,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model,
  elimination_kinetics = elimination_kinetics
)
auc_per_dose_MM <- threshold_AUC_DLT_dose1_PKdata_MM[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose1_PKdata_MM[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function(i) auc_per_dose_MM[which(tox[,i] == 1),i])

prior_tox_real_dose1 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose_MM[,i]))
prior_tox_real_dose1


# MTD = Dose 2 / Threshold = 20 -> p_tox = (0.0850, 0.2540, 0.5920, 0.8289) 

threshold_AUC_DLT_dose2_PKdata_MM <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 20, #28
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model,
  elimination_kinetics = elimination_kinetics
)
auc_per_dose_MM <- threshold_AUC_DLT_dose2_PKdata_MM[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose2_PKdata_MM[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function(i) auc_per_dose_MM[which(tox[,i] == 1),i])

prior_tox_real_dose2 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose_MM[,i]))
prior_tox_real_dose2


# MTD = Dose 3 / Threshold = 41 -> p_tox = (0.0113, 0.0584, 0.2542, 0.5296)

threshold_AUC_DLT_dose3_PKdata_MM <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 41,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model,
  elimination_kinetics = elimination_kinetics
)
auc_per_dose_MM <- threshold_AUC_DLT_dose3_PKdata_MM[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose3_PKdata_MM[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function(i) auc_per_dose_MM[which(tox[,i] == 1),i])

prior_tox_real_dose3 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose_MM[,i]))
prior_tox_real_dose3


# MTD = Dose 4 / Threshold = 75 -> p_tox = (0.0008, 0.0108, 0.0789, 0.2510)

threshold_AUC_DLT_dose4_PKdata_MM <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 75,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model,
  elimination_kinetics = elimination_kinetics
)
auc_per_dose_MM <- threshold_AUC_DLT_dose4_PKdata_MM[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose4_PKdata_MM[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function(i) auc_per_dose_MM[which(tox[,i] == 1),i])

prior_tox_real_dose4 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose_MM[,i]))
prior_tox_real_dose4


# MTD = Dose 0 / Threshold = 8 -> p_tox = (0.4168, 0.6878, 0.9140, 0.9788)

threshold_AUC_DLT_dose0_PKdata_MM <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 8,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = exposure_metric,
  AUC_method = AUC_method,
  compartmental_model = compartmental_model,
  elimination_kinetics = elimination_kinetics
)
auc_per_dose_MM <- threshold_AUC_DLT_dose0_PKdata_MM[[1]]$simu_exact_exposure
tox <- threshold_AUC_DLT_dose0_PKdata_MM[[1]]$toxicity
auc_with_DLT_per_dose <- lapply(1:length(doses), function(i) auc_per_dose_MM[which(tox[,i] == 1),i])

prior_tox_real_dose0 <- sapply(1:length(doses), function(i) length(auc_with_DLT_per_dose[[i]])/length(auc_per_dose_MM[,i]))
prior_tox_real_dose0



#####
#####
#####

# SENSITIVITY ANALYSIS FOR THE CHOICE OF TAUC THRESHOLD (Truncated AUC) ACCORDING TO THE SELECTED MTD COMPARED WITH 
# CURRENT AUC THRESHOLD USING HISTOGRAMS

N_simu <- 10000
threshold_TAUC_dose_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 6.6,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1
)
conc_per_dose <- list()
tauc_per_dose <- matrix(NA, ncol = length(doses), nrow = N)
Cl_ind <- threshold_TAUC_dose_PKdata[[1]]$individual_parameters[,2:(N+1)]
V_ind <- threshold_TAUC_dose_PKdata[[1]]$individual_parameters[,(N+2):(2*N+1)]
for (i in 1:length(doses)){
  conc_per_dose[[i]] <- threshold_TAUC_dose_PKdata[[1]]$simu_exact_conc[[i]]
  for (j in 1:N){
    integrand <- function(t) {
      (doses[i] / V_ind[j]) * (ka / (ka - (Cl_ind[j] / V_ind[j]))) * (exp(-(Cl_ind[j] / V_ind[j]) * t) - exp(-ka * t))
      }
    tauc_per_dose[j, i] <- integrate(integrand, lower = time_sampling[1], upper = tail(time_sampling, n = 1))$value
  }
}

hist(tauc_per_dose[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "TAUC", main = "Histogramme de TAUC de tous les patients simulés pour la dose 1")
lines(density(tauc_per_dose[,1]), lwd = 2, col = "chocolate")
hist(tauc_per_dose[,2], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "TAUC", main = "Histogramme de TAUC de tous les patients simulés pour la dose 2")
lines(density(tauc_per_dose[,2]), lwd = 2, col = "chocolate")
hist(tauc_per_dose[,3], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "TAUC", main = "Histogramme de TAUC de tous les patients simulés pour la dose 3")
lines(density(tauc_per_dose[,3]), lwd = 2, col = "chocolate")
hist(tauc_per_dose[,4], breaks = 50, col="peachpuff", border = "black", prob = TRUE, xlab = "TAUC", main = "Histogramme de TAUC de tous les patients simuléd pour la dose 4")
lines(density(tauc_per_dose[,4]), lwd = 2, col = "chocolate")

# Threshold at first glance (dose-dependent): (3-5), (6-8), (9-11), (9-11)

# SENSITIVITY ANALYSIS OF TAUC THRESHOLD ACCORDING TO DLT PROBABILITY FOR 1-CMT PK MODEL

N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
compartmental_model <- "1-cmt"
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20

# MTD = Dose 1 / Threshold = 4.6 -> p_tox = (0.2518, 0.4409, 0.6828, 0.8315)

threshold_TAUC_DLT_dose1_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 4.6,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC"
)
tauc_per_dose <- threshold_TAUC_DLT_dose1_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose1_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose1 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose1


# MTD = Dose 2 / Threshold = 7.6 -> p_tox = (0.1230, 0.2524, 0.4819, 0.6715) 

threshold_TAUC_DLT_dose2_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 7.6,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC"
)
tauc_per_dose <- threshold_TAUC_DLT_dose2_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose2_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose2 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose2


# MTD = Dose 3 / Threshold = 14 -> p_tox = (0.0348, 0.1023, 0.2537, 0.4294)

threshold_TAUC_DLT_dose3_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 14,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC"
)
tauc_per_dose <- threshold_TAUC_DLT_dose3_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose3_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose3 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose3


# MTD = Dose 4 / Threshold = 22.5 -> p_tox = (0.0094, 0.0378, 0.1301, 0.2531)

threshold_TAUC_DLT_dose4_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 22.5,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC"
)
tauc_per_dose <- threshold_TAUC_DLT_dose4_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose4_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose4 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose4


# MTD = Dose 0 / Threshold = 3 -> p_tox = (0.4105, 0.6104, 0.8187, 0.9140)

threshold_TAUC_DLT_dose0_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 3,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC"
)
tauc_per_dose <- threshold_TAUC_DLT_dose0_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose0_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose0 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose0



# SENSITIVITY ANALYSIS OF TAUC THRESHOLD ACCORDING TO DLT PROBABILITY FOR 2-CMT PK MODEL

N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
compartmental_model <- "2-cmt"
ka <- 2
Cl_pop <- 10
V_1_pop <- 100
Q_pop <- 5
V_2_pop <- 90
PK_parameters <- c(ka, Cl_pop, V_1_pop, Q_pop, V_2_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20

# MTD = Dose 1 / Threshold = 3.7 -> p_tox = (0.2533, 0.4613, 0.7139, 0.8536)

threshold_TAUC_DLT_dose1_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 3.7,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC",
  compartmental_model = "2-cmt"
)
tauc_per_dose <- threshold_TAUC_DLT_dose1_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose1_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose1 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose1


# MTD = Dose 2 / Threshold = 6.1 -> p_tox = (0.1135, 0.2553, 0.5093, 0.7042) 

threshold_TAUC_DLT_dose2_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 6.1,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC",
  compartmental_model = "2-cmt"
)
tauc_per_dose <- threshold_TAUC_DLT_dose2_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose2_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose2 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose2


# MTD = Dose 3 / Threshold = 11.3 -> p_tox = (0.0297, 0.0915, 0.2540, 0.4478)

threshold_TAUC_DLT_dose3_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 11.3,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC",
  compartmental_model = "2-cmt"
)
tauc_per_dose <- threshold_TAUC_DLT_dose3_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose3_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose3 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose3


# MTD = Dose 4 / Threshold = 18.1 -> p_tox = (0.0089, 0.0321, 0.1194, 0.2551)

threshold_TAUC_DLT_dose4_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 18.1,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC",
  compartmental_model = "2-cmt"
)
tauc_per_dose <- threshold_TAUC_DLT_dose4_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose4_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose4 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose4


# MTD = Dose 0 / Threshold = 2.5 -> p_tox = (0.4104, 0.6305, 0.8337, 0.9285)

threshold_TAUC_DLT_dose0_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = 2.5,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "TAUC",
  compartmental_model = "2-cmt"
)
tauc_per_dose <- threshold_TAUC_DLT_dose0_PKdata[[1]]$simu_exact_exposure
tox <- threshold_TAUC_DLT_dose0_PKdata[[1]]$toxicity
tauc_with_DLT_per_dose <- list()
for (i in 1:length(doses)){
  tauc_with_DLT_per_dose[[i]] <- tauc_per_dose[which(tox[,i] == 1),i]
}

prior_tox_real_dose0 <- sapply(1:length(doses), function(i) length(tauc_with_DLT_per_dose[[i]])/length(tauc_per_dose[,i]))
prior_tox_real_dose0



#####
#####
#####

# SENSIBILITY ANALYSIS FOR PRIOR DISTRIBUTIONS ASSIGNED TO INFERENCE PARAMETERS IN THE DIFFERENT PROSPECTIVE (PK) DOSE-FINDING 
# BAYESIAN DESIGNS (HISTOGRAMS, DENSITY PLOTS, ...)

# BLRM
N_simu <- 10000
# STANDARD PRIOR CHOICE : Log_alpha ~ Normal(-1.098612, 2) / Log_beta ~ Normal(0, 1)
log_alpha_prior_draws <- rnorm(N_simu, f_logit(targeted_tox), 2)
hist(log_alpha_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "alpha", main = "Histogram of alphas drawn from their prior distribution + Density of alphas")
lines(density(log_alpha_prior_draws), lwd = 2, col = "chocolate")
log_beta_prior_draws <- rnorm(N_simu, 0, 1)
hist(log_beta_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta", main = "Histogram of betas drawn from their prior distribution + Density of betas")
lines(density(log_beta_prior_draws), lwd = 2, col = "chocolate")
prior_probabability_toxicity <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)) {
  for (j in 1:N_simu) {
    prior_probabability_toxicity[j,i] <- f_inv_logit(log_alpha_prior_draws[j] + exp(log_beta_prior_draws[j])*log(doses[i]/ref_dose))
  }
}

mean(prior_probabability_toxicity[,1])
mean(prior_probabability_toxicity[,2])
mean(prior_probabability_toxicity[,3])
mean(prior_probabability_toxicity[,4])

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
# 2 components Beta mixture has better looking fit than univariate beta distribution for dose 4
#beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 2, constrain_gt1 = FALSE)
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita")
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita")
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita")
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita")
ESS_dose4

# PRIOR CHOICE #1 : Log_alpha ~ Normal(-1.098612, 0.8) / Log_beta ~ Normal(0, 1)
log_alpha_prior_draws <- rnorm(N_simu, f_logit(targeted_tox), 0.8)
hist(log_alpha_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "alpha", main = "Histogram of alphas drawn from their prior distribution + Density of alphas")
lines(density(log_alpha_prior_draws), lwd = 2, col = "chocolate")
log_beta_prior_draws <- rnorm(N_simu, 0, 1)
hist(log_beta_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta", main = "Histogram of betas drawn from their prior distribution + Density of betas")
lines(density(log_beta_prior_draws), lwd = 2, col = "chocolate")
prior_probabability_toxicity <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)) {
  for (j in 1:N_simu) {
    prior_probabability_toxicity[j,i] <- f_inv_logit(log_alpha_prior_draws[j] + exp(log_beta_prior_draws[j])*log(doses[i]/ref_dose))
  }
}

mean(prior_probabability_toxicity[,1])
mean(prior_probabability_toxicity[,2])
mean(prior_probabability_toxicity[,3])
mean(prior_probabability_toxicity[,4])

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
# 2 components Beta mixture has better looking fit than univariate beta distribution for dose 4
#beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 2, constrain_gt1 = FALSE)
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita")
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita")
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita")
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita")
ESS_dose4


# PRIOR CHOICE #2 : Log_alpha ~ Normal(-1.098612, 0.4) / Log_beta ~ Normal(0.6, 1)
log_alpha_prior_draws <- rnorm(N_simu, f_logit(targeted_tox), 0.4)
hist(log_alpha_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "alpha", main = "Histogram of alphas drawn from their prior distribution + Density of alphas")
lines(density(log_alpha_prior_draws), lwd = 2, col = "chocolate")
log_beta_prior_draws <- rnorm(N_simu, 0.6, 1)
hist(log_beta_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta", main = "Histogram of betas drawn from their prior distribution + Density of betas")
lines(density(log_beta_prior_draws), lwd = 2, col = "chocolate")
prior_probabability_toxicity <- matrix(NA, ncol = length(doses), nrow = N_simu)
for (i in 1:length(doses)) {
  for (j in 1:N_simu) {
    prior_probabability_toxicity[j,i] <- f_inv_logit(log_alpha_prior_draws[j] + exp(log_beta_prior_draws[j])*log(doses[i]/ref_dose))
  }
}

mean(prior_probabability_toxicity[,1])
mean(prior_probabability_toxicity[,2])
mean(prior_probabability_toxicity[,3])
mean(prior_probabability_toxicity[,4])

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
# 2 components Beta mixture has better looking fit than univariate beta distribution for dose 4
#beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 2, constrain_gt1 = FALSE)
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (BLRM)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita")
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita")
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita")
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita")
ESS_dose4



N_simu <- 10000
# PKLOGIT - Combination of 2 Bayesian model to represent a dose-exposure-toxicity relationship
set.seed(12345)
doses <- c(30.6, 50.69, 93.69, 150.37)
exposure_tox_threshold <- 18.2
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
prior_sensitivity_analysis_PKLOGIT_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)
simulated_log_auc <- log(prior_sensitivity_analysis_PKLOGIT_PKdata[[1]]$simu_exact_exposure)
simulated_varia_log_auc <- log(prior_sensitivity_analysis_PKLOGIT_PKdata[[1]]$sensitivity_auc)

hist(simulated_log_auc[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram of exact values for log_AUC for dose 1")
lines(density(simulated_log_auc[,1]), lwd = 2, col = "chocolate")
hist(simulated_log_auc[,2], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram of exact values for log_AUC for dose 2")
lines(density(simulated_log_auc[,2]), lwd = 2, col = "chocolate")
hist(simulated_log_auc[,3], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram of exact values for log_AUC for dose 3")
lines(density(simulated_log_auc[,3]), lwd = 2, col = "chocolate")
hist(simulated_log_auc[,4], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram of exact values for log_AUC for dose 4")
lines(density(simulated_log_auc[,4]), lwd = 2, col = "chocolate")

mean(simulated_log_auc[,1]) #-log(Cl_pop) + log(doses[1])
mean(simulated_log_auc[,2]) #-log(Cl_pop) + log(doses[2])
mean(simulated_log_auc[,3]) #-log(Cl_pop) + log(doses[3])
mean(simulated_log_auc[,4]) #-log(Cl_pop) + log(doses[4])

fixed_z_ij <- sapply(1:length(doses), function (i) rnorm(N_simu, -log(Cl_pop) + log(doses[i]), 1))
ref_log_auc <- mean(simulated_log_auc[,3])

# STANDARD PRIOR CHOICE FOR PK-LOGIT (See Ursino, 2017) - DOSE-AUC MODEL

nu_prior_draws <- rbeta(N_simu, 1, 1) # ESS = 2 - nu ~ Beta(1, 0.4285714) (WEAKLY INFORMATIVE) or Beta(2.33, 1) (MODERATELY INFORMATIVE) or Beta(10/3, 2) (HIGHLY INFORMATIVE) to have on average nu equal 0.7 (omega_IIV) while fixing alpha to be equal to 1
hist(nu_prior_draws, col="peachpuff", border = "black", probability = TRUE, xlab = "nu", main = "Histogram of nu drawn from their prior distribution")

### EVALUATION OF TARGET PRIOR MARGINAL DISTRIBUTION AS THE PRIOR MEAN FOR BETA USING THE CHAIN DRAW FROM THE DISTRIBUTION OF NU
beta0_mean <- -log(Cl_pop)
beta01_sd <- sqrt(1000)
beta0_prior_nu <- sapply(1:N_simu, function (i) rnorm(1, beta0_mean, beta01_sd*nu_prior_draws[i]))
beta1_prior_nu <- sapply(1:N_simu, function (i) rnorm(1, 1, beta01_sd*nu_prior_draws[i]))
z_ij_nu <- lapply(1:length(doses), function(i) {
  sapply(1:N_simu, function(k) rnorm(1, beta0_prior_nu[k] + beta1_prior_nu[k]*log(doses[i]), nu_prior_draws[k]))
})
norm_fitting_zi1 <- mixfit(z_ij_nu[[1]], type = "norm", Nc = 1)
hist(z_ij_nu[[1]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 1")
lines(y = dmix(norm_fitting_zi1, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi2 <- mixfit(z_ij_nu[[2]], type = "norm", Nc = 1)
hist(z_ij_nu[[2]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 2")
lines(y = dmix(norm_fitting_zi2, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi3 <- mixfit(z_ij_nu[[3]], type = "norm", Nc = 1)
hist(z_ij_nu[[3]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 3")
lines(y = dmix(norm_fitting_zi3, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi4 <- mixfit(z_ij_nu[[4]], type = "norm", Nc = 1)
hist(z_ij_nu[[4]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 4")
lines(y = dmix(norm_fitting_zi4, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")

ESS_dose1 <- ess(norm_fitting_zi1, method = "morita", sigma = posterior_mean_nu)
ESS_dose1
ESS_dose2 <- ess(norm_fitting_zi2, method = "morita", sigma = posterior_mean_nu)
ESS_dose2
ESS_dose3 <- ess(norm_fitting_zi3, method = "morita", sigma = posterior_mean_nu)
ESS_dose3
ESS_dose4 <- ess(norm_fitting_zi4, method = "morita", sigma = posterior_mean_nu)
ESS_dose4

### EVALUATION OF PRIOR DISTRIBUTION INDEPENDENTLY OF THE HYPERPRIOR LAYER FOR THE NU PARAMETER
beta0_prior_draws <- rnorm(N_simu, -log(Cl_pop), beta01_sd *omega_IIV)
hist(beta0_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta0", main = "Histogram + Density of beta0 drawn from the prior distribution")
#hist(beta0_prior_draws2, breaks = 50, probability = TRUE, lwd = 2, lty = 1, col = "red", add = T, ylim = c(0, 0.05))
beta1_prior_draws <- rnorm(N_simu, 1, beta01_sd*omega_IIV)
hist(beta1_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta0", main = "Histogram + Density of beta0 drawn from the prior distribution")
#hist(beta1_prior_draws2, breaks = 50, probability = TRUE, lwd = 2, lty = 1, col = "red", add = T, ylim = c(0, 0.05))

z_ij <- lapply(1:length(doses), function(i) {
  sapply(1:N_simu, function(k) rnorm(1, beta0_prior_draws[k] + beta1_prior_draws[k]*log(doses[i]), omega_IIV))
  })
norm_fitting_zi1 <- mixfit(z_ij[[1]], type = "norm", Nc = 1)
hist(z_ij[[1]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 1")
lines(y = dmix(norm_fitting_zi1, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi2 <- mixfit(z_ij[[2]], type = "norm", Nc = 1)
hist(z_ij[[2]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 2")
lines(y = dmix(norm_fitting_zi2, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi3 <- mixfit(z_ij[[3]], type = "norm", Nc = 1)
hist(z_ij[[3]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 3")
lines(y = dmix(norm_fitting_zi3, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi4 <- mixfit(z_ij[[4]], type = "norm", Nc = 1)
hist(z_ij[[4]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 4")
lines(y = dmix(norm_fitting_zi4, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")

ESS_dose1 <- ess(norm_fitting_zi1, method = "morita", sigma = posterior_mean_nu)
ESS_dose1
ESS_dose2 <- ess(norm_fitting_zi2, method = "morita", sigma = posterior_mean_nu)
ESS_dose2
ESS_dose3 <- ess(norm_fitting_zi3, method = "morita", sigma = posterior_mean_nu)
ESS_dose3
ESS_dose4 <- ess(norm_fitting_zi4, method = "morita", sigma = posterior_mean_nu)
ESS_dose4


# CHOICE #1 FOR PK-LOGIT - DOSE-AUC MODEL
nu_prior_draws <- rbeta(N_simu, 1, 1) # ESS = 2 - nu ~ Beta(1, 0.4285714) (WEAKLY INFORMATIVE) or Beta(2.33, 1) (MODERATELY INFORMATIVE) or Beta(10/3, 2) (HIGHLY INFORMATIVE) to have on average nu equal 0.7 (omega_IIV) while fixing alpha to be equal to 1
hist(nu_prior_draws, col="peachpuff", border = "black", probability = TRUE, xlab = "nu", main = "Histogram of nu drawn from their prior distribution")

### EVALUATION OF TARGET PRIOR MARGINAL DISTRIBUTION AS THE PRIOR MEAN FOR BETA USING THE CHAIN DRAW FROM THE DISTRIBUTION OF NU
beta0_mean <- -log(Cl_pop)
beta01_sd <- 10
beta0_prior_nu <- sapply(1:N_simu, function (i) rnorm(1, beta0_mean, beta01_sd*nu_prior_draws[i]))
beta1_prior_nu <- sapply(1:N_simu, function (i) rnorm(1, 1, beta01_sd*nu_prior_draws[i]))
z_ij_nu <- lapply(1:length(doses), function(i) {
  sapply(1:N_simu, function(k) rnorm(1, beta0_prior_nu[k] + beta1_prior_nu[k]*log(doses[i]), nu_prior_draws[k]))
})
norm_fitting_zi1 <- mixfit(z_ij_nu[[1]], type = "norm", Nc = 1)
hist(z_ij_nu[[1]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 1")
lines(y = dmix(norm_fitting_zi1, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi2 <- mixfit(z_ij_nu[[2]], type = "norm", Nc = 1)
hist(z_ij_nu[[2]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 2")
lines(y = dmix(norm_fitting_zi2, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi3 <- mixfit(z_ij_nu[[3]], type = "norm", Nc = 1)
hist(z_ij_nu[[3]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 3")
lines(y = dmix(norm_fitting_zi3, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi4 <- mixfit(z_ij_nu[[4]], type = "norm", Nc = 1)
hist(z_ij_nu[[4]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 4")
lines(y = dmix(norm_fitting_zi4, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")

ESS_dose1 <- ess(norm_fitting_zi1, method = "morita", sigma = omega_IIV)
ESS_dose1
ESS_dose2 <- ess(norm_fitting_zi2, method = "morita", sigma = omega_IIV)
ESS_dose2
ESS_dose3 <- ess(norm_fitting_zi3, method = "morita", sigma = omega_IIV)
ESS_dose3
ESS_dose4 <- ess(norm_fitting_zi4, method = "morita", sigma = omega_IIV)
ESS_dose4

### EVALUATION OF PRIOR DISTRIBUTION INDEPENDENTLY OF THE HYPERPRIOR LAYER FOR THE NU PARAMETER
beta0_prior_draws <- rnorm(N_simu, -log(Cl_pop), beta01_sd *omega_IIV)
hist(beta0_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta0", main = "Histogram + Density of beta0 drawn from the prior distribution")
#hist(beta0_prior_draws2, breaks = 50, probability = TRUE, lwd = 2, lty = 1, col = "red", add = T, ylim = c(0, 0.05))
beta1_prior_draws <- rnorm(N_simu, 1, beta01_sd*omega_IIV)
hist(beta1_prior_draws, breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta0", main = "Histogram + Density of beta0 drawn from the prior distribution")
#hist(beta1_prior_draws2, breaks = 50, probability = TRUE, lwd = 2, lty = 1, col = "red", add = T, ylim = c(0, 0.05))

z_ij <- lapply(1:length(doses), function(i) {
  sapply(1:N_simu, function(k) rnorm(1, beta0_prior_draws[k] + beta1_prior_draws[k]*log(doses[i]), omega_IIV))
})
norm_fitting_zi1 <- mixfit(z_ij[[1]], type = "norm", Nc = 1)
hist(z_ij[[1]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 1")
lines(y = dmix(norm_fitting_zi1, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi2 <- mixfit(z_ij[[2]], type = "norm", Nc = 1)
hist(z_ij[[2]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 2")
lines(y = dmix(norm_fitting_zi2, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi3 <- mixfit(z_ij[[3]], type = "norm", Nc = 1)
hist(z_ij[[3]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 3")
lines(y = dmix(norm_fitting_zi3, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")
norm_fitting_zi4 <- mixfit(z_ij[[4]], type = "norm", Nc = 1)
hist(z_ij[[4]], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_AUC", main = "Histogram + Density of log_AUC for dose 4")
lines(y = dmix(norm_fitting_zi4, seq(-1000, 1000, 0.1)), x = seq(-1000, 1000, 0.1), lwd = 2, lty = 1, col = "red")

ESS_dose1 <- ess(norm_fitting_zi1, method = "morita", sigma = omega_IIV)
ESS_dose1
ESS_dose2 <- ess(norm_fitting_zi2, method = "morita", sigma = omega_IIV)
ESS_dose2
ESS_dose3 <- ess(norm_fitting_zi3, method = "morita", sigma = omega_IIV)
ESS_dose3
ESS_dose4 <- ess(norm_fitting_zi4, method = "morita", sigma = omega_IIV)
ESS_dose4



# STANDARD PRIOR CHOICE FOR PK-LOGIT (See Ursino, 2017) - AUC-TOXICITY MODEL
beta2_prior_draws <- runif(N_simu, 0, 20)
beta3_prior_draws <- runif(N_simu, 0, 10)

simu_sample_logAUC <- matrix(sapply(1:length(doses), function(i) rnorm(N_simu, -log(Cl_pop) + log(doses[i]), omega_IIV)),
                             nrow = N_simu,
                             ncol = length(doses),
                             byrow = FALSE)
simu_ref_log_AUC <- mean(simu_sample_logAUC[, which(doses == ref_dose)])
prior_probabability_toxicity <- matrix(sapply(1:length(doses), function(i) f_inv_logit(-beta2_prior_draws + exp(log_beta3_prior_draws)*(simu_sample_logAUC[,i]-simu_ref_log_AUC))),
                                       nrow = N_simu,
                                       ncol = length(doses),
                                       byrow = FALSE)
mean_prior_probabability_toxicity <- sapply(1:length(doses), function(i) mean(prior_probabability_toxicity[,i]))


mean_prior_probabability_toxicity[1]
mean_prior_probabability_toxicity[2]
mean_prior_probabability_toxicity[3] 
mean_prior_probabability_toxicity[4] 

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita", s = 1000)
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita", s = 1000) 
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita", s = 1000) 
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita", s = 1000)
ESS_dose4


# PRIOR CHOICE BASED ON OPTIM FOR PK-LOGIT - DOSE-TOXICITY MODEL
f_to_be_optim <- function(beta2_logbeta3_means) {
  beta2 <- rnorm(N_simu, beta2_logbeta3_means[1], 5)
  log_beta3 <- rnorm(N_simu, beta2_logbeta3_means[2], 1)
  simu_sample_logAUC <- matrix(sapply(1:length(doses), function(i) rnorm(N_simu, -log(Cl_pop) + log(doses[i]), omega_IIV)),
                                nrow = N_simu,
                                ncol = length(doses),
                                byrow = FALSE)
  simu_ref_log_AUC <- mean(simu_sample_logAUC[, which(doses == ref_dose)])
  prior_probabability_toxicity <- matrix(sapply(1:length(doses), function(i) f_inv_logit(-beta2 + exp(log_beta3)*(simu_sample_logAUC[,i]-simu_ref_log_AUC))),
                                  nrow = N_simu,
                                  ncol = length(doses),
                                  byrow = FALSE)
  mean_prior_probabability_toxicity <- sapply(1:length(doses), function(i) mean(prior_probabability_toxicity[,i]))
  return(abs(mean_prior_probabability_toxicity[2]-0.10) + abs(mean_prior_probabability_toxicity[3]-0.25) + abs(mean_prior_probabability_toxicity[4]-0.40))
}
optim_priors_pklogit <- optim(c(2, 0), f_to_be_optim, method = "BFGS")
print(optim_priors_pklogit)
# 4.392648 1.537353 - CONVERGENCE = 0 - "L-BFGS-B" "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
# 3.638781 -1.666915 - CONVERGENCE = 0 - "BFGS"
# 5.960688 1.967642 - CONVERGENCE = 0 - "BFGS" => FINAL CHOICE

beta2_prior_draws <- rnorm(N_simu, 6, 5)
log_beta3_prior_draws <- rnorm(N_simu, 2, 1) #2.15 => IDEAL CHOICE BUT CAN BE ROUND TO 2
simu_sample_logAUC <- matrix(sapply(1:length(doses), function(i) rnorm(N_simu, -log(Cl_pop) + log(doses[i]), omega_IIV)),
                             nrow = N_simu,
                             ncol = length(doses),
                             byrow = FALSE)
simu_ref_log_AUC <- mean(simu_sample_logAUC[, which(doses == ref_dose)])
prior_probabability_toxicity <- matrix(sapply(1:length(doses), function(i) f_inv_logit(-beta2_prior_draws + exp(log_beta3_prior_draws)*(simu_sample_logAUC[,i]-simu_ref_log_AUC))),
                                nrow = N_simu,
                                ncol = length(doses),
                                byrow = FALSE)
mean_prior_probabability_toxicity <- sapply(1:length(doses), function(i) mean(prior_probabability_toxicity[,i]))

mean_prior_probabability_toxicity[1]
mean_prior_probabability_toxicity[2]
mean_prior_probabability_toxicity[3]
mean_prior_probabability_toxicity[4]

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (PKLOGIT)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita", s = 1000)
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita", s = 1000) 
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita", s = 1000) 
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita", s = 1000)
ESS_dose4



# ED-EWOC / ED (Parameters: BETA=(beta0, beta1)) - ED (See Micallef, 2022)
set.seed(12345)
doses <- c(30.6, 50.69, 93.69, 150.37)
exposure_tox_threshold <- 18.2
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
prior_sensitivity_analysis_ED_EWOC_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)
simulated_auc <- prior_sensitivity_analysis_ED_EWOC_PKdata[[1]]$simu_exact_exposure
simulated_ref_auc <- ref_AUC_estim_ed_ewoc(ref_dose = ref_dose, PK_parameters = PK_parameters, cv = cv, time_sampling = time_sampling)

# STANDARD PRIOR CHOICE FOR ED-EWOC/ED TOXICITY ESTIMATION  (See Micallef, 2022)
beta_paramsmean <- c(-1.808, 0.032)
beta_paramscov <- matrix(c(1.172, 0.142, 0.142, 0.365), 2, 2)
mvt_beta_prior_draws <- mvrnorm(N_simu, beta_paramsmean, beta_paramscov)
hist(mvt_beta_prior_draws[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta_2", main = "Histogram of betas_2 drawn from their prior distribution + Density of betas_2")
lines(density(mvt_beta_prior_draws[,1]), lwd = 2, col = "chocolate")
hist(mvt_beta_prior_draws[,2], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_beta_3", main = "Histogram of log_betas_3 drawn from their prior distribution + Density of log_betas_3")
lines(density(mvt_beta_prior_draws[,2]), lwd = 2, col = "chocolate")

prior_probabability_toxicity <- matrix(NA, ncol = length(doses), nrow = N_simu)
for(i in 1:length(doses)) {
  prior_probabability_toxicity[,i] <- f_inv_logit(mvt_beta_prior_draws[,1] + exp(mvt_beta_prior_draws[,2])*(fixed_z_ij[,i] - ref_log_auc))
}

mean(prior_probabability_toxicity[,1])
mean(prior_probabability_toxicity[,2])
mean(prior_probabability_toxicity[,3])
mean(prior_probabability_toxicity[,4])

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita", s = 1000) # Same as doing alpha + beta but using package function
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita", s = 1000) 
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita", s = 1000) 
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita", s = 1000)
ESS_dose4

# PRIOR CHOICE #1 FOR ED-EWOC / ED
beta_paramsmean <- c(-2, 0) # APPROXIMATION OF f_logit(targeted_tox)-2
beta_paramscov <- matrix(c(5, 1, 1, 2), 2, 2)
mvt_beta_prior_draws <- mvrnorm(N_simu, beta_paramsmean, beta_paramscov)
hist(mvt_beta_prior_draws[,1], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "beta_2", main = "Histogram of betas_2 drawn from their prior distribution + Density of betas_2")
lines(density(mvt_beta_prior_draws[,1]), lwd = 2, col = "chocolate")
hist(mvt_beta_prior_draws[,2], breaks = 50, col="peachpuff", border = "black", probability = TRUE, xlab = "log_beta_3", main = "Histogram of log_betas_3 drawn from their prior distribution + Density of log_betas_3")
lines(density(mvt_beta_prior_draws[,2]), lwd = 2, col = "chocolate")

prior_probabability_toxicity <- matrix(NA, ncol = length(doses), nrow = N_simu)
for(i in 1:length(doses)) {
  prior_probabability_toxicity[,i] <- f_inv_logit(mvt_beta_prior_draws[,1] + exp(mvt_beta_prior_draws[,2])*(fixed_z_ij[,i] - ref_log_auc))
}

mean(prior_probabability_toxicity[,1])
mean(prior_probabability_toxicity[,2])
mean(prior_probabability_toxicity[,3])
mean(prior_probabability_toxicity[,4])

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita", s = 1000) # Same as doing alpha + beta but using package function
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita", s = 1000) 
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita", s = 1000) 
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita", s = 1000)
ESS_dose4


# PRIOR CHOICE #2 FOR ED-EWOC / ED
# Using parameters obtained for PKLOGIT using optimization
beta0_prior_draws <- rnorm(N_simu, -optim_priors_pklogit$par[1], 5)
log_beta1_prior_draws <- rnorm(N_simu, optim_priors_pklogit$par[2], 1)
#beta_paramsmean <- c(-4.559307, 1.703042) 
#beta_paramscov <- matrix(c(5, 0, 0, 1), 2, 2)
#mvt_beta_prior_draws <- mvrnorm(N_simu, beta_paramsmean, beta_paramscov)
prior_probabability_toxicity <- matrix(NA, ncol = length(doses), nrow = N_simu)
for(i in 1:length(doses)) {
  prior_probabability_toxicity[,i] <- f_inv_logit(beta0_prior_draws + exp(log_beta1_prior_draws)*(fixed_z_ij[,i] - ref_log_auc))
}

mean(prior_probabability_toxicity[,1])
mean(prior_probabability_toxicity[,2])
mean(prior_probabability_toxicity[,3])
mean(prior_probabability_toxicity[,4])

beta_fitting_dose1 <- mixfit(prior_probabability_toxicity[,1], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,1], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 1 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose1, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose2 <- mixfit(prior_probabability_toxicity[,2], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,2], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 2 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose2, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose3 <- mixfit(prior_probabability_toxicity[,3], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,3], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 3 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose3, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")
beta_fitting_dose4 <- mixfit(prior_probabability_toxicity[,4], type = "beta", Nc = 1, constrain_gt1 = FALSE)
hist(prior_probabability_toxicity[,4], col="peachpuff", border = "black", probability = TRUE, xlab = "Prior probability of toxicity", main = "Histogram of prior toxicity probabilities for dose 4 (ED-EWOC/ED)", xlim = c(0,1))
lines(y = dmix(beta_fitting_dose4, seq(0, 1, 0.001)), x = seq(0, 1, 0.001), lwd = 2, lty = 2, col = "red")

ESS_dose1 <- ess(beta_fitting_dose1, method = "morita", s = 1000) # Same as doing alpha + beta but using package function
ESS_dose1
ESS_dose2 <- ess(beta_fitting_dose2, method = "morita", s = 1000) 
ESS_dose2
ESS_dose3 <- ess(beta_fitting_dose3, method = "morita", s = 1000) 
ESS_dose3
ESS_dose4 <- ess(beta_fitting_dose4, method = "morita", s = 1000)
ESS_dose4



# PK-CRM (Parameters: alpha, beta)
# ELIMINATED FROM PK MODEL SELECTION FOR FURTHER SELECTION USING SIMULATIONS
set.seed(12345)
doses <- c(30.6, 50.69, 93.69, 150.37)
exposure_tox_threshold <- 18.2
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
prior_sensitivity_analysis_PK_CRM_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)
simulated_auc <- prior_sensitivity_analysis_PK_CRM_PKdata[[1]]$simu_exact_exposure
seq_doses <- as.vector(sapply(1:length(doses), function(i) rep(doses[i], N_simu/length(doses))))
seq_doses_num <- as.vector(sapply(1:length(doses), function(i) rep(i, N_simu/length(doses))))
conc_without_error <- prior_sensitivity_analysis_ED_EWOC_PKdata[[1]]$simu_exact_conc
auc_s <- as.vector(sapply(1:N_simu, function(i) prior_sensitivity_analysis_PK_CRM_PKdata[[1]]$simu_exact_exposure[i, seq_doses_num[i]]))

alpha_prior <- 10
beta_prior <- 5
# Performing linear regression model for auc-dose relationship as stated in the model to estimate the parameters (b0, b1)
reg_auc_dose <- lm(auc_s ~ log(doses[x]))
estimated_b0_b1 <- as.vector(reg_auc_dose$coefficients)
# Use the estimation for the parameters of the frequentist linear regression to estimate the AUC for each dose
log_estimated_auc_dose <- estimated_b0_b1[1] + estimated_b0_b1[2]*log(doses) # Arbitrary solution to solve the problem of negative AUC

# Computation of the standardized adjustement AUC (x_ij) at each dose level using the previously estimated AUC
sd_adj_auc <- log_estimated_auc_dose - sum(log_estimated_auc_dose)/length(doses)
sd_adj_auc_ind <- sd_adj_auc[x]


# TITE-PK / NAIVE AND INFORMED VERSIONS (Parameters: beta)
# PRIOR DISTRIBUTION FOR BETA FROM THE PAPER KEPT AS IT IS : GOOD ENOUGH!



#####
#####
#####

# SUMMARY PLOT FOR SCENARIO B1-E5 RESULTS USING 95% CLOPPER-PEARSON CONFIDENCE INTERVALS FOR THE PERCENTAGE OF CORRECT MTD SELECTION IN EACH SCENARIOS FOR ALL METHODS

# Data importation for all selected scenarios (Results for scenarios B1-E5)
scenario_id_start <- 6
scenario_id_end <- 25
set_name_list <- c(rep("B", 5), rep("C", 5), rep("D", 5), rep("E", 5))
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
MTD_id_scenarios <- rep(c(1, 2, 3, 4, 0), 6)
clopper_pearson_95CI_all_scenarios <- list()
clopper_pearson_95CI_plots <- list()
for(i in 1:(length(doses) + 1)){
  standard_scenario_name <- paste0("A", i)
  load(paste0(folder_path, "run_scen/scenario_", standard_scenario_name, ".Rdata"))
  current_set_name <- set_name_list[i]
  current_scenario_num <- ifelse(MTD_id_scenarios[i - scenario_id_start + 1] == 0, 5, MTD_id_scenarios[i - scenario_id_start + 1])
  for(j in res_model_name) {
    ifelse(i == 5, assign(paste0("res_", j, "_scenario", i), length(which(eval(parse(text = paste0("res_", j, "@MTD_each_trial"))) == 0))), assign(paste0("res_", j, "_scenario", i), length(which(eval(parse(text = paste0("res_", j, "@MTD_each_trial"))) == i))))
    assign(paste0("clopper_pearson_95CI_", j, "_scenario", i), BinomCI(x = eval(parse(text = paste0("res_", j, "_scenario", i))), n = n_trials, conf.level = 0.95, sides = "two.sided", method = "clopper-pearson", rand = seed))
    assign(paste0("estim_95CI_", j, "_scenario", i), as.numeric(eval(parse(text = paste0("clopper_pearson_95CI_", j, "_scenario", i)))[,"est"]))
    assign(paste0("lower_95CI_", j, "_scenario", i), as.numeric(eval(parse(text = paste0("clopper_pearson_95CI_", j, "_scenario", i)))[,"lwr.ci"]))
    assign(paste0("upper_95CI_", j, "_scenario", i), as.numeric(eval(parse(text = paste0("clopper_pearson_95CI_", j, "_scenario", i)))[,"upr.ci"]))
  }
}
for (i in (scenario_id_start:scenario_id_end)) {
  if(1 <= i & i <= 5) {
    scenario_name <- paste0("A", i)
  } else if (6 <= i & i <= 10) {
    scenario_name <- paste0("B", i - (length(doses) + 1))
  } else if (11 <= i & i <= 15) {
    scenario_name <- paste0("C", i - 2*(length(doses) + 1))
  } else if (16 <= i & i <= 20) {
    scenario_name <- paste0("D", i - 3*(length(doses) + 1))
  } else if (21 <= i & i <= 25) {
    scenario_name <- paste0("E", i - 4*(length(doses) + 1))
  }
  load(paste0(folder_path, "run_scen/scenario_", scenario_name, ".Rdata"))
  number_correct_MTD_selection <- c()
  estim_95CI <- c()
  lower_95CI <- c()
  upper_95CI <- c()
  scenario_type_95CI <- c()
  current_set_name <- set_name_list[i - scenario_id_start + 1]
  current_scenario_num <- ifelse(MTD_id_scenarios[i - scenario_id_start + 1] == 0, 5, MTD_id_scenarios[i - scenario_id_start + 1])
  for (j in 1:length(res_model_name)) {
    number_correct_MTD_selection <- append(number_correct_MTD_selection, length(which(eval(parse(text = paste0("res_", res_model_name[j], "@MTD_each_trial"))) == MTD_id_scenarios[i - scenario_id_start + 1])))
    clopper_pearson_95CI <- BinomCI(x = number_correct_MTD_selection[j], n = n_trials, conf.level = 0.95, sides = "two.sided", method = "clopper-pearson", rand = seed)
    if (MTD_id_scenarios[i - scenario_id_start + 1] == 0) {
      estim_95CI <- append(estim_95CI, c(as.numeric(clopper_pearson_95CI[,"est"]), eval(parse(text = paste0("estim_95CI_", res_model_name[j], "_scenario5")))))
      lower_95CI <- append(lower_95CI, c(as.numeric(clopper_pearson_95CI[,"lwr.ci"]), eval(parse(text = paste0("lower_95CI_", res_model_name[j], "_scenario5")))))
      upper_95CI <- append(upper_95CI, c(as.numeric(clopper_pearson_95CI[,"upr.ci"]), eval(parse(text = paste0("upper_95CI_", res_model_name[j], "_scenario5")))))
      scenario_type_95CI <- append(scenario_type_95CI, c("Current scenario", "Corresponding scenario in set A"))
    } else {
      estim_95CI <- append(estim_95CI, c(as.numeric(clopper_pearson_95CI[,"est"]), eval(parse(text = paste0("estim_95CI_", res_model_name[j], "_scenario", MTD_id_scenarios[i - scenario_id_start + 1])))))
      lower_95CI <- append(lower_95CI, c(as.numeric(clopper_pearson_95CI[,"lwr.ci"]), eval(parse(text = paste0("lower_95CI_", res_model_name[j], "_scenario", MTD_id_scenarios[i - scenario_id_start + 1])))))
      upper_95CI <- append(upper_95CI, c(as.numeric(clopper_pearson_95CI[,"upr.ci"]), eval(parse(text = paste0("upper_95CI_", res_model_name[j], "_scenario", MTD_id_scenarios[i - scenario_id_start + 1])))))
      scenario_type_95CI <- append(scenario_type_95CI, c("Current scenario", "Corresponding scenario in set A"))
    }
  }
  clopper_pearson_95CI_scenario <- data.frame(model_name = as.vector(sapply(1:length(res_model_name), function(k) rep(res_model_name[k], 2))), scenario_type_95CI = scenario_type_95CI, estim_95CI = estim_95CI, lower_95CI = lower_95CI, upper_95CI = upper_95CI)
  clopper_pearson_95CI_scenario$model_name <- factor(clopper_pearson_95CI_scenario$model_name, levels = res_model_name)
  clopper_pearson_95CI_scenario$scenario_type_95CI <- factor(clopper_pearson_95CI_scenario$scenario_type_95CI, levels = c("Current scenario", "Corresponding scenario in set A"))
  clopper_pearson_95CI_all_scenarios[[paste0("Scenario_", scenario_name, sep = "")]] <- clopper_pearson_95CI_scenario
  clopper_pearson_95CI_plots[[paste0("Scenario_", scenario_name, sep = "")]] <- local({
    ggplot(data = eval(parse(text = paste0("clopper_pearson_95CI_all_scenarios$Scenario_", scenario_name))), aes(x = model_name, y = estim_95CI, family = scenario_type_95CI, color = scenario_type_95CI)) +
      geom_errorbar(aes(ymin = lower_95CI, ymax = upper_95CI), width = 0.2, position = position_dodge(0.8)) + #0.3
      geom_point(position = position_dodge(0.8)) + #0.3
      scale_y_continuous(limits = c(0,1)) +
      scale_x_discrete(labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "N.TITE-PK", "informed_tite_pk" = "Inf.TITE-PK")) +
      scale_color_grey(labels = c("Current scenario" = "Current scenario", "Corresponding scenario in set A" = "Corresponding scenario in set A")) +
      labs(title = paste("Scenario", paste0(current_set_name, current_scenario_num)), y = "% of correct MTD selection", x = "Dose-finding designs", color = "Scenarios") +
      scale_fill_brewer(palette="Dark2") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5, angle = 90), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5), legend.title = element_text(face = "bold", size = 10))
  })
}
setB <- ggplot() + annotate(geom = 'text', x=1, y=1, label="B", fontface = "bold", size = 4) + theme_void() 
setC <- ggplot() + annotate(geom = 'text', x=1, y=1, label="C", fontface = "bold", size = 4) + theme_void()
setD <- ggplot() + annotate(geom = 'text', x=1, y=1, label="D", fontface = "bold", size = 4) + theme_void() 
setE <- ggplot() + annotate(geom = 'text', x=1, y=1, label="E", fontface = "bold", size = 4) + theme_void() 
MTDdose1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MTD = Dose 1", fontface = "bold", size = 4) + theme_void() 
MTDdose2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MTD = Dose 2", fontface = "bold", size = 4) + theme_void()
MTDdose3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MTD = Dose 3", fontface = "bold", size = 4) + theme_void() 
MTDdose4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="MTD = Dose 4", fontface = "bold", size = 4) + theme_void()
noMTD <- ggplot() + annotate(geom = 'text', x=1, y=1, label="No MTD", fontface = "bold", size = 4) + theme_void()
fill <- ggplot() + annotate(geom = 'text', x=1, y=1, label=" ", fontface = "bold", size = 4) + theme_void()

clopper_pearson_95CI_figure <- ggarrange(fill,
                                         MTDdose1,
                                         MTDdose2,
                                         MTDdose3,
                                         MTDdose4,
                                         noMTD,
                                         setB, 
                                         clopper_pearson_95CI_plots$Scenario_B1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_B2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_B3 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_B4 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_B5 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                         setC,
                                         clopper_pearson_95CI_plots$Scenario_C1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_C2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                         clopper_pearson_95CI_plots$Scenario_C3 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_C4 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_C5 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                         setD,
                                         clopper_pearson_95CI_plots$Scenario_D1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"),
                                         clopper_pearson_95CI_plots$Scenario_D2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                         clopper_pearson_95CI_plots$Scenario_D3 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                         clopper_pearson_95CI_plots$Scenario_D4 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                         clopper_pearson_95CI_plots$Scenario_D5 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                         setE,
                                         clopper_pearson_95CI_plots$Scenario_E1, 
                                         clopper_pearson_95CI_plots$Scenario_E2 + rremove("ylab"), 
                                         clopper_pearson_95CI_plots$Scenario_E3 + rremove("ylab"), 
                                         clopper_pearson_95CI_plots$Scenario_E4 + rremove("ylab"), 
                                         clopper_pearson_95CI_plots$Scenario_E5 + rremove("ylab"), 
                                         ncol = 6, nrow = 5, widths = c(0.1, 1.05, 1, 1, 1, 1), heights = c(0.1, 1, 1, 1, 1.25), common.legend = TRUE, legend = "bottom")
ggsave(paste0("plots/clopper_pearson_95CI_figure.png"), plot = wrap_elements(clopper_pearson_95CI_figure), dpi = "retina", units = "mm", width = 210, height = 230)



#####
#####
#####


# Investigation on simulation study results for PKLOGIT method with a simulated clinical trial for 10,000 patients
set.seed(12345)
N_simu <- 10000
doses <- c(30.6, 50.69, 93.69, 150.37)
exposure_tox_threshold <- 18.2
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
PKLOGIT_results_investigation_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)
simulated_log_auc <- log(PKLOGIT_results_investigation_PKdata[[1]]$simu_exact_exposure)
hist(simulated_log_auc, breaks = 50)
simulated_noisy_log_auc <- log(PKLOGIT_results_investigation_PKdata[[1]]$sensitivity_auc)
hist(simulated_noisy_log_auc, breaks = 50)
simulated_log_auc_noncompartmental <- matrix(sapply(1:N_simu, function (i) {
  log_auc_noncompartmental_patient_i <- rep(NA, length(doses))
  for (j in 1:length(doses)) {
    conc_patient_i <- PKLOGIT_results_investigation_PKdata[[1]]$simu_exact_conc[[j]][i,real_sampling]
    log_auc_noncompartmental_patient_i[j] <- log(PK::auc(conc = conc_patient_i, time = time_sampling[real_sampling], design = "complete")$est)
  }
  return(log_auc_noncompartmental_patient_i)
}),nrow = N_simu, ncol = length(doses), byrow = TRUE)
hist(simulated_log_auc_noncompartmental, breaks = 50)
# Small statistical check using Kolgomorov-Smirnov test to verigy if both samples of log(AUC) were drawn from the same continuous distribution
ks.test(x = as.vector(simulated_log_auc_noncompartmental), y = as.vector(simulated_log_auc), alternative = "two.sided")


# Running test clinical trial for 40 patients on PKLOGIT method by giving one time, the true/noisy AUC values obtained with the compartmental
# approach and the other time, the AUC values obtained with the trapezoidal rule using the non-compartmental approach
set.seed(12345)
N_simu <- 40
doses <- c(30.6, 50.69, 93.69, 150.37)
cohort_size <- 2
exposure_tox_threshold <- 6.6
ka <- 2
Cl_pop <- 10
V_pop <- 100
PK_parameters <- c(ka, Cl_pop, V_pop)
omega_IIV <- 0.7
omega_alpha <- 0.8
n_sampling_timepoints <- 10
real_sampling <- c(2:6, round(seq(9, 48, ((48-9)/4))))
cv <- 0.20
prior_tox_real <- pnorm((log(doses) - log(exposure_tox_threshold) - log(PK_parameters[2])) / sqrt(omega_IIV^2 + omega_alpha^2))
PKLOGIT_40Trial_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)
conc <- PKLOGIT_40Trial_PKdata[[1]]$simu_exact_conc
conc_with_error <- PKLOGIT_40Trial_PKdata[[1]]$simu_conc_with_error
exact_auc <- PKLOGIT_40Trial_PKdata[[1]]$simu_exact_exposure
sensitivity_auc <- PKLOGIT_40Trial_PKdata[[1]]$sensitivity_auc
tox <- PKLOGIT_40Trial_PKdata[[1]]$toxicity
Cl_ind <- PKLOGIT_40Trial_PKdata[[1]]$individual_parameters[2:(N+1)]
V_ind <- PKLOGIT_40Trial_PKdata[[1]]$individual_parameters[(N+2):((2*N)+1)]
alpha <- PKLOGIT_40Trial_PKdata[[1]]$alpha

### Assigning the patients of the trial to different doses to have 10 patients receiving each dose and extract the corresponding toxicity outcome
### as well as the AUC corresponding to the administered dose for each patient
x <- c(rep(1, N_simu/length(doses)), rep(4, N_simu/length(doses)), rep(2, N_simu/length(doses)), rep(3, N_simu/length(doses)))
y <- tox[cbind(1:length(x), x)]

# True AUC values with compartmental approach
#auc_s <- sapply(1:N_simu, function(i) exact_auc[i, x[i]])
# Noisy AUC values with compartmental approach
#auc_s <- sapply(1:N_simu, function(i) sensitivity_auc[i, x[i]])
# AUC values computed with the non-compartmental approach using true concentration values
#conc_all <- matrix(sapply(1:N_simu, function(i) as.vector(conc[[x[i]]][i, real_sampling])), nrow = N_simu, ncol = length(real_sampling), byrow = TRUE)
#auc_s <- sapply(1:N_simu, function(i) PK::auc(conc = conc_all[i,], time = time_sampling[real_sampling], design = "complete")$est)
# AUC values computed with the non-compartmental approach using noisy concentration values (data used in clinical trials for simulation study)
conc_all <- matrix(sapply(1:N_simu, function(i) as.vector(conc_with_error[[x[i]]][i, real_sampling])), nrow = N_simu, ncol = length(real_sampling), byrow = TRUE)
auc_s <- sapply(1:N_simu, function(i) PK::auc(conc = conc_all[i,], time = time_sampling[real_sampling], design = "complete")$est)

### Setting up the appropriate parameters to instantiate the PKLOGIT method
beta0_mean <- -log(Cl_pop)
beta01_sd <- 10
beta2_mean <- 6
beta2_sd <- 5
log_beta3_mean <- 2
log_beta3_sd <- 1

### Model estimation using Bayesian inference

# Stan for AUC computation using normal prior distribution
doses_AUC <- cbind(rep(1, N_simu), log(doses[x]))
data_stan <- list(N = N_simu, log_AUC = log(auc_s), doses = doses_AUC, beta0_mean = beta0_mean, 
                  beta01_sd = beta01_sd)
reg_AUC <- sampling(model_PKLOGIT_AUC,
                    data = data_stan, iter = options$n_iter, chains = options$n_chains,
                    control = list(adapt_delta = options$n_adapt),
                    cores = options$n_cores
)
res_AUC <- get_posterior_mean(reg_AUC)
sampl_AUC <- extract(reg_AUC)
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
data_stan <- list(N = N_simu, y = y, log_AUC = log(auc_s), ref_log_AUC = simu_ref_log_AUC, 
                  beta2_mean = beta2_mean, beta2_sd = beta2_sd, 
                  log_beta3_mean = log_beta3_mean, log_beta3_sd = log_beta3_sd)
reg_logit <- sampling(model_PKLOGIT_logit,
                      data = data_stan, iter = options$n_iter, chains = options$n_chains,
                      control = list(adapt_delta = options$n_adapt),
                      cores = options$n_cores
)
res_logit <- get_posterior_mean(reg_logit)
sampl_logit <- extract(reg_logit)
beta2 <- res_logit[1, options$n_chains + 1]
log_beta3 <- res_logit[2, options$n_chains + 1]

# Computation of predictive probability of toxicity to determine the recommended dose to be given to the next cohort OR the MTD using the
# final posterior probabilities of toxicity for the trial if all cohorts have been included
p_estim_dose_tox_MCMC <<- matrix(sapply(1:length(doses), function(i) f_inv_logit(-res_logit[1, options$n_chains + 1] + exp(res_logit[2, options$n_chains + 1])*(simu_sample_logAUC[,i] - simu_ref_log_AUC))),
                                 nrow = n_simu,
                                 ncol = length(doses),
                                 byrow = FALSE)
p_estim_dose_mean <- sapply(1:length(doses), function(i) mean(p_estim_dose_tox_MCMC[, i]))

### The posterior probabilities of toxicity
centered_sample_logAUC <<- matrix(sapply(simu_sample_logAUC, function(i) i - simu_ref_log_AUC), nrow = n_simu, ncol = length(doses), byrow = FALSE)
beta2_mc <- sampl_logit$beta2
log_beta3_mc <- sampl_logit$log_beta3
f_estim_sum <- function(o) {
  p_estim_simu_sample_logAUC <- f_inv_logit(-beta2_mc[o] + exp(log_beta3_mc[o])*centered_sample_logAUC)
  return(colMeans(p_estim_simu_sample_logAUC))
}
p_estim_sum <- matrix(sapply(1:((options$n_chains*options$n_iter)/2), f_estim_sum), nrow = options$n_chains*options$n_iter/2, ncol = length(doses), byrow = TRUE)

### Check the safety rule(s) to see if the trial need to be stopped before allocating the recommended dose to the next cohort
proba_stop <- checking_stop_rule(p_estim_sum[1, ], target = targeted_tox, error = 0)
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



#####
#####
#####

# SUMMARY PLOT FOR SCENARIO 6-25 RESULTS USING BOXPLOTS FOR THE ESTIMATION OF THE PROBABILITY OF 
# TOXICITY FOR DOSE 3 (REFERENCE DOSE) BY EACH DOSE-FINDING METHOD

# Data importation for all selected scenarios (Results for scenarios 6-25)
scenario_id_start <- 6 #2
scenario_id_end <- 25 #26
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
estimated_proba_tox_dose3_all_scenarios <- list()
estimated_proba_tox_dose3_plots <- list()
real_proba_tox_dose3 <- list()
for (i in (scenario_id_start:scenario_id_end)) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_dose3[[paste0("Scenario_", i, sep = "")]] <- data.frame(real_proba_tox_dose3 = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real"))))[3])
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_alldoses_scenario <- data.frame(model_name =  model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_alldoses_scenario$model_name <- factor(estimated_proba_tox_alldoses_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_alldoses_scenario$dose_num <- factor(estimated_proba_tox_alldoses_scenario$dose_num)
  estimated_proba_tox_dose3_scenario <- estimated_proba_tox_alldoses_scenario[estimated_proba_tox_alldoses_scenario$dose_num == "dose3",]
  estimated_proba_tox_dose3_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_dose3_scenario
  estimated_proba_tox_dose3_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_dose3_all_scenarios$Scenario_", i))), aes(x = model_name, y = estimated_proba_tox, fill = model_name)) +
      geom_hline(data = eval(parse(text = paste0("real_proba_tox_dose3$Scenario_", i))), aes(yintercept = real_proba_tox_dose3,  linetype = "Prior probability of toxicity"), color = "blue") +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "N.TITE-PK", "informed_tite_pk" = "Inf.TITE-PK")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(name = "", values = c(2, 1), guide = guide_legend(override.aes = list(color = c("blue", "red")))) +
      labs(title = paste("Scenario", i), y = "Probability of toxicity", x = "Dose-finding designs", fill = NULL, color = NULL, group = "to") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5, angle = 90), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_dose3_figure <- ggarrange(estimated_proba_tox_dose3_plots$Scenario_6 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_7 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_8 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_9 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_10 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_11 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_12 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_13 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_14 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_15 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_16 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_17 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_dose3_plots$Scenario_18 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_19 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_20 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"),
                                              estimated_proba_tox_dose3_plots$Scenario_21, 
                                              estimated_proba_tox_dose3_plots$Scenario_22 + rremove("ylab"), 
                                              estimated_proba_tox_dose3_plots$Scenario_23 + rremove("ylab"), 
                                              estimated_proba_tox_dose3_plots$Scenario_24 + rremove("ylab"), 
                                              estimated_proba_tox_dose3_plots$Scenario_25 + rremove("ylab"), 
                                              common.legend = TRUE, legend = "bottom", ncol = 5, nrow = 4, heights = c(1, 1, 1, 1.35))
#estimated_proba_tox_dose3_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for dose 3", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_dose3_figure.png"), plot = wrap_elements(estimated_proba_tox_dose3_figure), dpi = "retina", units = "in", width = 12, height = 8)



# SUMMARY PLOT BY SET OF SCENARIOS WHEN THE MTD IS LOCATED ON DOSE 1 USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF 
# TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD

plotted_scenarios_MTDdose1 <- c("A1", "B1", "C1", "D1", "E1")
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
set_name_list <- c("A", "B", "C", "D", "E")
estimated_proba_tox_MTDdose1_all_scenarios <- list()
estimated_proba_tox_MTDdose1_plots <- list()
real_proba_tox_MTDdose1 <- list()
for (i in plotted_scenarios_MTDdose1) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_MTDdose1[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_MTDdose1 = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  current_set_name <- set_name_list[set_name_cursor]
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_MTDdose1_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_MTDdose1_scenario$model_name <- factor(estimated_proba_tox_MTDdose1_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_MTDdose1_scenario$dose_num <- factor(estimated_proba_tox_MTDdose1_scenario$dose_num)
  estimated_proba_tox_MTDdose1_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_MTDdose1_scenario
  estimated_proba_tox_MTDdose1_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_MTDdose1_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_MTDdose1$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_MTDdose1, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", paste0(current_set_name, "1")), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "Dose-finding designs") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_MTDdose1_figure <- ggarrange(estimated_proba_tox_MTDdose1_plots$Scenario_A1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_MTDdose1_plots$Scenario_B1 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_MTDdose1_plots$Scenario_C1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                              estimated_proba_tox_MTDdose1_plots$Scenario_D1 + rremove("ylab"), 
                                              estimated_proba_tox_MTDdose1_plots$Scenario_E1,
                                              common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_MTDdose1_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_MTDdose1_figure.png"), plot = wrap_elements(estimated_proba_tox_MTDdose1_figure), dpi = "retina", units = "mm", width = 210, height = 230) #A4 Format: 297Heightx210Width


# SUMMARY PLOT BY SET OF SCENARIOS WHEN THE MTD IS LOCATED ON DOSE 2 USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF 
# TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD

plotted_scenarios_MTDdose2 <- c("A2", "B2", "C2", "D2", "E2")
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
set_name_list <- c("A", "B", "C", "D", "E")
estimated_proba_tox_MTDdose2_all_scenarios <- list()
estimated_proba_tox_MTDdose2_plots <- list()
real_proba_tox_MTDdose2 <- list()
for (i in plotted_scenarios_MTDdose2) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_MTDdose2[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_MTDdose2 = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  current_set_name <- set_name_list[set_name_cursor]
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_MTDdose2_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_MTDdose2_scenario$model_name <- factor(estimated_proba_tox_MTDdose2_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_MTDdose2_scenario$dose_num <- factor(estimated_proba_tox_MTDdose2_scenario$dose_num)
  estimated_proba_tox_MTDdose2_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_MTDdose2_scenario
  estimated_proba_tox_MTDdose2_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_MTDdose2_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_MTDdose2$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_MTDdose2, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", paste0(current_set_name, "2")), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "Dose-finding designs") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_MTDdose2_figure <- ggarrange(estimated_proba_tox_MTDdose2_plots$Scenario_A2 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose2_plots$Scenario_B2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose2_plots$Scenario_C2 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose2_plots$Scenario_D2 + rremove("ylab"), 
                                                 estimated_proba_tox_MTDdose2_plots$Scenario_E2,
                                                 common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_MTDdose2_figure <- annotate_figure(estimated_proba_tox_dose2_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_MTDdose2_figure.png"), plot = wrap_elements(estimated_proba_tox_MTDdose2_figure), dpi = "retina", units = "mm", width = 210, height = 230) #A4 Format: 297Heightx210Width


# SUMMARY PLOT BY SET OF SCENARIOS WHEN THE MTD IS LOCATED ON DOSE 3 USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF 
# TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD

plotted_scenarios_MTDdose3 <- c("A3", "B3", "C3", "D3", "E3")
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
set_name_list <- c("A", "B", "C", "D", "E", "F")
estimated_proba_tox_MTDdose3_all_scenarios <- list()
estimated_proba_tox_MTDdose3_plots <- list()
real_proba_tox_MTDdose3 <- list()
for (i in plotted_scenarios_MTDdose3) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_MTDdose3[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_MTDdose3 = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  current_set_name <- set_name_list[set_name_cursor]
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_MTDdose3_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_MTDdose3_scenario$model_name <- factor(estimated_proba_tox_MTDdose3_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_MTDdose3_scenario$dose_num <- factor(estimated_proba_tox_MTDdose3_scenario$dose_num)
  estimated_proba_tox_MTDdose3_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_MTDdose3_scenario
  estimated_proba_tox_MTDdose3_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_MTDdose3_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_MTDdose3$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_MTDdose3, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", paste0(current_set_name, "3")), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "Dose-finding designs") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_MTDdose3_figure <- ggarrange(estimated_proba_tox_MTDdose3_plots$Scenario_A3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose3_plots$Scenario_B3 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose3_plots$Scenario_C3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose3_plots$Scenario_D3 + rremove("ylab"), 
                                                 estimated_proba_tox_MTDdose3_plots$Scenario_E3,
                                                 common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_MTDdose3_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_MTDdose3_figure.png"), plot = wrap_elements(estimated_proba_tox_MTDdose3_figure), dpi = "retina", units = "mm", width = 210, height = 230) #A4 Format: 297Heightx210Width


# SUMMARY PLOT BY SET OF SCENARIOS WHEN THE MTD IS LOCATED ON DOSE 4 USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF 
# TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD

plotted_scenarios_MTDdose4 <- c("A4", "B4", "C4", "D4", "E4")
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
set_name_list <- c("A", "B", "C", "D", "E")
estimated_proba_tox_MTDdose4_all_scenarios <- list()
estimated_proba_tox_MTDdose4_plots <- list()
real_proba_tox_MTDdose4 <- list()
for (i in plotted_scenarios_MTDdose4) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_MTDdose4[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_MTDdose4 = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  current_set_name <- set_name_list[set_name_cursor]
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_MTDdose4_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_MTDdose4_scenario$model_name <- factor(estimated_proba_tox_MTDdose4_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_MTDdose4_scenario$dose_num <- factor(estimated_proba_tox_MTDdose4_scenario$dose_num)
  estimated_proba_tox_MTDdose4_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_MTDdose4_scenario
  estimated_proba_tox_MTDdose4_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_MTDdose4_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_MTDdose4$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_MTDdose4, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", paste0(current_set_name, "4")), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "Dose-finding designs") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_MTDdose4_figure <- ggarrange(estimated_proba_tox_MTDdose4_plots$Scenario_A4 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose4_plots$Scenario_B4 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose4_plots$Scenario_C4 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose4_plots$Scenario_D4 + rremove("ylab"), 
                                                 estimated_proba_tox_MTDdose4_plots$Scenario_E4,
                                                 common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_MTDdose4_figure <- annotate_figure(estimated_proba_tox_dose4_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_MTDdose4_figure.png"), plot = wrap_elements(estimated_proba_tox_MTDdose4_figure), dpi = "retina", units = "mm", width = 210, height = 230) #A4 Format: 297Heightx210Width


# SUMMARY PLOT BY SET OF SCENARIOS WHEN ALL DOSES ARE TOO TOXIC AND TRIAL SHOULD BE STOPPED EARLY USING BOXPLOTS
# FOR THE ESTIMATED PROBABILITIES OF TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD

plotted_scenarios_MTDdose0 <- c("A5", "B5", "C5", "D5", "E5")
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
set_name_list <- c("A", "B", "C", "D", "E")
estimated_proba_tox_MTDdose0_all_scenarios <- list()
estimated_proba_tox_MTDdose0_plots <- list()
real_proba_tox_MTDdose0 <- list()
for (i in plotted_scenarios_MTDdose0) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_MTDdose0[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_MTDdose0 = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  current_set_name <- set_name_list[set_name_cursor]
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_MTDdose0_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_MTDdose0_scenario$model_name <- factor(estimated_proba_tox_MTDdose0_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_MTDdose0_scenario$dose_num <- factor(estimated_proba_tox_MTDdose0_scenario$dose_num)
  estimated_proba_tox_MTDdose0_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_MTDdose0_scenario
  estimated_proba_tox_MTDdose0_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_MTDdose0_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_MTDdose0$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_MTDdose0, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", paste0(current_set_name, "5")), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "Dose-finding designs") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_MTDdose0_figure <- ggarrange(estimated_proba_tox_MTDdose0_plots$Scenario_A5 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose0_plots$Scenario_B5 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose0_plots$Scenario_C5 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_MTDdose0_plots$Scenario_D5 + rremove("ylab"), 
                                                 estimated_proba_tox_MTDdose0_plots$Scenario_E5,
                                                 common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_MTDdose0_figure <- annotate_figure(estimated_proba_tox_dose0_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_MTDdose0_figure.png"), plot = wrap_elements(estimated_proba_tox_MTDdose0_figure), dpi = "retina", units = "mm", width = 210, height = 230) #A4 Format: 297Heightx210Width


# SUMMARY PLOT WHEN THE MTD IS AT DOSE 1, 2, 3, 4 OR THERE IS NOT MTD AT ALL USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF 
# TOXICITY FOR ALL DOSES AND DOSE-FINDING METHODS BY PICKING ONE SCENARIO IN EACH SET OF SCENARIOS

plotted_scenarios_eachMTD <- c("A1", "B2", "C3", "D4", "E5")
MTD_id_scenarios <- c(1, 2, 3, 4, 0)
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
set_name_list <- c("A", "B", "C", "D", "E")
estimated_proba_tox_eachMTD_all_scenarios <- list()
estimated_proba_tox_eachMTD_plots <- list()
real_proba_tox_eachMTD <- list()
for (i in plotted_scenarios_eachMTD) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_eachMTD[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_eachMTD = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  current_set_name <- set_name_list[set_name_cursor]
  current_scenario_num <- ifelse(MTD_id_scenarios[set_name_cursor] == 0, 5, MTD_id_scenarios[set_name_cursor])
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_eachMTD_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_eachMTD_scenario$model_name <- factor(estimated_proba_tox_eachMTD_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_eachMTD_scenario$dose_num <- factor(estimated_proba_tox_eachMTD_scenario$dose_num)
  estimated_proba_tox_eachMTD_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_eachMTD_scenario
  estimated_proba_tox_eachMTD_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_eachMTD_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_eachMTD$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_eachMTD, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", paste0(current_set_name, current_scenario_num)), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_eachMTD_figure <- ggarrange(estimated_proba_tox_eachMTD_plots$Scenario_A1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_eachMTD_plots$Scenario_B2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_eachMTD_plots$Scenario_C3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                 estimated_proba_tox_eachMTD_plots$Scenario_D4 + rremove("xylab"), 
                                                 estimated_proba_tox_eachMTD_plots$Scenario_E5,
                                                 common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_eachMTD_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_eachMTD_figure.png"), plot = wrap_elements(estimated_proba_tox_eachMTD_figure), dpi = "retina", units = "mm", width = 210, height = 230) #230Heightx210Width


# SUMMARY PLOT FOR SET A USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD
plotted_scenarios_setA <- c("A1", "A2", "A3", "A4", "A5")
MTD_id_scenarios <- c(1, 2, 3, 4, 0)
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
estimated_proba_tox_setA_all_scenarios <- list()
estimated_proba_tox_setA_plots <- list()
real_proba_tox_setA <- list()
for (i in plotted_scenarios_setA) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_setA[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_setA = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_setA_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_setA_scenario$model_name <- factor(estimated_proba_tox_setA_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_setA_scenario$dose_num <- factor(estimated_proba_tox_setA_scenario$dose_num)
  estimated_proba_tox_setA_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_setA_scenario
  estimated_proba_tox_setA_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_setA_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_setA$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_setA, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", i), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_setA_figure <- ggarrange(estimated_proba_tox_setA_plots$Scenario_A1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                estimated_proba_tox_setA_plots$Scenario_A2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                                estimated_proba_tox_setA_plots$Scenario_A3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                                estimated_proba_tox_setA_plots$Scenario_A4 + rremove("xylab"), 
                                                estimated_proba_tox_setA_plots$Scenario_A5,
                                                common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_setA_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_setA_figure.png"), plot = wrap_elements(estimated_proba_tox_setA_figure), dpi = "retina", units = "mm", width = 210, height = 230) #230Heightx210Width


# SUMMARY PLOT FOR SET B USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD
plotted_scenarios_setB <- c("B1", "B2", "B3", "B4", "B5")
MTD_id_scenarios <- c(1, 2, 3, 4, 0)
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
estimated_proba_tox_setB_all_scenarios <- list()
estimated_proba_tox_setB_plots <- list()
real_proba_tox_setB <- list()
for (i in plotted_scenarios_setB) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_setB[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_setB = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_setB_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_setB_scenario$model_name <- factor(estimated_proba_tox_setB_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_setB_scenario$dose_num <- factor(estimated_proba_tox_setB_scenario$dose_num)
  estimated_proba_tox_setB_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_setB_scenario
  estimated_proba_tox_setB_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_setB_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_setB$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_setB, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", i), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_setB_figure <- ggarrange(estimated_proba_tox_setB_plots$Scenario_B1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setB_plots$Scenario_B2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setB_plots$Scenario_B3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setB_plots$Scenario_B4 + rremove("xylab"), 
                                             estimated_proba_tox_setB_plots$Scenario_B5,
                                             common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_setB_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_setB_figure.png"), plot = wrap_elements(estimated_proba_tox_setB_figure), dpi = "retina", units = "mm", width = 210, height = 230) #230Heightx210Width


# SUMMARY PLOT FOR SET C USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD
plotted_scenarios_setC <- c("C1", "C2", "C3", "C4", "C5")
MTD_id_scenarios <- c(1, 2, 3, 4, 0)
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
estimated_proba_tox_setC_all_scenarios <- list()
estimated_proba_tox_setC_plots <- list()
real_proba_tox_setC <- list()
for (i in plotted_scenarios_setC) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_setC[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_setC = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_setC_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_setC_scenario$model_name <- factor(estimated_proba_tox_setC_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_setC_scenario$dose_num <- factor(estimated_proba_tox_setC_scenario$dose_num)
  estimated_proba_tox_setC_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_setC_scenario
  estimated_proba_tox_setC_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_setC_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_setC$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_setC, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", i), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_setC_figure <- ggarrange(estimated_proba_tox_setC_plots$Scenario_C1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setC_plots$Scenario_C2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setC_plots$Scenario_C3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setC_plots$Scenario_C4 + rremove("xylab"), 
                                             estimated_proba_tox_setC_plots$Scenario_C5,
                                             common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_setC_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_setC_figure.png"), plot = wrap_elements(estimated_proba_tox_setC_figure), dpi = "retina", units = "mm", width = 210, height = 230) #230Heightx210Width


# SUMMARY PLOT FOR SET D USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD
plotted_scenarios_setD <- c("D1", "D2", "D3", "D4", "D5")
MTD_id_scenarios <- c(1, 2, 3, 4, 0)
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
estimated_proba_tox_setD_all_scenarios <- list()
estimated_proba_tox_setD_plots <- list()
real_proba_tox_setD <- list()
for (i in plotted_scenarios_setD) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_setD[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_setD = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_setD_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_setD_scenario$model_name <- factor(estimated_proba_tox_setD_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_setD_scenario$dose_num <- factor(estimated_proba_tox_setD_scenario$dose_num)
  estimated_proba_tox_setD_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_setD_scenario
  estimated_proba_tox_setD_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_setD_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_setD$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_setD, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", i), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_setD_figure <- ggarrange(estimated_proba_tox_setD_plots$Scenario_D1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setD_plots$Scenario_D2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setD_plots$Scenario_D3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setD_plots$Scenario_D4 + rremove("xylab"), 
                                             estimated_proba_tox_setD_plots$Scenario_D5,
                                             common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_setD_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_setD_figure.png"), plot = wrap_elements(estimated_proba_tox_setD_figure), dpi = "retina", units = "mm", width = 210, height = 230) #230Heightx210Width


# SUMMARY PLOT FOR SET E USING BOXPLOTS FOR THE ESTIMATED PROBABILITIES OF TOXICITY FOR ALL DOSES BY EACH DOSE-FINDING METHOD
plotted_scenarios_setE <- c("E1", "E2", "E3", "E4", "E5")
MTD_id_scenarios <- c(1, 2, 3, 4, 0)
res_model_name <- c("blrm", "pklogit", "ed_ewoc", "ed", "naive_tite_pk", "informed_tite_pk")
set_name_cursor <- 0
estimated_proba_tox_setE_all_scenarios <- list()
estimated_proba_tox_setE_plots <- list()
real_proba_tox_setE <- list()
for (i in plotted_scenarios_setE) {
  load(paste0(folder_path, "run_scen/scenario_", i, ".Rdata"))
  model_name <- c()
  dose_num <- c()
  estimated_proba_tox <- c()
  real_proba_tox_setE[[paste0("Scenario_", i, sep = "")]] <- data.frame(dose_num = c("dose1", "dose2", "dose3", "dose4"), real_proba_tox_setE = as.vector(eval(parse(text = paste0("res_blrm", "@prior_tox_real")))))
  set_name_cursor <- set_name_cursor + 1
  for (j in 1:length(res_model_name)) {
    model_name <- append(model_name, rep(res_model_name[j], length(doses)*n_trials))
    dose_num <- append(dose_num, rep(c("dose1", "dose2", "dose3", "dose4"), n_trials))
    estimated_proba_tox <- append(estimated_proba_tox,  as.vector(matrix(unlist(eval(parse(text = paste0("res_", res_model_name[j], "@p_estim")))), nrow = length(doses)*n_trials, ncol = 1, byrow = TRUE, dimnames = list(c(1:(length(doses)*n_trials)), "proba_tox"))))
  }
  estimated_proba_tox_setE_scenario <- data.frame(model_name = model_name, dose_num = dose_num, estimated_proba_tox = estimated_proba_tox)
  estimated_proba_tox_setE_scenario$model_name <- factor(estimated_proba_tox_setE_scenario$model_name, levels = res_model_name)
  estimated_proba_tox_setE_scenario$dose_num <- factor(estimated_proba_tox_setE_scenario$dose_num)
  estimated_proba_tox_setE_all_scenarios[[paste0("Scenario_", i, sep = "")]] <- estimated_proba_tox_setE_scenario
  estimated_proba_tox_setE_plots[[paste0("Scenario_", i, sep = "")]] <- local({
    ggplot() +
      geom_boxplot(data = eval(parse(text = paste0("estimated_proba_tox_setE_all_scenarios$Scenario_", i))), aes(x = dose_num, y = estimated_proba_tox, fill = model_name)) +
      geom_point(data = eval(parse(text = paste0("real_proba_tox_setE$Scenario_", i))), aes(x = dose_num, y = real_proba_tox_setE, group = 1, color = "Real probability of toxicity"), shape = 19, size = 3) +
      geom_hline(aes(yintercept = targeted_tox, linetype = "Target toxicity"), color = "red") +
      scale_x_discrete(labels = c("dose1" = "Dose 1", "dose2" = "Dose 2", "dose3" = "Dose 3", "dose4" = "Dose 4")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_brewer(palette="Dark2", labels = c("blrm" = "BLRM", "pklogit" = "PKLOGIT", "ed_ewoc" = "ED-EWOC", "ed" = "ED", "naive_tite_pk" = "Naive TITE-PK", "informed_tite_pk" = "Informed TITE-PK")) +
      scale_linetype_manual(values = 1, guide = guide_legend(override.aes = list(color = "red"))) +
      scale_color_manual(values = "blue", guide = guide_legend(override.aes = list(color = "blue"))) +
      labs(title = paste("Scenario", i), y = "Probability of toxicity", x = "Panel of doses", linetype = "", color = "", fill = "") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold"), axis.title.x = element_text(size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=8), axis.text.y = element_text(size=5))
  })
}
estimated_proba_tox_setE_figure <- ggarrange(estimated_proba_tox_setE_plots$Scenario_E1 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setE_plots$Scenario_E2 + rremove("xylab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setE_plots$Scenario_E3 + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
                                             estimated_proba_tox_setE_plots$Scenario_E4 + rremove("xylab"), 
                                             estimated_proba_tox_setE_plots$Scenario_E5,
                                             common.legend = TRUE, legend = "none", ncol = 2, nrow = 3, heights = c(1, 1, 1.1))
#estimated_proba_tox_setE_figure <- annotate_figure(estimated_proba_tox_dose3_figure, top = text_grob("Estimated probability of toxicity by all dose-finding method for all doses", color = "black", face = "bold", size = 14))
ggsave(paste0("plots/estimated_proba_tox_setE_figure.png"), plot = wrap_elements(estimated_proba_tox_setE_figure), dpi = "retina", units = "mm", width = 210, height = 230) #230Heightx210Width


# PLOT OF MULTIPLE DRAWS OF FIXED PK PARAMETERS PROVIDED TO THE PSEUDO-PK MODEL IN THE NAIVE VERSION
# OF TITE-PK DESIGN TO REPRESENT THE POPULATION CONCENTRATION CURVES USED BY THE DESIGN
# AND COMPARE THESE CURVES TO THE TRUE CONCENTRATION CURVE FROM THE SIMULATED PK DATA
library(tidyr)
set.seed(12345)
n_curves <- 50
time_sampling <- seq(0, 24, 0.1)
## T_e: half elimination rate constant
## Calculated from PK analysis (mean of the estimated T_e values)
Cl <- runif(n_curves, 0, 50)
V <- runif(n_curves, 50, 200)
k <- Cl / V #Cl_pop/V_pop
T_e <- log(2) / k
## k_eq: kinetic constant which govern delay between concentration in central compartment and effect site
## Calculated from PK analysis
k_eq <- runif(n_curves, 0, 10)
naive_TITEPK_pseudoPKparameters <- c("k" = k, "k_e" = k_e)
N_simu <- 1000
sim_PKdata <- simu_data(
  PK_parameters = PK_parameters,
  omega_IIV = omega_IIV,
  omega_alpha = omega_alpha,
  cv = cv,
  doses = doses,
  exposure_tox_threshold = exposure_tox_threshold,
  n_sampling_timepoints = n_sampling_timepoints,
  time_sampling = time_sampling,
  N = N_simu,
  n_trials = 1,
  exposure_metric = "AUC",
  AUC_method = "compartmental",
  compartmental_model = "1-cmt"
)
sim_PKdata_pop_curve <- lapply(seq_along(doses), function(x) colMeans(sim_PKdata[[1]]$simu_exact_conc[[x]]))
original_naive_TITEPK_pop_curves <- lapply(seq_along(doses), function(x) {
  return(matrix(sapply(1:n_curves, function(y) doses[x] * (k_eq[y]/(k_eq[y] - k[y])) * (exp(-k[y]*time_sampling) - exp(-k_eq[y]*time_sampling))), nrow = n_curves, ncol = length(time_sampling), byrow = TRUE))
})
corrected_naive_TITEPK_pop_curves <- lapply(seq_along(doses), function(x) {
  return(matrix(sapply(1:n_curves, function(y) (doses[x]/V[y]) * (k_eq[y]/(k_eq[y] - k[y])) * (exp(-k[y]*time_sampling) - exp(-k_eq[y]*time_sampling))), nrow = n_curves, ncol = length(time_sampling), byrow = TRUE))
})
prepare_ggplot_data <- function(doses, sim_PKdata_pop_curve, naive_TITEPK_pop_curves, time_sampling) {
  plot_data <- data.frame()
  
  for (i in seq_along(doses)) {
    df_pop <- data.frame(
      Time = time_sampling,
      Concentration = sim_PKdata_pop_curve[[i]],
      Dose = factor(doses[i]),
      Type = "Simulated PKdata"
    )
    
    df_naive <- as.data.frame(t(naive_TITEPK_pop_curves[[i]]))
    colnames(df_naive) <- paste0("Curve_", seq_len(ncol(df_naive)))
    df_naive$Time <- time_sampling
    df_naive <- df_naive %>%
      pivot_longer(cols = starts_with("Curve_"), names_to = "Curve", values_to = "Concentration") %>%
      mutate(Dose = factor(doses[i]), Type = "Naive TITE-PK")
    
    # Fusionner les données
    plot_data <- bind_rows(plot_data, df_pop, df_naive)
    plot_data$Type <- factor(plot_data$Type)
  }
  
  return(plot_data)
}
original_plot_data <- prepare_ggplot_data(doses, sim_PKdata_pop_curve, original_naive_TITEPK_pop_curves, time_sampling)
corrected_plot_data <- prepare_ggplot_data(doses, sim_PKdata_pop_curve, corrected_naive_TITEPK_pop_curves, time_sampling)
PKsim_pop_curve_vs_original_naive_TITEPK_pop_curves_dose1_figure <- ggplot(original_plot_data, aes(x = Time, y = Concentration, group = Type, color = Type)) +
  geom_line(data = filter(original_plot_data, Type == "Naive TITE-PK", Dose == doses[1]), aes(group = Curve), alpha = 0.5) +
  geom_line(data = filter(original_plot_data, Type == "Simulated PKdata", Dose == doses[1]), size = 1) +
  scale_color_manual(values = c("Naive TITE-PK" = "grey", "Simulated PKdata" = "red")) +
  theme_bw() +
  labs(title = "Baseline naive TITE-PK curves", x = "Time (hours)", y = "Concentration (mg/L)", color = "Type") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
PKsim_pop_curve_vs_original_naive_TITEPK_pop_curves_dose1_figure
PKsim_pop_curve_vs_corrected_naive_TITEPK_pop_curves_dose1_figure <- ggplot(corrected_plot_data, aes(x = Time, y = Concentration, group = Type, color = Type)) +
  geom_line(data = filter(corrected_plot_data, Type == "Naive TITE-PK", Dose == doses[1]), aes(group = Curve), alpha = 0.5) +
  geom_line(data = filter(corrected_plot_data, Type == "Simulated PKdata", Dose == doses[1]), size = 1) +
  scale_color_manual(values = c("Naive TITE-PK" = "grey", "Simulated PKdata" = "red")) +
  theme_bw() +
  labs(title = "Standardized by V", x = "Time (hours)", y = "Concentration (mg/L)", color = "Type") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
PKsim_pop_curve_vs_corrected_naive_TITEPK_pop_curves_dose1_figure
PKsim_pop_curve_vs_naive_TITEPK_pop_curves_dose1_final_figure <- ggarrange(PKsim_pop_curve_vs_original_naive_TITEPK_pop_curves_dose1_figure,
                                                                           PKsim_pop_curve_vs_corrected_naive_TITEPK_pop_curves_dose1_figure,
                                                                           common.legend = TRUE, legend = "none", ncol = 2, nrow = 1)
PKsim_pop_curve_vs_naive_TITEPK_pop_curves_dose1_final_figure
ggsave(paste0("plots/PKsim_pop_curve_vs_naive_TITEPK_pop_curves_dose1_figure.png"), plot = wrap_elements(PKsim_pop_curve_vs_naive_TITEPK_pop_curves_dose1_final_figure), dpi = "retina", units = "mm", width = 210, height = 230)





