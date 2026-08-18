# ==============================================================================
# 0. LIBRARIES & DIRECTORIES
# ==============================================================================
source("load_packages.R")

dir_source  <- here("auxiliary_source")
dir_results <- here("results", "Results_Interval", "sd")
dir.shiny2  <- here("results", "Shiny_App", "datasets_shiny_app", "IFE")
dir.d       <- here("data", "updated_data")

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)
if (!dir.exists(dir.shiny2))  dir.create(dir.shiny2, recursive = TRUE)

# ==============================================================================
# 1. DATASET SETUP & LOADING
# ==============================================================================
year <- 1975:2023
n_year <- length(year)
age <- 0:110
n_age <- length(age)
n_prefectures <- 47
n_populations <- 2

# Row partition
part_list <- list()
for (ik in 1:n_prefectures) {
  part_list[[ik]] <- (n_year * ik - (n_year - 1)):(n_year * ik)
}

# Column partition
part_list_c <- list()
for (ik in 1:n_populations) {
  part_list_c[[ik]] <- (n_age * ik - (n_age - 1)):(n_age * ik)
}

state <- c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
           "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
           "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
           "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
           "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
           "Tokushima", "Kagawa", "Ehime", "Kochi",
           "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

female_prefecture_qx <- male_prefecture_qx <- total_prefecture_qx <- list()

for (ik in seq_along(state)) {
  female_data <- read.table(paste0(dir.d, "female_prefecture_", ik, ".txt"), header = TRUE, skip = 2)
  female_data <- female_data[female_data$Year >= 1975, ]
  female_prefecture_qx[[ik]] <- t(matrix(female_data$qx, n_age, n_year))
  
  male_data <- read.table(paste0(dir.d, "male_prefecture_", ik, ".txt"), header = TRUE, skip = 2)
  male_data <- male_data[male_data$Year >= 1975, ]
  male_prefecture_qx[[ik]] <- t(matrix(male_data$qx, n_age, n_year))
  
  total_data <- read.table(paste0(dir.d, "total_prefecture_", ik, ".txt"), header = TRUE, skip = 2)
  total_data <- total_data[total_data$Year >= 1975, ]
  total_prefecture_qx[[ik]] <- t(matrix(total_data$qx, n_age, n_year))
}

female_prefecture_dx <- male_prefecture_dx <- total_prefecture_dx <- list()
for (iw in 1:length(state)) {
  female_prefecture_dum <- male_prefecture_dum <- total_prefecture_dum <- matrix(NA, n_year, n_age)
  for (ij in 1:n_year) {
    start_pop_female <- start_pop_male <- start_pop_total <- 10^5
    for (ik in 1:n_age) {
      female_prefecture_dum[ij, ik] <- (female_prefecture_qx[[iw]])[ij, ik] * start_pop_female
      start_pop_female <- start_pop_female - female_prefecture_dum[ij, ik]
      
      male_prefecture_dum[ij, ik] <- (male_prefecture_qx[[iw]])[ij, ik] * start_pop_male
      start_pop_male <- start_pop_male - male_prefecture_dum[ij, ik]
      
      total_prefecture_dum[ij, ik] <- (total_prefecture_qx[[iw]])[ij, ik] * start_pop_total
      start_pop_total <- start_pop_total - total_prefecture_dum[ij, ik]
    }
  }
  female_prefecture_dx[[iw]] <- t(female_prefecture_dum)
  male_prefecture_dx[[iw]]   <- t(male_prefecture_dum)
  total_prefecture_dx[[iw]]  <- t(total_prefecture_dum)
}

All_Japan_female_qx <- female_prefecture_dx
All_Japan_male_qx   <- male_prefecture_dx
n_states <- length(state)

# ==============================================================================
# 2. CDF TRANSFORMATION & FANOVA DECOMPOSITION
# ==============================================================================
source(file.path(dir_source, "CDF_transformation.R"))

Transformed_female <- list()
Transformed_male <- list()
for (i in 1:n_states) {
  Transformed_female[[i]] <- t(cdf_transformation(t(All_Japan_female_qx[[i]]) / 10^5, year))
  Transformed_male[[i]]   <- t(cdf_transformation(t(All_Japan_male_qx[[i]]) / 10^5, year))
}

all_unconstrained_female <- t(list.cbind(Transformed_female))
all_unconstrained_male   <- t(list.cbind(Transformed_male))

new_age <- 0:109
FANOVA_means <- hdftsa::Two_way_mean(data_pop1 = t(all_unconstrained_male), data_pop2 = t(all_unconstrained_female), year, new_age, n_states, n_populations)

Residuals_means <- hdftsa::Two_way_mean_residuals(data_pop1 = t(all_unconstrained_male), data_pop2 = t(all_unconstrained_female), year, new_age, n_states, n_populations)

Res1_means <- Residuals_means$residuals1_mean
Res2_means <- Residuals_means$residuals2_mean

Fixed_part_means <- Residuals_means$Fixed_comp_mean
nn_age <- length(new_age)
Fixed_part_means_1 <- Fixed_part_means[, 1:nn_age]
Fixed_part_means_2 <- Fixed_part_means[, (nn_age + 1):(2 * nn_age)]

# Split per prefecture
male_prefecture_res_means   <- lapply(1:length(part_list), function(k) { Res1_means[part_list[[k]], ] })
female_prefecture_res_means <- lapply(1:length(part_list), function(k) { Res2_means[part_list[[k]], ] })

male_prefecture_fixed_means   <- lapply(1:length(part_list), function(k) { Fixed_part_means_1[part_list[[k]], ] })
female_prefecture_fixed_means <- lapply(1:length(part_list), function(k) { Fixed_part_means_2[part_list[[k]], ] })

# ==============================================================================
# 3. AUXILIARY FUNCTIONS FOR FORECASTING & INTERVAL TUNING
# ==============================================================================
tune_para_find_function <- function(tune_para, resi_mat, sd_val_input, alpha_level, PI_type) {
  n_age <- nrow(resi_mat)
  if (PI_type == "pointwise") {
    ind <- matrix(NA, n_age, ncol(resi_mat))
    for (iw in 1:ncol(resi_mat)) {
      ind[, iw] <- ifelse(between(resi_mat[, iw], -tune_para * sd_val_input, tune_para * sd_val_input), 1, 0)
    }
    ecp <- sum(ind) / (n_age * ncol(resi_mat))
  } else if (PI_type == "uniform") {
    ind <- vector("numeric", ncol(resi_mat))
    for (iw in 1:ncol(resi_mat)) {
      ind[iw] <- ifelse(all(between(resi_mat[, iw], -tune_para * sd_val_input, tune_para * sd_val_input)), 1, 0)
    }
    ecp <- sum(ind) / ncol(resi_mat)
  }
  return(abs(ecp - alpha_level))
}

interval_score <- function(holdout, lb, ub, alpha) {
  lb_ind <- ifelse(holdout < lb, 1, 0)
  ub_ind <- ifelse(holdout > ub, 1, 0)
  score  <- (ub - lb) + 2 / alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover  <- 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1))) / length(holdout)
  cpd    <- abs(cover - (1 - alpha))
  return(c(cover, cpd, mean(score)))
}

select_k <- function(tau, eigenvalue) {
  k_max <- length(eigenvalue)
  k_all <- rep(0, k_max - 1)
  for (k in 1:(k_max - 1)) {
    k_all[k] <- (eigenvalue[k + 1] / eigenvalue[k]) * ifelse(eigenvalue[k] / eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k] / eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

Pref_forecast_curves <- function(fixed_com, Residuals_f, est_method = c("lrc", "cov"), fh = 30, B = 1000, prediction_method = c("ARIMA", "VAR", "ets"), select_K = c("Fixed", "EVR"), K = 6) {
  med_polish_resi <- t(Residuals_f)
  if (est_method == "lrc") {
    med_polish_resi_lrc <- long_run_covariance_estimation(med_polish_resi)
  } else if (est_method == "cov") {
    med_polish_resi_lrc <- cov(t(med_polish_resi))
  }
  med_polish_resi_eigen <- eigen(med_polish_resi_lrc)
  
  if (select_K == "Fixed") {
    retain_component <- K
  } else if (select_K == "EVR") {
    retain_component <- select_k(tau = 10^-2, eigenvalue = med_polish_resi_eigen$values)
  }
  
  var_total_variations <- (sum(med_polish_resi_eigen$values[1:retain_component]) / sum(med_polish_resi_eigen$values)) * 100
  med_polish_resi_basis <- as.matrix(med_polish_resi_eigen$vectors[, 1:retain_component])
  med_polish_resi_score <- crossprod(med_polish_resi, med_polish_resi_basis)
  
  med_polish_resi_score_forecast <- matrix(NA, retain_component, fh)
  if (prediction_method == "ARIMA") {
    for (ip in 1:retain_component) {
      dum <- forecast_Arima(object = auto.arima(med_polish_resi_score[, ip]), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip, ] <- dum$mean
    }
  } else if (prediction_method == "ets") {
    for (ip in 1:retain_component) {
      dum <- forecast:::forecast.ets(object = ets(y = as.vector(med_polish_resi_score[, ip])), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip, ] <- as.vector(dum$mean)
    }
  }
  
  med_polish_resi_forecast <- med_polish_resi_basis %*% med_polish_resi_score_forecast
  Fixed <- t(fixed_com)[, 1:fh]
  med_polish_curve_forecast <- med_polish_resi_forecast + Fixed
  
  return(list(med_polish_curve_forecast = med_polish_curve_forecast, med_polish_resi_forecast = med_polish_resi_forecast, TV = var_total_variations))
}

fore_FANOVA_cdf <- function(Fixed, Residuals, ncomp_method, fh, fmethod, est_method) {
  data_cumsum_logit_fore <- Pref_forecast_curves(fixed_com = Fixed, Residuals_f = Residuals, est_method = est_method, fh = fh, B = 1000, prediction_method = fmethod, select_K = ncomp_method, K = 6)$med_polish_curve_forecast
  data_cumsum_logit_fore_add <- c(invlogit(data_cumsum_logit_fore[, fh]), 1)
  data_cumsum_logit_fore_add_diff <- c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
  return(data_cumsum_logit_fore_add_diff * 10^5)
}

interval_fore_FANOVA_cdf <- function(Fixed, Residuals, Holdout, method_ncomp, est_method, horizon, fore_method = "CDF", uni_fore_method, level_sig) {
  fdata <- t(Holdout)
  fore_validation <- matrix(NA, ncol(fdata), (18 - horizon))
  
  for (ij in 1:(18 - horizon)) {
    fore_validation[, ij] <- fore_FANOVA_cdf(Fixed = Fixed[1:(15 + ij), ], Residuals = Residuals[1:(15 + ij), ], ncomp_method = method_ncomp, fh = horizon, fmethod = uni_fore_method, est_method = est_method)
  }
  
  holdout_validation_dum <- t(matrix(fdata[(16 + horizon):33, ], length((16 + horizon):33), ncol(fdata)))
  resi_mat <- holdout_validation_dum - fore_validation
  sd_val_input <- apply(resi_mat, 1, sd)
  
  tune_para_find_val_1 <- optimise(f = tune_para_find_function, interval = c(0, 10), resi_mat = resi_mat, sd_val_input = sd_val_input, alpha_level = level_sig, PI_type = "pointwise")
  tune_para_find_val_2 <- optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B", resi_mat = resi_mat, sd_val_input = sd_val_input, alpha_level = level_sig, PI_type = "pointwise")
  tune_para_find_val_3 <- optim(par = 1, fn = tune_para_find_function, method = "Nelder-Mead", resi_mat = resi_mat, sd_val_input = sd_val_input, alpha_level = level_sig, PI_type = "pointwise")
  
  obj_val <- c(tune_para_find_val_1$objective, tune_para_find_val_2$value, tune_para_find_val_3$value)
  tune_para_find <- c(tune_para_find_val_1$minimum, tune_para_find_val_2$par, tune_para_find_val_3$par)[which.min(obj_val)]
  
  fore_val <- fore_val_lb <- fore_val_ub <- matrix(NA, ncol(fdata), (17 - horizon))
  for (ij in 1:(17 - horizon)) {
    fore_val[, ij] <- fore_FANOVA_cdf(Fixed = Fixed[1:(32 + ij), ], Residuals = Residuals[1:(32 + ij), ], ncomp_method = method_ncomp, fh = horizon, fmethod = uni_fore_method, est_method = est_method)
    fore_val_lb[, ij] <- fore_val[, ij] - tune_para_find * sd_val_input
    fore_val_ub[, ij] <- fore_val[, ij] + tune_para_find * sd_val_input
  }
  
  holdout_val_dum <- t(matrix(fdata[(33 + horizon):49, ], length((33 + horizon):49), ncol(fdata)))
  int_err <- interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, alpha = 1 - level_sig)
  
  return(list(int_err = int_err, tune_para_find = tune_para_find, tune_para_find_obj = min(obj_val)))
}

# ==============================================================================
# 4. LEVEL_SIG = 0.80
# ==============================================================================
# --- EVR ---
int_fore_subnational_err_F_EVR_ETS <- int_fore_subnational_err_M_EVR_ETS <- array(NA, dim = c(47, 16, 3), dimnames = list(state, 1:16, c("ECP", "CPD", "score")))
int_fore_subnational_err_F_EVR_ETS_tune_para <- int_fore_subnational_err_F_EVR_ETS_tune_para_obj <- 
  int_fore_subnational_err_M_EVR_ETS_tune_para <- int_fore_subnational_err_M_EVR_ETS_tune_para_obj <- matrix(NA, 47, 16)

for (ij in 1:47) {
  for (iw in 1:16) {
    ## (F)
    dum <- interval_fore_FANOVA_cdf(Fixed = female_prefecture_fixed_means[[ij]], Residuals = female_prefecture_res_means[[ij]],
                                    Holdout = female_prefecture_dx[[ij]], method_ncomp = "EVR", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.80)
    int_fore_subnational_err_F_EVR_ETS[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_F_EVR_ETS_tune_para[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_F_EVR_ETS_tune_para_obj[ij, iw] <- dum$tune_para_find_obj
    rm(dum)
    
    ## (M)
    dum <- interval_fore_FANOVA_cdf(Fixed = male_prefecture_fixed_means[[ij]], Residuals = male_prefecture_res_means[[ij]],
                                    Holdout = male_prefecture_dx[[ij]], method_ncomp = "EVR", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.80)
    int_fore_subnational_err_M_EVR_ETS[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_M_EVR_ETS_tune_para[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_M_EVR_ETS_tune_para_obj[ij, iw] <- dum$tune_para_find_obj
    rm(dum); rm(iw)
  }
  print(ij); rm(ij)
}

horizon_specific_int_fore_subnational_err_F_EVR_ETS <- apply(int_fore_subnational_err_F_EVR_ETS, c(2, 3), mean)
horizon_specific_int_fore_subnational_err_M_EVR_ETS <- apply(int_fore_subnational_err_M_EVR_ETS, c(2, 3), mean)

# --- K = 6 ---
int_fore_subnational_err_F_ETS_K6 <- int_fore_subnational_err_M_ETS_K6 <- array(NA, dim = c(47, 16, 3), dimnames = list(state, 1:16, c("ECP", "CPD", "score")))
int_fore_subnational_err_F_ETS_tune_para_K6 <- int_fore_subnational_err_F_ETS_tune_para_obj_K6 <- 
  int_fore_subnational_err_M_ETS_tune_para_K6 <- int_fore_subnational_err_M_ETS_tune_para_obj_K6 <- matrix(NA, 47, 16)

for (ij in 1:47) {
  for (iw in 1:16) {
    ## (F)
    dum <- interval_fore_FANOVA_cdf(Fixed = female_prefecture_fixed_means[[ij]], Residuals = female_prefecture_res_means[[ij]],
                                    Holdout = female_prefecture_dx[[ij]], method_ncomp = "Fixed", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.80)
    int_fore_subnational_err_F_ETS_K6[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_F_ETS_tune_para_K6[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_F_ETS_tune_para_obj_K6[ij, iw] <- dum$tune_para_find_obj
    rm(dum)
    
    ## (M)
    dum <- interval_fore_FANOVA_cdf(Fixed = male_prefecture_fixed_means[[ij]], Residuals = male_prefecture_res_means[[ij]],
                                    Holdout = male_prefecture_dx[[ij]], method_ncomp = "Fixed", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.80)
    int_fore_subnational_err_M_ETS_K6[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_M_ETS_tune_para_K6[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_M_ETS_tune_para_obj_K6[ij, iw] <- dum$tune_para_find_obj
    rm(dum); rm(iw)
  }
  print(ij); rm(ij)
}

horizon_specific_int_fore_subnational_err_F_ETS_K6 <- apply(int_fore_subnational_err_F_ETS_K6, c(2, 3), mean)
horizon_specific_int_fore_subnational_err_M_ETS_K6 <- apply(int_fore_subnational_err_M_ETS_K6, c(2, 3), mean)

# ==============================================================================
# 5. LEVEL_SIG = 0.95
# ==============================================================================
# --- EVR ---
int_fore_subnational_err_F_EVR_ETS_alpha_0.95 <- int_fore_subnational_err_M_EVR_ETS_alpha_0.95 <- array(NA, dim = c(47, 16, 3), dimnames = list(state, 1:16, c("ECP", "CPD", "score")))
int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95 <- int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95 <- 
  int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95 <- int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95 <- matrix(NA, 47, 16)

for (ij in 1:47) {
  for (iw in 1:16) {
    ## (F)
    dum <- interval_fore_FANOVA_cdf(Fixed = female_prefecture_fixed_means[[ij]], Residuals = female_prefecture_res_means[[ij]],
                                    Holdout = female_prefecture_dx[[ij]], method_ncomp = "EVR", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.95)
    int_fore_subnational_err_F_EVR_ETS_alpha_0.95[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95[ij, iw] <- dum$tune_para_find_obj
    rm(dum)
    
    ## (M) - Corrected inputs to male_prefecture_*
    dum <- interval_fore_FANOVA_cdf(Fixed = male_prefecture_fixed_means[[ij]], Residuals = male_prefecture_res_means[[ij]],
                                    Holdout = male_prefecture_dx[[ij]], method_ncomp = "EVR", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.95)
    int_fore_subnational_err_M_EVR_ETS_alpha_0.95[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95[ij, iw] <- dum$tune_para_find_obj
    rm(dum); rm(iw)
  }
  print(ij); rm(ij)
}

horizon_specific_int_fore_subnational_err_F_EVR_ETS_alpha_0.95 <- apply(int_fore_subnational_err_F_EVR_ETS_alpha_0.95, c(2, 3), mean)
horizon_specific_int_fore_subnational_err_M_EVR_ETS_alpha_0.95 <- apply(int_fore_subnational_err_M_EVR_ETS_alpha_0.95, c(2, 3), mean)

# --- K = 6 ---
int_fore_subnational_err_F_ETS_K6_alpha_0.95 <- int_fore_subnational_err_M_ETS_K6_alpha_0.95 <- array(NA, dim = c(47, 16, 3), dimnames = list(state, 1:16, c("ECP", "CPD", "score")))
int_fore_subnational_err_F_ETS_tune_para_K6_alpha_0.95 <- int_fore_subnational_err_F_ETS_tune_para_obj_K6_alpha_0.95 <- 
  int_fore_subnational_err_M_ETS_tune_para_K6_alpha_0.95 <- int_fore_subnational_err_M_ETS_tune_para_obj_K6_alpha_0.95 <- matrix(NA, 47, 16)

for (ij in 1:47) {
  for (iw in 1:16) {
    ## (F)
    dum <- interval_fore_FANOVA_cdf(Fixed = female_prefecture_fixed_means[[ij]], Residuals = female_prefecture_res_means[[ij]],
                                    Holdout = female_prefecture_dx[[ij]], method_ncomp = "Fixed", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.95)
    int_fore_subnational_err_F_ETS_K6_alpha_0.95[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_F_ETS_tune_para_K6_alpha_0.95[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_F_ETS_tune_para_obj_K6_alpha_0.95[ij, iw] <- dum$tune_para_find_obj
    rm(dum)
    
    ## (M)
    dum <- interval_fore_FANOVA_cdf(Fixed = male_prefecture_fixed_means[[ij]], Residuals = male_prefecture_res_means[[ij]],
                                    Holdout = male_prefecture_dx[[ij]], method_ncomp = "Fixed", est_method = "cov",
                                    horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = 0.95)
    int_fore_subnational_err_M_ETS_K6_alpha_0.95[ij, iw, ] <- dum$int_err
    int_fore_subnational_err_M_ETS_tune_para_K6_alpha_0.95[ij, iw] <- dum$tune_para_find
    int_fore_subnational_err_M_ETS_tune_para_obj_K6_alpha_0.95[ij, iw] <- dum$tune_para_find_obj
    rm(dum); rm(iw)
  }
  print(ij); rm(ij)
}

horizon_specific_int_fore_subnational_err_F_ETS_K6_alpha_0.95 <- apply(int_fore_subnational_err_F_ETS_K6_alpha_0.95, c(2, 3), mean)
horizon_specific_int_fore_subnational_err_M_ETS_K6_alpha_0.95 <- apply(int_fore_subnational_err_M_ETS_K6_alpha_0.95, c(2, 3), mean)

# ==============================================================================
# 6. SAVE ALL OBJECTS TO RESULTS DIRECTORY
# ==============================================================================
# --- 80% Level ---
saveRDS(int_fore_subnational_err_F_EVR_ETS, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS_tune_para.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS_tune_para.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para_obj, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS_tune_para_obj.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para_obj, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS_tune_para_obj.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_F_EVR_ETS, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_F_EVR_ETS.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_M_EVR_ETS, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_M_EVR_ETS.rds"))

saveRDS(int_fore_subnational_err_F_ETS_K6, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_ETS_K6.rds"))
saveRDS(int_fore_subnational_err_M_ETS_K6, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_ETS_K6.rds"))
saveRDS(int_fore_subnational_err_F_ETS_tune_para_K6, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_ETS_tune_para_K6.rds"))
saveRDS(int_fore_subnational_err_M_ETS_tune_para_K6, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_ETS_tune_para_K6.rds"))
saveRDS(int_fore_subnational_err_F_ETS_tune_para_obj_K6, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_ETS_tune_para_obj_K6.rds"))
saveRDS(int_fore_subnational_err_M_ETS_tune_para_obj_K6, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_ETS_tune_para_obj_K6.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_F_ETS_K6, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_F_ETS_K6.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_M_ETS_K6, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_M_ETS_K6.rds"))

# --- 95% Level ---
saveRDS(int_fore_subnational_err_F_EVR_ETS_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_F_EVR_ETS_alpha_0.95, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_F_EVR_ETS_alpha_0.95.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_M_EVR_ETS_alpha_0.95, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_M_EVR_ETS_alpha_0.95.rds"))

saveRDS(int_fore_subnational_err_F_ETS_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_ETS_K6_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_ETS_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_ETS_K6_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_ETS_tune_para_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_ETS_tune_para_K6_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_ETS_tune_para_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_ETS_tune_para_K6_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_ETS_tune_para_obj_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_F_ETS_tune_para_obj_K6_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_ETS_tune_para_obj_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_int_fore_subnational_err_M_ETS_tune_para_obj_K6_alpha_0.95.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_F_ETS_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_F_ETS_K6_alpha_0.95.rds"))
saveRDS(horizon_specific_int_fore_subnational_err_M_ETS_K6_alpha_0.95, file = file.path(dir_results, "FANOVA_horizon_specific_int_fore_subnational_err_M_ETS_K6_alpha_0.95.rds"))

# ==============================================================================
# 7. EXPORT TO SHINY APP DIRECTORY
# ==============================================================================
arrays_list <- list(
  F_EVR_80 = int_fore_subnational_err_F_EVR_ETS,
  M_EVR_80 = int_fore_subnational_err_M_EVR_ETS,
  F_K6_80  = int_fore_subnational_err_F_ETS_K6,
  M_K6_80  = int_fore_subnational_err_M_ETS_K6,
  F_EVR_95 = int_fore_subnational_err_F_EVR_ETS_alpha_0.95,
  M_EVR_95 = int_fore_subnational_err_M_EVR_ETS_alpha_0.95,
  F_K6_95  = int_fore_subnational_err_F_ETS_K6_alpha_0.95,
  M_K6_95  = int_fore_subnational_err_M_ETS_K6_alpha_0.95
)

for (name in names(arrays_list)) {
  arr    <- arrays_list[[name]]
  gender <- ifelse(grepl("^F", name), "female", "male")
  comp   <- ifelse(grepl("EVR", name), "EVR", "K")
  level  <- ifelse(grepl("80", name), "80", "95")
  
  file_name <- paste0("FANOVA_", gender, "_", comp, "_", level, ".rds")
  saveRDS(arr, file = file.path(dir.shiny2, file_name))
}
