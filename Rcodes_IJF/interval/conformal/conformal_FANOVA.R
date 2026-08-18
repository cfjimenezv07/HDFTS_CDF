# ==============================================================================
# Directory Setup & Source Dependencies
# ==============================================================================
source("load_packages.R")
source(here("auxiliary_source", "auxiliary_interval.R"))
source(here("auxiliary_source", "CDF_transformation.R"))

# --- Define Directories ---
dir.d      <- here("data", "updated_data")
dir.r      <- here("results", "Results_Interval", "Conformal")
dir.shiny2 <- here("results", "Shiny_App", "datasets_shiny_app", "IFE")

# Create directories if they do not exist
if (!dir.exists(dir.r))      dir.create(dir.r, recursive = TRUE)
if (!dir.exists(dir.shiny2)) dir.create(dir.shiny2, recursive = TRUE)


# ==============================================================================
# Data Loading and Preprocessing
# ==============================================================================

year         <- 1975:2023
n_year       <- length(year)
age          <- 0:110
n_age        <- length(age)
n_prefectures<- 47
n_populations<- 2

# Row & Column Partitions
part_list   <- lapply(1:n_prefectures, function(ik) (n_year * ik - (n_year - 1)):(n_year * ik))
part_list_c <- lapply(1:n_populations, function(ik) (n_age * ik - (n_age - 1)):(n_age * ik))

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

# Convert qx to dx death counts
female_prefecture_dx <- male_prefecture_dx <- total_prefecture_dx <- list()

for (iw in seq_along(state)) {
  female_dum <- male_dum <- total_dum <- matrix(NA, n_year, n_age)
  for (ij in 1:n_year) {
    pop_F <- pop_M <- pop_T <- 10^5
    for (ik in 1:n_age) {
      female_dum[ij, ik] <- female_prefecture_qx[[iw]][ij, ik] * pop_F
      pop_F              <- pop_F - female_dum[ij, ik]
      
      male_dum[ij, ik]   <- male_prefecture_qx[[iw]][ij, ik] * pop_M
      pop_M              <- pop_M - male_dum[ij, ik]
      
      total_dum[ij, ik]  <- total_prefecture_qx[[iw]][ij, ik] * pop_T
      pop_T              <- pop_T - total_dum[ij, ik]
    }
  }
  female_prefecture_dx[[iw]] <- t(female_dum)
  male_prefecture_dx[[iw]]   <- t(male_dum)
  total_prefecture_dx[[iw]]  <- t(total_dum)
}

# ==============================================================================
# Functional Mean ANOVA (FM-ANOVA) Transformation & Decomposition
# ==============================================================================

Transformed_female <- Transformed_male <- list()
for (i in 1:n_prefectures) {
  Transformed_female[[i]] <- t(cdf_transformation(t(female_prefecture_dx[[i]]) / 10^5, year))
  Transformed_male[[i]]   <- t(cdf_transformation(t(male_prefecture_dx[[i]]) / 10^5, year))
}

all_unconstrained_female <- t(list.cbind(Transformed_female))
all_unconstrained_male   <- t(list.cbind(Transformed_male))

new_age <- 0:109
nn_age  <- length(new_age)

# FM-ANOVA Decomposition
FANOVA_means    <- hdftsa::Two_way_mean(data_pop1 = t(all_unconstrained_male), data_pop2 = t(all_unconstrained_female), year, new_age, n_prefectures, n_populations)
Residuals_means <- hdftsa::Two_way_mean_residuals(data_pop1 = t(all_unconstrained_male), data_pop2 = t(all_unconstrained_female), year, new_age, n_prefectures, n_populations)

Res1_means <- Residuals_means$residuals1_mean
Res2_means <- Residuals_means$residuals2_mean

Fixed_part_means   <- Residuals_means$Fixed_comp_mean
Fixed_part_means_1 <- Fixed_part_means[, 1:nn_age]
Fixed_part_means_2 <- Fixed_part_means[, (nn_age + 1):(2 * nn_age)]

# Split fixed parts and residuals per prefecture unit
male_prefecture_res_means     <- lapply(1:length(part_list), function(k) Res1_means[part_list[[k]], ])
female_prefecture_res_means   <- lapply(1:length(part_list), function(k) Res2_means[part_list[[k]], ])

male_prefecture_fixed_means   <- lapply(1:length(part_list), function(k) Fixed_part_means_1[part_list[[k]], ])
female_prefecture_fixed_means <- lapply(1:length(part_list), function(k) Fixed_part_means_2[part_list[[k]], ])


# ==============================================================================
# Forecasting Functions
# ==============================================================================

select_k <- function(tau, eigenvalue) {
  k_max <- length(eigenvalue)
  k_all <- rep(0, k_max - 1)
  for (k in 1:(k_max - 1)) {
    k_all[k] <- (eigenvalue[k + 1] / eigenvalue[k]) * ifelse(eigenvalue[k] / eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k] / eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

Pref_forecast_curves <- function(fixed_com, Residuals_f, est_method = c("lrc", "cov"), fh = 30, B = 1000, 
                                 prediction_method = c("ARIMA", "VAR", "ets"), select_K = c("Fixed", "EVR"), K = 6) {
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
  
  var_total_variations   <- (sum(med_polish_resi_eigen$values[1:retain_component]) / sum(med_polish_resi_eigen$values)) * 100
  med_polish_resi_basis <- as.matrix(med_polish_resi_eigen$vectors[, 1:retain_component])
  med_polish_resi_score <- crossprod(med_polish_resi, med_polish_resi_basis)
  
  med_polish_resi_score_forecast <- matrix(NA, retain_component, fh)
  
  if (prediction_method == "ARIMA") {
    for (ip in 1:retain_component) {
      dum <- forecast_Arima(object = auto.arima(med_polish_resi_score[, ip]), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip, ] <- dum$mean
    }
  } else if (prediction_method == "VAR") {
    object <- med_polish_resi_score
    colnames(object) <- 1:dim(object)[2]
    lag <- VARselect(y = object, type = "const")$selection[1]
    model_VAR <- VAR(y = object, type = "const", ic = "AIC", p = lag)
    pred <- predict(model_VAR, n.ahead = fh)$fcst
    for (ip in 1:retain_component) {
      med_polish_resi_score_forecast[ip, ] <- pred[[ip]][, 1]
    }
  } else if (prediction_method == "ets") {
    for (ip in 1:retain_component) {
      dum <- forecast:::forecast.ets(object = ets(y = as.vector(med_polish_resi_score[, ip])), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip, ] <- as.vector(dum$mean)
    }
  }
  
  med_polish_resi_forecast  <- med_polish_resi_basis %*% med_polish_resi_score_forecast
  Fixed                     <- t(fixed_com)[, 1:fh]
  med_polish_curve_forecast <- med_polish_resi_forecast + Fixed
  
  return(list(med_polish_curve_forecast = med_polish_curve_forecast, 
              med_polish_resi_forecast  = med_polish_resi_forecast, 
              TV                        = var_total_variations))
}

fore_FANOVA_cdf <- function(Fixed, Residuals, ncomp_method, fh, fmethod = "ets", est_method) {
  data_cumsum_logit_fore <- Pref_forecast_curves(
    fixed_com = Fixed, Residuals_f = Residuals, est_method = est_method,
    fh = fh, B = 1000, prediction_method = fmethod, select_K = ncomp_method, K = 6
  )$med_polish_curve_forecast
  
  data_cumsum_logit_fore_add      <- c(invlogit(data_cumsum_logit_fore[, fh]), 1)
  data_cumsum_logit_fore_add_diff <- c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
  
  return(data_cumsum_logit_fore_add_diff * 10^5)
}


# ==============================================================================
# Main Split Conformal Prediction Function (FANOVA)
# ==============================================================================

interval_fore_FANOVA_cdf_conformal <- function(Fixed, Residuals, Holdout, method_ncomp, est_method, horizon, 
                                               fore_method = "CDF", uni_fore_method, level_sig) {
  fdata      <- t(Holdout)
  n_age      <- ncol(fdata)
  val_steps  <- 18 - horizon
  test_steps <- 17 - horizon
  
  fore_validation <- matrix(NA, n_age, val_steps)
  
  # Validation Phase
  if (fore_method == "CDF") {
    for (ij in 1:val_steps) {
      fore_validation[, ij] <- fore_FANOVA_cdf(
        Fixed = Fixed[1:(15 + ij), ], Residuals = Residuals[1:(15 + ij), ],
        ncomp_method = method_ncomp, fh = horizon, fmethod = uni_fore_method, est_method = est_method
      )
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:val_steps) {
      fore_validation[, ij] <- as.numeric(clr_fun(
        fdata = fdata[1:(15 + ij), ], ncomp_selection = method_ncomp, fh = horizon, fore_method = "ETS"
      )$fore_count)
    }
  } else {
    stop("Forecasting method must either be 'CDF' or 'CLR'.")
  }
  
  # Holdout validation residuals
  holdout_validation_dum <- t(fdata[(16 + horizon):33, , drop = FALSE])
  resi_mat               <- holdout_validation_dum - fore_validation
  quantile_resid         <- apply(resi_mat, 1, function(x) quantile(abs(x), probs = level_sig))
  
  # Testing Phase
  fore_val <- fore_val_lb <- fore_val_ub <- matrix(NA, n_age, test_steps)
  
  if (fore_method == "CDF") {
    for (ij in 1:test_steps) {
      fore_val[, ij]    <- fore_FANOVA_cdf(
        Fixed = Fixed[1:(32 + ij), ], Residuals = Residuals[1:(32 + ij), ],
        ncomp_method = method_ncomp, fh = horizon, fmethod = uni_fore_method, est_method = est_method
      )
      fore_val_lb[, ij] <- fore_val[, ij] - quantile_resid
      fore_val_ub[, ij] <- fore_val[, ij] + quantile_resid
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:test_steps) {
      fore_val[, ij]    <- as.numeric(clr_fun(
        fdata = fdata[1:(32 + ij), ], ncomp_selection = method_ncomp, fh = horizon, fore_method = "ETS"
      )$fore_count)
      fore_val_lb[, ij] <- fore_val[, ij] - quantile_resid
      fore_val_ub[, ij] <- fore_val[, ij] + quantile_resid
    }
  }
  
  # Interval Evaluation
  holdout_val_dum <- t(fdata[(33 + horizon):49, , drop = FALSE])
  int_err         <- interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, alpha = (1 - level_sig))
  
  return(list(int_err = int_err, dimn = ncol(resi_mat)))
}


# ==============================================================================
# Simulation Routine Wrapper
# ==============================================================================

run_fanova_conformal_experiment <- function(method_ncomp, level_sig) {
  err_F <- err_M <- array(NA, dim = c(16, 3, 47))
  
  for (ij in 1:47) {
    for (iw in 1:16) {
      # Female
      dum_F <- interval_fore_FANOVA_cdf_conformal(
        Fixed = female_prefecture_fixed_means[[ij]], Residuals = female_prefecture_res_means[[ij]],
        Holdout = female_prefecture_dx[[ij]], method_ncomp = method_ncomp, est_method = "cov",
        horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = level_sig
      )
      err_F[iw, , ij] <- dum_F$int_err
      
      # Male
      dum_M <- interval_fore_FANOVA_cdf_conformal(
        Fixed = male_prefecture_fixed_means[[ij]], Residuals = male_prefecture_res_means[[ij]],
        Holdout = male_prefecture_dx[[ij]], method_ncomp = method_ncomp, est_method = "cov",
        horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = level_sig
      )
      err_M[iw, , ij] <- dum_M$int_err
    }
    cat("Completed prefecture unit:", ij, "\n")
  }
  
  mean_F <- apply(err_F, c(1, 2), mean)
  mean_M <- apply(err_M, c(1, 2), mean)
  
  colnames(mean_F) <- colnames(mean_M) <- c("ECP", "CPD", "score")
  rownames(mean_F) <- rownames(mean_M) <- 1:16
  
  return(list(err_F = err_F, err_M = err_M, mean_F = mean_F, mean_M = mean_M))
}


# ==============================================================================
# Execute Experiments
# ==============================================================================

# Level of significance = 0.80
res_EVR_080 <- run_fanova_conformal_experiment(method_ncomp = "EVR",   level_sig = 0.80)
res_K6_080  <- run_fanova_conformal_experiment(method_ncomp = "Fixed", level_sig = 0.80)

# Level of significance = 0.95
res_EVR_095 <- run_fanova_conformal_experiment(method_ncomp = "EVR",   level_sig = 0.95)
res_K6_095  <- run_fanova_conformal_experiment(method_ncomp = "Fixed", level_sig = 0.95)


# ==============================================================================
# Save Outputs to dir.r (Prefixed with FANOVA_)
# ==============================================================================

save_fanova_results <- function(res, method_label, suffix) {
  # Arrays
  saveRDS(res$err_F, file.path(dir.r, paste0("FANOVA_int_fore_subnational_err_F_", method_label, "_conformal", suffix, ".rds")))
  saveRDS(res$err_M, file.path(dir.r, paste0("FANOVA_int_fore_subnational_err_M_", method_label, "_conformal", suffix, ".rds")))
  
  # Means
  saveRDS(res$mean_F, file.path(dir.r, paste0("FANOVA_int_fore_subnational_err_F_", method_label, "_conformal", suffix, "_mean.rds")))
  saveRDS(res$mean_M, file.path(dir.r, paste0("FANOVA_int_fore_subnational_err_M_", method_label, "_conformal", suffix, "_mean.rds")))
}

# --- level_sig = 0.80 ---
save_fanova_results(res_EVR_080, "EVR", "")
save_fanova_results(res_K6_080,  "K6",  "")

# --- level_sig = 0.95 ---
save_fanova_results(res_EVR_095, "EVR", "_alpha_0.95")
save_fanova_results(res_K6_095,  "K6",  "_alpha_0.95")


# ==============================================================================
# Export to Shiny Directory (Prefixed with FANOVA_)
# ==============================================================================

arrays_list_shiny <- list(
  FANOVA_female_EVR_80 = res_EVR_080$err_F,
  FANOVA_male_EVR_80   = res_EVR_080$err_M,
  FANOVA_female_K_80   = res_K6_080$err_F,
  FANOVA_male_K_80     = res_K6_080$err_M,
  
  FANOVA_female_EVR_95 = res_EVR_095$err_F,
  FANOVA_male_EVR_95   = res_EVR_095$err_M,
  FANOVA_female_K_95   = res_K6_095$err_F,
  FANOVA_male_K_95     = res_K6_095$err_M
)

for (name in names(arrays_list_shiny)) {
  file_name <- paste0(name, ".rds")
  saveRDS(arrays_list_shiny[[name]], file = file.path(dir.shiny2, file_name))
}

