# Load helper sources
source("load_packages.R")
source(here("auxiliary_source", "auxiliary_point.R"))
source(here("auxiliary_source", "CDF_transformation.R"))
source(here("auxiliary_source", "forecast_errors.R"))

# Centralized base paths using project-relative paths
dir_data    <- here("data")
dir_results <- here("results", "Results_Point")
dir_shiny   <- here("results", "Shiny_App")

dirs_to_create <- c(dir_data, dir_results, dir_shiny)
for (d in dirs_to_create) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# Dataset parameters
year         <- 1975:2023
n_year       <- length(year)
age          <- 0:110
n_age        <- length(age)
n_prefectures<- 47
n_populations<- 2

state <- c(
  "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
  "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa",
  "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", 
  "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", 
  "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
  "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", 
  "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa"
)
n_states <- length(state)

# Load processed mortality data
female_prefecture_dx <- readRDS(file.path(dir_data, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(file.path(dir_data, "male_prefecture_dx.rds"))

# CDF transformation
Transformed_female <- list()
Transformed_male   <- list()

for (i in 1:n_states) {
  Transformed_female[[i]] <- t(cdf_transformation((female_prefecture_dx[[i]]) / 10^5, year))
  Transformed_male[[i]]   <- t(cdf_transformation((male_prefecture_dx[[i]]) / 10^5, year))
}

all_unconstrained_female <- t(list.cbind(Transformed_female))
all_unconstrained_male   <- t(list.cbind(Transformed_male))

# FM-ANOVA Decomposition
new_age  <- 0:109
nn_age   <- length(new_age)

FANOVA_means <- hdftsa::Two_way_mean(
  data_pop1     = t(all_unconstrained_male),
  data_pop2     = t(all_unconstrained_female),
  year          = year,
  age           = new_age,
  n_prefectures = n_states,
  n_populations = n_populations
)

Residuals_means <- hdftsa::Two_way_mean_residuals(
  data_pop1     = t(all_unconstrained_male),
  data_pop2     = t(all_unconstrained_female),
  year          = year,
  age           = new_age,
  n_prefectures = n_states,
  n_populations = n_populations
)

Res1_means     <- Residuals_means$residuals1_mean
Res2_means     <- Residuals_means$residuals2_mean
Residuals_mean <- cbind(Res1_means, Res2_means)

Fixed_part_means <- Residuals_means$Fixed_comp_mean

# Component selection helper
select_k <- function(tau, eigenvalue) {
  k_max <- length(eigenvalue)
  k_all <- rep(0, k_max - 1)
  for (k in 1:(k_max - 1)) {
    k_all[k] <- (eigenvalue[k + 1] / eigenvalue[k]) * ifelse(eigenvalue[k] / eigenvalue[1] > tau, 1, 0) +
      ifelse(eigenvalue[k] / eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

# Curve Forecasting function
Pref_forecast_curves <- function(fixed_com, Residuals_f, est_method = c("lrc", "cov"),
                                 fh = 30, B = 1000, prediction_method = c("ARIMA", "VAR", "ets"),
                                 select_K = c("Fixed", "EVR"), K = 6) {
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
    lambda_val <- med_polish_resi_eigen$values
    retain_component <- select_k(tau = 10^-2, eigenvalue = lambda_val)
  }
  
  var_total_variations <- (sum(med_polish_resi_eigen$values[1:retain_component]) / sum(med_polish_resi_eigen$values)) * 100
  med_polish_resi_basis <- as.matrix(med_polish_resi_eigen$vectors[, 1:retain_component])
  med_polish_resi_score <- crossprod(med_polish_resi, med_polish_resi_basis)
  
  med_polish_resi_score_forecast <- matrix(NA, retain_component, fh)
  
  if (prediction_method == "ets") {
    for (ip in 1:retain_component) {
      dum <- forecast:::forecast.ets(object = ets(y = as.vector(med_polish_resi_score[, ip])), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip, ] <- as.vector(dum$mean)
    }
  }
  
  med_polish_resi_forecast  <- med_polish_resi_basis %*% med_polish_resi_score_forecast
  Fixed                     <- t(fixed_com)[, 1:fh]
  med_polish_curve_forecast <- med_polish_resi_forecast + Fixed
  
  return(list(
    med_polish_curve_forecast = med_polish_curve_forecast,
    med_polish_resi_forecast  = med_polish_resi_forecast,
    TV                        = var_total_variations
  ))
}

# Forecasting execution loop
Boot_female_EVR <- list()
Boot_female_K   <- list()
Boot_male_EVR   <- list()
Boot_male_K     <- list()

max_h   <- 17
no_core <- detectCores() - 2

for (i in 1:n_states) {
  pref           <- (n_year * i - (n_year - 1)):(n_year * i)
  n_training_ini <- length(pref) - max_h
  fixed_com1     <- Fixed_part_means[pref, ]
  Residuals_f1   <- Residuals_mean[pref, ]
  
  registerDoMC(no_core)
  
  boot_EVR <- foreach(iwk = 1:max_h) %dopar% Pref_forecast_curves(
    fixed_com = fixed_com1[1:(n_training_ini - 1 + iwk), ],
    Residuals_f = Residuals_f1[1:(n_training_ini - 1 + iwk), ],
    est_method = "cov", fh = (max_h - iwk + 1), B = 1000,
    prediction_method = "ets", select_K = "EVR", K = 6
  )$med_polish_curve_forecast
  
  boot_K <- foreach(iwk = 1:max_h) %dopar% Pref_forecast_curves(
    fixed_com = fixed_com1[1:(n_training_ini - 1 + iwk), ],
    Residuals_f = Residuals_f1[1:(n_training_ini - 1 + iwk), ],
    est_method = "cov", fh = (max_h - iwk + 1), B = 1000,
    prediction_method = "ets", select_K = "Fixed", K = 6
  )$med_polish_curve_forecast
  
  Forecast_female_EVR <- list()
  Forecast_male_EVR   <- list()
  Forecast_female_K   <- list()
  Forecast_male_K     <- list()
  
  for (h in seq_along(boot_K)) {
    Forecast_female_EVR[[h]] <- boot_EVR[[h]][111:220, ]
    Forecast_male_EVR[[h]]   <- boot_EVR[[h]][1:110, ]
    Forecast_female_K[[h]]   <- boot_K[[h]][111:220, ]
    Forecast_male_K[[h]]     <- boot_K[[h]][1:110, ]
  }
  
  Forecast_female_transformed_EVR <- list()
  Forecast_male_transformed_EVR   <- list()
  Forecast_female_transformed_K   <- list()
  Forecast_male_transformed_K     <- list()
  
  for (h in seq_along(Forecast_female_EVR)) {
    fore_val_EVR_1 <- as.matrix(Forecast_female_EVR[[h]])
    fore_val_EVR_2 <- as.matrix(Forecast_male_EVR[[h]])
    fore_val_K_1   <- as.matrix(Forecast_female_K[[h]])
    fore_val_K_2   <- as.matrix(Forecast_male_K[[h]])
    
    n_age_sub <- nrow(fore_val_EVR_1)
    n_cols    <- ncol(fore_val_EVR_1)
    
    f_x_t_star_fore_EVR_1 <- f_x_t_star_fore_EVR_2 <- d_x_t_star_fore_EVR_1 <- d_x_t_star_fore_EVR_2 <- matrix(NA, (n_age_sub + 1), n_cols)
    f_x_t_star_fore_K_1   <- f_x_t_star_fore_K_2   <- d_x_t_star_fore_K_1   <- d_x_t_star_fore_K_2   <- matrix(NA, (n_age_sub + 1), n_cols)
    
    for (ik in 1:n_cols) {
      f_x_t_star_fore_EVR_1[, ik] <- c(invlogit(fore_val_EVR_1[, ik]), 1)
      f_x_t_star_fore_EVR_2[, ik] <- c(invlogit(fore_val_EVR_2[, ik]), 1)
      d_x_t_star_fore_EVR_1[, ik] <- c(f_x_t_star_fore_EVR_1[1, ik], diff(f_x_t_star_fore_EVR_1[, ik]))
      d_x_t_star_fore_EVR_2[, ik] <- c(f_x_t_star_fore_EVR_2[1, ik], diff(f_x_t_star_fore_EVR_2[, ik]))
      
      f_x_t_star_fore_K_1[, ik]   <- c(invlogit(fore_val_K_1[, ik]), 1)
      f_x_t_star_fore_K_2[, ik]   <- c(invlogit(fore_val_K_2[, ik]), 1)
      d_x_t_star_fore_K_1[, ik]   <- c(f_x_t_star_fore_K_1[1, ik], diff(f_x_t_star_fore_K_1[, ik]))
      d_x_t_star_fore_K_2[, ik]   <- c(f_x_t_star_fore_K_2[1, ik], diff(f_x_t_star_fore_K_2[, ik]))
    }
    
    Forecast_female_transformed_EVR[[h]] <- d_x_t_star_fore_EVR_1 * 10^5
    Forecast_male_transformed_EVR[[h]]   <- d_x_t_star_fore_EVR_2 * 10^5
    Forecast_female_transformed_K[[h]]   <- d_x_t_star_fore_K_1 * 10^5
    Forecast_male_transformed_K[[h]]     <- d_x_t_star_fore_K_2 * 10^5
  }
  
  Boot_female_EVR[[i]] <- Forecast_female_transformed_EVR
  Boot_female_K[[i]]   <- Forecast_female_transformed_K
  Boot_male_EVR[[i]]   <- Forecast_male_transformed_EVR
  Boot_male_K[[i]]     <- Forecast_male_transformed_K
}

# Error Metrics Evaluation
compute_error <- function(forecast, true) {
  c(KLD(forecast, true), JSD_geom(forecast, true))
}

pairwise_errors <- function(j, forecast, true) {
  forecast <- forecast[, j]
  if (is.null(dim(true))) true <- matrix(true, ncol = 1)
  true <- true[, j]
  compute_error(forecast, true)
}

create_new_list <- function(A) {
  num_columns <- 17
  B <- vector("list", num_columns)
  for (col in 1:num_columns) {
    B[[col]] <- sapply(A, function(mat) {
      if (ncol(mat) >= col) mat[, col] else rep(NA, nrow(mat))
    })
  }
  return(B)
}

errors_KLD_female_EVR <- errors_KLD_male_EVR <- matrix(NA, 17, 47)
errors_JSD_female_EVR <- errors_JSD_male_EVR <- matrix(NA, 17, 47)
errors_KLD_female_K   <- errors_KLD_male_K   <- matrix(NA, 17, 47)
errors_JSD_female_K   <- errors_JSD_male_K   <- matrix(NA, 17, 47)

for (i in 1:n_states) {
  forecast_female_transformed_EVR <- Boot_female_EVR[[i]]
  forecast_female_transformed_K   <- Boot_female_K[[i]]
  forecast_male_transformed_EVR   <- Boot_male_EVR[[i]]
  forecast_male_transformed_K     <- Boot_male_K[[i]]
  
  Holdout_female <- t(female_prefecture_dx[[i]])
  Holdout_male   <- t(male_prefecture_dx[[i]])
  
  Holdout_female[Holdout_female == 0] <- 10^-5
  Holdout_male[Holdout_male == 0]     <- 10^-5
  
  errors_female_EVR_list <- errors_male_EVR_list <- vector("list", 17)
  errors_female_K_list   <- errors_male_K_list   <- vector("list", 17)
  
  for (h in seq_along(forecast_female_transformed_EVR)) {
    forecast_female_EVR <- as.matrix(forecast_female_transformed_EVR[[h]])
    forecast_male_EVR   <- as.matrix(forecast_male_transformed_EVR[[h]])
    forecast_female_K   <- as.matrix(forecast_female_transformed_K[[h]])
    forecast_male_K     <- as.matrix(forecast_male_transformed_K[[h]])
    
    holdout_cols <- (33 + (h - 1)):49
    true_female  <- as.matrix(Holdout_female[, holdout_cols])
    true_male    <- as.matrix(Holdout_male[, holdout_cols])
    
    errors_female_EVR_list[[h]] <- sapply(seq_len(ncol(forecast_female_EVR)), function(j) pairwise_errors(j, forecast_female_EVR, true_female))
    errors_male_EVR_list[[h]]   <- sapply(seq_len(ncol(forecast_male_EVR)), function(j) pairwise_errors(j, forecast_male_EVR, true_male))
    errors_female_K_list[[h]]   <- sapply(seq_len(ncol(forecast_female_K)), function(j) pairwise_errors(j, forecast_female_K, true_female))
    errors_male_K_list[[h]]     <- sapply(seq_len(ncol(forecast_male_K)), function(j) pairwise_errors(j, forecast_male_K, true_male))
  }
  
  new_female_EVR <- create_new_list(errors_female_EVR_list)
  new_female_K   <- create_new_list(errors_female_K_list)
  new_male_EVR   <- create_new_list(errors_male_EVR_list)
  new_male_K     <- create_new_list(errors_male_K_list)
  
  errors_female_EVR_final <- errors_male_EVR_final <- matrix(NA, 2, 17)
  errors_female_K_final   <- errors_male_K_final   <- matrix(NA, 2, 17)
  
  for (h in 1:17) {
    errors_female_EVR_final[, h] <- rowMeans(new_female_EVR[[h]], na.rm = TRUE)
    errors_male_EVR_final[, h]   <- rowMeans(new_male_EVR[[h]], na.rm = TRUE)
    errors_female_K_final[, h]   <- rowMeans(new_female_K[[h]], na.rm = TRUE)
    errors_male_K_final[, h]     <- rowMeans(new_male_K[[h]], na.rm = TRUE)
  }
  
  errors_KLD_female_EVR[, i] <- errors_female_EVR_final[1, ]
  errors_KLD_male_EVR[, i]   <- errors_male_EVR_final[1, ]
  errors_JSD_female_EVR[, i] <- errors_female_EVR_final[2, ]
  errors_JSD_male_EVR[, i]   <- errors_male_EVR_final[2, ]
  
  errors_KLD_female_K[, i]   <- errors_female_K_final[1, ]
  errors_KLD_male_K[, i]     <- errors_male_K_final[1, ]
  errors_JSD_female_K[, i]   <- errors_female_K_final[2, ]
  errors_JSD_male_K[, i]     <- errors_male_K_final[2, ]
}

errors_pref_EVR <- list(errors_KLD_female_EVR, errors_JSD_female_EVR, errors_KLD_male_EVR, errors_JSD_male_EVR)
errors_pref_K   <- list(errors_KLD_female_K, errors_JSD_female_K, errors_KLD_male_K, errors_JSD_male_K)

errors_EVR <- cbind(apply(errors_KLD_female_EVR, 1, mean), apply(errors_JSD_female_EVR, 1, mean), apply(errors_KLD_male_EVR, 1, mean), apply(errors_JSD_male_EVR, 1, mean))
colnames(errors_EVR) <- c("KLD_Female", "JSD_Female", "KLD_Male", "JSD_Male")

errors_K <- cbind(apply(errors_KLD_female_K, 1, mean), apply(errors_JSD_female_K, 1, mean), apply(errors_KLD_male_K, 1, mean), apply(errors_JSD_male_K, 1, mean))
colnames(errors_K) <- c("KLD_Female", "JSD_Female", "KLD_Male", "JSD_Male")

# Save primary output matrices (Prefixed with FANOVA_)
saveRDS(errors_KLD_female_EVR, file = file.path(dir_results, "FANOVA_errors_KLD_female_EVR.rds"))
saveRDS(errors_JSD_female_EVR, file = file.path(dir_results, "FANOVA_errors_JSD_female_EVR.rds"))
saveRDS(errors_KLD_male_EVR,   file = file.path(dir_results, "FANOVA_errors_KLD_male_EVR.rds"))
saveRDS(errors_JSD_male_EVR,   file = file.path(dir_results, "FANOVA_errors_JSD_male_EVR.rds"))

saveRDS(errors_KLD_female_K,   file = file.path(dir_results, "FANOVA_errors_KLD_female_K.rds"))
saveRDS(errors_JSD_female_K,   file = file.path(dir_results, "FANOVA_errors_JSD_female_K.rds"))
saveRDS(errors_KLD_male_K,     file = file.path(dir_results, "FANOVA_errors_KLD_male_K.rds"))
saveRDS(errors_JSD_male_K,     file = file.path(dir_results, "FANOVA_errors_JSD_male_K.rds"))

saveRDS(errors_pref_EVR, file = file.path(dir_results, "FANOVA_errors_pref_EVR_list.rds"))
saveRDS(errors_pref_K,   file = file.path(dir_results, "FANOVA_errors_pref_K_list.rds"))

saveRDS(errors_EVR, file = file.path(dir_results, "FANOVA_errors_EVR_matrix.rds"))
saveRDS(errors_K,   file = file.path(dir_results, "FANOVA_errors_K_matrix.rds"))

saveRDS(colMeans(errors_EVR), file = file.path(dir_results, "FANOVA_errors_EVR_colMeans.rds"))
saveRDS(apply(errors_EVR, 2, median), file = file.path(dir_results, "FANOVA_errors_EVR_median.rds"))

saveRDS(colMeans(errors_K), file = file.path(dir_results, "FANOVA_errors_K_colMeans.rds"))
saveRDS(apply(errors_K, 2, median), file = file.path(dir_results, "FANOVA_errors_K_median.rds"))

# Prepare interactive Shiny App exports
fanova_list <- list(
  errors_KLD_female_EVR = errors_KLD_female_EVR,
  errors_JSD_female_EVR = errors_JSD_female_EVR,
  errors_KLD_male_EVR   = errors_KLD_male_EVR,
  errors_JSD_male_EVR   = errors_JSD_male_EVR,
  errors_KLD_female_K   = errors_KLD_female_K,
  errors_JSD_female_K   = errors_JSD_female_K,
  errors_KLD_male_K     = errors_KLD_male_K,
  errors_JSD_male_K     = errors_JSD_male_K
)

for (name in names(fanova_list)) {
  mat        <- fanova_list[[name]]
  err_metric <- ifelse(grepl("KLD", name), "KLD", "JSD")
  gender     <- ifelse(grepl("female", name), "female", "male")
  comp       <- ifelse(grepl("EVR", name), "EVR", "K")
  
  file_name  <- paste0("FANOVA_", err_metric, "_", gender, "_", comp, ".rds")
  saveRDS(mat, file = file.path(dir_shiny, file_name))
}

