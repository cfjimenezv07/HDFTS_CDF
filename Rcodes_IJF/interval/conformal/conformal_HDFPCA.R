# ==============================================================================
# Directory Setup & Source Dependencies
# ==============================================================================
source("load_packages.R")
source(here("auxiliary_source", "auxiliary_interval.R"))

# --- Define Output Directories ---
dir.results <- here("results", "Results_Interval", "Conformal")
dir.shiny2  <- here("results", "Shiny_App", "datasets_shiny_app", "IFE")

# Create directories if they do not exist
if (!dir.exists(dir.results)) dir.create(dir.results, recursive = TRUE)
if (!dir.exists(dir.shiny2))  dir.create(dir.shiny2, recursive = TRUE)



# ==============================================================================
# Main Split Conformal Prediction Function (HDFPCA)
# ==============================================================================

hdfpca_interval_fore_subnational_cdf_conformal <- function(fdata_F_data, fdata_M_data, fore_method, horizon, level_sig) {
  n_age <- ncol(fdata_F_data)
  val_steps  <- 18 - horizon
  test_steps <- 17 - horizon
  
  # ----------------------------------------------------------------------------
  # Validation Phase
  # ----------------------------------------------------------------------------
  if (fore_method == "CDF") {
    dum <- hdfpca_fun_fore(
      fdata_F = fdata_F_data[1:33, ],
      fdata_M = fdata_M_data[1:33, ],
      horizon = horizon,
      first_order = 6,
      second_order = 2,
      transformation = "CDF"
    )
    forecast_validation_F <- dum$forecast_val_F
    forecast_validation_M <- dum$forecast_val_M
  } else if (fore_method == "CLR") {
    forecast_validation_F <- forecast_validation_M <- matrix(NA, n_age, val_steps)
    for (ij in 1:val_steps) {
      dum <- clr_MFTS_fun(
        fdata_F = fdata_F_data[1:(15 + ij), ],
        fdata_M = fdata_M_data[1:(15 + ij), ],
        ncomp_selection = "EVR",
        fh = horizon,
        fore_method = "ets"
      )
      forecast_validation_F[, ij] <- dum$MFTS_res_fore_F
      forecast_validation_M[, ij] <- dum$MFTS_res_fore_M
    }
  } else {
    stop("Forecasting method must either be 'CDF' or 'CLR'.")
  }
  
  # Holdout validation data & Residuals Computation
  holdout_validation_F <- t(fdata_F_data[(16 + horizon):33, , drop = FALSE])
  holdout_validation_M <- t(fdata_M_data[(16 + horizon):33, , drop = FALSE])
  
  resi_mat_F <- holdout_validation_F - forecast_validation_F
  resi_mat_M <- holdout_validation_M - forecast_validation_M
  
  quantile_resid_F <- apply(resi_mat_F, 1, function(x) quantile(abs(x), probs = level_sig))
  quantile_resid_M <- apply(resi_mat_M, 1, function(x) quantile(abs(x), probs = level_sig))
  
  # ----------------------------------------------------------------------------
  # Testing Phase
  # ----------------------------------------------------------------------------
  if (fore_method == "CDF") {
    dum <- hdfpca_fun_fore(
      fdata_F = fdata_F_data,
      fdata_M = fdata_M_data,
      horizon = horizon,
      first_order = 6,
      second_order = 2,
      transformation = "CDF_test"
    )
    forecast_test_F <- dum$forecast_val_F_test
    forecast_test_F_lb <- forecast_test_F - quantile_resid_F
    forecast_test_F_ub <- forecast_test_F + quantile_resid_F
    
    forecast_test_M <- dum$forecast_val_M_test
    forecast_test_M_lb <- forecast_test_M - quantile_resid_M
    forecast_test_M_ub <- forecast_test_M + quantile_resid_M
  } else if (fore_method == "CLR") {
    forecast_test_F <- forecast_test_F_lb <- forecast_test_F_ub <- matrix(NA, n_age, test_steps)
    forecast_test_M <- forecast_test_M_lb <- forecast_test_M_ub <- matrix(NA, n_age, test_steps)
    
    for (ij in 1:test_steps) {
      dum <- clr_MFTS_fun(
        fdata_F = fdata_F_data[1:(31 + ij), ],
        fdata_M = fdata_M_data[1:(31 + ij), ],
        ncomp_selection = "EVR",
        fh = horizon,
        fore_method = "ets"
      )
      
      forecast_test_F[, ij]    <- dum$MFTS_res_fore_F
      forecast_test_F_lb[, ij] <- forecast_test_F[, ij] - quantile_resid_F
      forecast_test_F_ub[, ij] <- forecast_test_F[, ij] + quantile_resid_F
      
      forecast_test_M[, ij]    <- dum$MFTS_res_fore_M
      forecast_test_M_lb[, ij] <- forecast_test_M[, ij] - quantile_resid_M
      forecast_test_M_ub[, ij] <- forecast_test_M[, ij] + quantile_resid_M
    }
  }
  
  # Interval Evaluation via holdout testing data
  holdout_val_F <- t(fdata_F_data[(33 + horizon):49, , drop = FALSE])
  holdout_val_M <- t(fdata_M_data[(33 + horizon):49, , drop = FALSE])
  
  int_F_err <- interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
  int_M_err <- interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
  
  return(list(int_F_err = int_F_err, int_M_err = int_M_err))
}


# ==============================================================================
# Simulation Routine Wrapper
# ==============================================================================

run_hdfpca_conformal_experiment <- function(level_sig) {
  err_F <- err_M <- array(NA, dim = c(16, 3, 47))
  
  for (ij in 1:47) {
    for (iw in 1:16) {
      dum <- hdfpca_interval_fore_subnational_cdf_conformal(
        fdata_F_data = female_prefecture_dx[[ij]],
        fdata_M_data = male_prefecture_dx[[ij]],
        fore_method  = "CDF",
        horizon      = iw,
        level_sig    = level_sig
      )
      
      err_F[iw, , ij] <- dum$int_F_err
      err_M[iw, , ij] <- dum$int_M_err
    }
    cat("Completed prefecture unit:", ij, "\n")
  }
  
  # Calculate aggregated means
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
res_HDFPCA_080 <- run_hdfpca_conformal_experiment(level_sig = 0.80)

# Level of significance = 0.95
res_HDFPCA_095 <- run_hdfpca_conformal_experiment(level_sig = 0.95)


# ==============================================================================
# Save Outputs to dir.results
# ==============================================================================

save_hdfpca_results <- function(res, suffix = "") {
  # Arrays
  saveRDS(res$err_F, file.path(dir.results, paste0("hdfpca_int_fore_subnational_err_F_EVR_conformal", suffix, ".rds")))
  saveRDS(res$err_M, file.path(dir.results, paste0("hdfpca_int_fore_subnational_err_M_EVR_conformal", suffix, ".rds")))
  
  # Means
  saveRDS(res$mean_F, file.path(dir.results, paste0("hdfpca_int_fore_subnational_err_F_EVR_conformal", suffix, "_mean.rds")))
  saveRDS(res$mean_M, file.path(dir.results, paste0("hdfpca_int_fore_subnational_err_M_EVR_conformal", suffix, "_mean.rds")))
}

# --- Save level_sig = 0.80 ---
save_hdfpca_results(res_HDFPCA_080, suffix = "")

# --- Save level_sig = 0.95 ---
save_hdfpca_results(res_HDFPCA_095, suffix = "_alpha_0.95")


# ==============================================================================
# Export to Shiny Directory
# ==============================================================================

arrays_list_shiny <- list(
  HDFPCA_female_EVR_80 = res_HDFPCA_080$err_F,
  HDFPCA_male_EVR_80   = res_HDFPCA_080$err_M,
  HDFPCA_female_EVR_95 = res_HDFPCA_095$err_F,
  HDFPCA_male_EVR_95   = res_HDFPCA_095$err_M
)

for (name in names(arrays_list_shiny)) {
  file_name <- paste0(name, ".rds")
  saveRDS(arrays_list_shiny[[name]], file = file.path(dir.shiny2, file_name))
}
