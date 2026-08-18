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
# Main Split Conformal Prediction Function (MLFTS)
# ==============================================================================

interval_fore_subnational_cdf_MLFTS_conformal <- function(fdata_F, fdata_M, fdata_common, fore_method, horizon,
                                                          way_ncomp, level_sig, type_PI = NULL) {
  n_age <- ncol(fdata_F)
  val_steps  <- 18 - horizon
  test_steps <- 17 - horizon
  
  forecast_validation_F <- forecast_validation_M <- matrix(NA, n_age, val_steps)
  
  # ----------------------------------------------------------------------------
  # Validation Phase
  # ----------------------------------------------------------------------------
  if (fore_method == "CDF") {
    for (ij in 1:val_steps) {
      dum <- fore_national_cdf_MLFTS(
        data_set_F = fdata_F[1:(15 + ij), ],
        data_set_M = fdata_M[1:(15 + ij), ],
        aux_variable = if (!is.null(fdata_common)) fdata_common[1:(15 + ij), ] else NULL,
        fh = horizon,
        fmethod = "ets",
        method_ncomp = way_ncomp
      )
      forecast_validation_F[, ij] <- dum$mlfts_fore_F
      forecast_validation_M[, ij] <- dum$mlfts_fore_M
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:val_steps) {
      dum <- clr_MLFTS_fun(
        fdata_F = fdata_F[1:(15 + ij), ],
        fdata_M = fdata_M[1:(15 + ij), ],
        ncomp_selection = way_ncomp,
        fh = horizon
      )
      forecast_validation_F[, ij] <- dum$MLFTS_res_fore_F
      forecast_validation_M[, ij] <- dum$MLFTS_res_fore_M
    }
  } else {
    stop("Forecasting method must either be 'CDF' or 'CLR'.")
  }
  
  rownames(forecast_validation_F) <- rownames(forecast_validation_M) <- 1:n_age
  colnames(forecast_validation_F) <- colnames(forecast_validation_M) <- 1:val_steps
  
  # Holdout validation data & Residuals Computation
  holdout_validation_F <- t(fdata_F[(16 + horizon):33, , drop = FALSE])
  holdout_validation_M <- t(fdata_M[(16 + horizon):33, , drop = FALSE])
  
  resi_mat_F <- holdout_validation_F - forecast_validation_F
  resi_mat_M <- holdout_validation_M - forecast_validation_M
  
  quantile_resid_F <- apply(resi_mat_F, 1, function(x) quantile(abs(x), probs = level_sig))
  quantile_resid_M <- apply(resi_mat_M, 1, function(x) quantile(abs(x), probs = level_sig))
  
  # ----------------------------------------------------------------------------
  # Testing Phase
  # ----------------------------------------------------------------------------
  forecast_test_F <- forecast_test_F_lb <- forecast_test_F_ub <- matrix(NA, n_age, test_steps)
  forecast_test_M <- forecast_test_M_lb <- forecast_test_M_ub <- matrix(NA, n_age, test_steps)
  
  if (fore_method == "CDF") {
    for (ij in 1:test_steps) {
      dum <- fore_national_cdf_MLFTS(
        data_set_F = fdata_F[1:(32 + ij), ],
        data_set_M = fdata_M[1:(32 + ij), ],
        aux_variable = if (!is.null(fdata_common)) fdata_common[1:(32 + ij), ] else NULL,
        fh = horizon,
        fmethod = "ets",
        method_ncomp = way_ncomp
      )
      
      forecast_test_F[, ij]    <- dum$mlfts_fore_F
      forecast_test_F_lb[, ij] <- forecast_test_F[, ij] - quantile_resid_F
      forecast_test_F_ub[, ij] <- forecast_test_F[, ij] + quantile_resid_F
      
      forecast_test_M[, ij]    <- dum$mlfts_fore_M
      forecast_test_M_lb[, ij] <- forecast_test_M[, ij] - quantile_resid_M
      forecast_test_M_ub[, ij] <- forecast_test_M[, ij] + quantile_resid_M
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:test_steps) {
      dum <- clr_MLFTS_fun(
        fdata_F = fdata_F[1:(32 + ij), ],
        fdata_M = fdata_M[1:(32 + ij), ],
        ncomp_selection = way_ncomp,
        fh = horizon
      )
      
      forecast_test_F[, ij]    <- dum$MLFTS_res_fore_F
      forecast_test_F_lb[, ij] <- forecast_test_F[, ij] - quantile_resid_F
      forecast_test_F_ub[, ij] <- forecast_test_F[, ij] + quantile_resid_F
      
      forecast_test_M[, ij]    <- dum$MLFTS_res_fore_M
      forecast_test_M_lb[, ij] <- forecast_test_M[, ij] - quantile_resid_M
      forecast_test_M_ub[, ij] <- forecast_test_M[, ij] + quantile_resid_M
    }
  }
  
  # Interval Evaluation via holdout testing data
  holdout_val_F <- t(fdata_F[(33 + horizon):49, , drop = FALSE])
  holdout_val_M <- t(fdata_M[(33 + horizon):49, , drop = FALSE])
  
  int_F_err <- interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
  int_M_err <- interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
  
  return(list(int_F_err = int_F_err, int_M_err = int_M_err))
}


# ==============================================================================
# Simulation Routine Wrapper
# ==============================================================================

run_mlfts_conformal_experiment <- function(way_ncomp, level_sig) {
  err_F <- err_M <- array(NA, dim = c(16, 3, 47))
  
  for (ij in 1:47) {
    for (iw in 1:16) {
      dum <- interval_fore_subnational_cdf_MLFTS_conformal(
        fdata_F = female_prefecture_dx[[ij]],
        fdata_M = male_prefecture_dx[[ij]],
        fdata_common = NULL,
        fore_method = "CDF",
        horizon = iw,
        way_ncomp = way_ncomp,
        level_sig = level_sig
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
res_EVR_080 <- run_mlfts_conformal_experiment(way_ncomp = "EVR",     level_sig = 0.80)
res_K6_080  <- run_mlfts_conformal_experiment(way_ncomp = "provide", level_sig = 0.80)

# Level of significance = 0.95
res_EVR_095 <- run_mlfts_conformal_experiment(way_ncomp = "EVR",     level_sig = 0.95)
res_K6_095  <- run_mlfts_conformal_experiment(way_ncomp = "provide", level_sig = 0.95)


# ==============================================================================
# Save Outputs to dir.results
# ==============================================================================

save_mlfts_results <- function(res, method_label, suffix) {
  # Arrays
  saveRDS(res$err_F, file.path(dir.results, paste0("MLFTS_int_fore_subnational_err_F_", method_label, "_conformal", suffix, ".rds")))
  saveRDS(res$err_M, file.path(dir.results, paste0("MLFTS_int_fore_subnational_err_M_", method_label, "_conformal", suffix, ".rds")))
  
  # Means
  saveRDS(res$mean_F, file.path(dir.results, paste0("MLFTS_int_fore_subnational_err_F_", method_label, "_conformal", suffix, "_mean.rds")))
  saveRDS(res$mean_M, file.path(dir.results, paste0("MLFTS_int_fore_subnational_err_M_", method_label, "_conformal", suffix, "_mean.rds")))
}

# --- level_sig = 0.80 ---
save_mlfts_results(res_EVR_080, "EVR", "")
save_mlfts_results(res_K6_080,  "K6",  "")

# --- level_sig = 0.95 ---
save_mlfts_results(res_EVR_095, "EVR", "_alpha_0.95")
save_mlfts_results(res_K6_095,  "K6",  "_alpha_0.95")


# ==============================================================================
# Export to Shiny Directory
# ==============================================================================

arrays_list_shiny <- list(
  MLFTS_female_EVR_80 = res_EVR_080$err_F,
  MLFTS_male_EVR_80   = res_EVR_080$err_M,
  MLFTS_female_K_80   = res_K6_080$err_F,
  MLFTS_male_K_80     = res_K6_080$err_M,
  
  MLFTS_female_EVR_95 = res_EVR_095$err_F,
  MLFTS_male_EVR_95   = res_EVR_095$err_M,
  MLFTS_female_K_95   = res_K6_095$err_F,
  MLFTS_male_K_95     = res_K6_095$err_M
)

for (name in names(arrays_list_shiny)) {
  file_name <- paste0(name, ".rds")
  saveRDS(arrays_list_shiny[[name]], file = file.path(dir.shiny2, file_name))
}