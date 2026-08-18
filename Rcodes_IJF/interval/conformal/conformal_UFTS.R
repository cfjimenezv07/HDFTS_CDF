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
# Main Split Conformal Prediction Function
# ==============================================================================

interval_fore_subnational_cdf_conformal <- function(fdata, method_ncomp, horizon, fore_method, uni_fore_method, level_sig) {
  n_age <- ncol(fdata)
  val_steps <- 18 - horizon
  test_steps <- 17 - horizon
  
  fore_validation <- matrix(NA, n_age, val_steps)
  
  # ----------------------------------------------------------------------------
  # Validation Phase
  # ----------------------------------------------------------------------------
  if (fore_method == "CDF") {
    for (ij in 1:val_steps) {
      fore_validation[, ij] <- fore_national_cdf(data_set = fdata[1:(15 + ij), ],
                                                 ncomp_method = method_ncomp,
                                                 fh = horizon, fmethod = uni_fore_method)
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:val_steps) {
      fore_validation[, ij] <- as.numeric(clr_fun(fdata = fdata[1:(15 + ij), ],
                                                  ncomp_selection = method_ncomp,
                                                  fh = horizon, fore_method = "ETS")$fore_count)
    }
  } else {
    stop("fore_method must be either 'CDF' or 'CLR'.")
  }
  
  # Holdout validation data & Residuals Computation
  holdout_validation_dum <- t(fdata[(16 + horizon):33, , drop = FALSE])
  resi_mat <- holdout_validation_dum - fore_validation
  
  # Extract absolute empirical quantiles
  quantile_resid <- apply(resi_mat, 1, function(x) quantile(abs(x), probs = level_sig))
  
  # ----------------------------------------------------------------------------
  # Testing Phase
  # ----------------------------------------------------------------------------
  fore_val <- fore_val_lb <- fore_val_ub <- matrix(NA, n_age, test_steps)
  
  if (fore_method == "CDF") {
    for (ij in 1:test_steps) {
      fore_val[, ij] <- fore_national_cdf(data_set = fdata[1:(32 + ij), ], ncomp_method = method_ncomp,
                                          fh = horizon, fmethod = uni_fore_method)
      fore_val_lb[, ij] <- fore_val[, ij] - quantile_resid
      fore_val_ub[, ij] <- fore_val[, ij] + quantile_resid
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:test_steps) {
      fore_val[, ij] <- as.numeric(clr_fun(fdata = fdata[1:(32 + ij), ], ncomp_selection = method_ncomp,
                                           fh = horizon, fore_method = "ETS")$fore_count)
      fore_val_lb[, ij] <- fore_val[, ij] - quantile_resid
      fore_val_ub[, ij] <- fore_val[, ij] + quantile_resid
    }
  }
  
  # Interval Evaluation via holdout testing data
  holdout_val_dum <- t(fdata[(33 + horizon):49, , drop = FALSE])
  int_err <- interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, alpha = (1 - level_sig))
  
  return(list(int_err = int_err, dimn = ncol(resi_mat)))
}


# ==============================================================================
# Simulation Routine Wrapper
# ==============================================================================

run_conformal_experiment <- function(method_ncomp, level_sig) {
  err_F <- err_M <- array(NA, dim = c(16, 3, 47))
  
  for (ij in 1:47) {
    for (iw in 1:16) {
      # Female Eval
      res_F <- interval_fore_subnational_cdf_conformal(
        fdata = female_prefecture_dx[[ij]], method_ncomp = method_ncomp,
        horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = level_sig
      )
      err_F[iw, , ij] <- res_F$int_err
      
      # Male Eval
      res_M <- interval_fore_subnational_cdf_conformal(
        fdata = male_prefecture_dx[[ij]], method_ncomp = method_ncomp,
        horizon = iw, fore_method = "CDF", uni_fore_method = "ets", level_sig = level_sig
      )
      err_M[iw, , ij] <- res_M$int_err
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
# Execute Experiments (CDF)
# ==============================================================================

# Note: method_ncomp = "provide" is equivalent to the K = 6 condition

# Alpha = 0.80
res_EVR_080 <- run_conformal_experiment(method_ncomp = "EVR", level_sig = 0.80)
res_K6_080  <- run_conformal_experiment(method_ncomp = "provide", level_sig = 0.80)

# Alpha = 0.95
res_EVR_095 <- run_conformal_experiment(method_ncomp = "EVR", level_sig = 0.95)
res_K6_095  <- run_conformal_experiment(method_ncomp = "provide", level_sig = 0.95)


# ==============================================================================
# Save Outputs to dir.results
# ==============================================================================

save_results <- function(res, method_label, suffix) {
  # Arrays
  saveRDS(res$err_F, file.path(dir.results, paste0("int_fore_subnational_err_F_", method_label, "_conformal", suffix, ".rds")))
  saveRDS(res$err_M, file.path(dir.results, paste0("int_fore_subnational_err_M_", method_label, "_conformal", suffix, ".rds")))
  
  # Means
  saveRDS(res$mean_F, file.path(dir.results, paste0("int_fore_subnational_err_F_", method_label, "_conformal", suffix, "_mean.rds")))
  saveRDS(res$mean_M, file.path(dir.results, paste0("int_fore_subnational_err_M_", method_label, "_conformal", suffix, "_mean.rds")))
}

# --- level_sig = 0.80 (No alpha suffix in your original logic) ---
save_results(res_EVR_080, "EVR", "")
save_results(res_K6_080,  "K6",  "")

# --- level_sig = 0.95 ---
save_results(res_EVR_095, "EVR", "_alpha_0.95")
save_results(res_K6_095,  "K6",  "_alpha_0.95")


# ==============================================================================
# (Optional) Export to Shiny Directory if needed for UFTS
# ==============================================================================

arrays_list_shiny <- list(
  UFTS_female_EVR_80 = res_EVR_080$err_F,
  UFTS_male_EVR_80   = res_EVR_080$err_M,
  UFTS_female_K_80   = res_K6_080$err_F,
  UFTS_male_K_80     = res_K6_080$err_M,
  
  UFTS_female_EVR_95 = res_EVR_095$err_F,
  UFTS_male_EVR_95   = res_EVR_095$err_M,
  UFTS_female_K_95   = res_K6_095$err_F,
  UFTS_male_K_95     = res_K6_095$err_M
)

for (name in names(arrays_list_shiny)) {
  file_name <- paste0(name, ".rds")
  saveRDS(arrays_list_shiny[[name]], file = file.path(dir.shiny2, file_name))
}

