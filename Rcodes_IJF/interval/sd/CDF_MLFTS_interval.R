# ==============================================================================
# 0. SETUP & DEPENDENCIES
# ==============================================================================
source("load_packages.R")

# Directories
dir_source  <- here("auxiliary_source")
dir_results <- here("results", "Results_Interval", "sd")
dir.shiny2  <- here("results", "Shiny_App", "datasets_shiny_app", "IFE")

# Source auxiliary scripts
source(file.path(dir_source, "auxiliary_interval.R"))
source(file.path(dir_source, "auxiliary_point.R"))

# Create output directories if they don't exist
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(dir.shiny2, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

#' Optimize Tuning Parameter across standard intervals and methods
optimize_tuning_parameter <- function(resi_mat, sd_val, alpha_level, PI_type) {
  opt1 <- optimise(tune_para_find_function, interval = c(0, 1), 
                   resi_mat = resi_mat, sd_val_input = sd_val, alpha_level = alpha_level, PI_type = PI_type)
  opt2 <- optimise(tune_para_find_function, interval = c(0, 5), 
                   resi_mat = resi_mat, sd_val_input = sd_val, alpha_level = alpha_level, PI_type = PI_type)
  opt3 <- optimise(tune_para_find_function, interval = c(0, 10), 
                   resi_mat = resi_mat, sd_val_input = sd_val, alpha_level = alpha_level, PI_type = PI_type)
  opt4 <- optimise(tune_para_find_function, interval = c(0, 20), 
                   resi_mat = resi_mat, sd_val_input = sd_val, alpha_level = alpha_level, PI_type = PI_type)
  
  opt5 <- optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B", 
                resi_mat = resi_mat, sd_val_input = sd_val, alpha_level = alpha_level, PI_type = PI_type)
  opt6 <- optim(par = 1, fn = tune_para_find_function, 
                resi_mat = resi_mat, sd_val_input = sd_val, alpha_level = alpha_level, PI_type = PI_type)
  
  objectives <- c(opt1$objective, opt2$objective, opt3$objective, opt4$objective, opt5$value, opt6$value)
  candidates <- c(opt1$minimum, opt2$minimum, opt3$minimum, opt4$minimum, opt5$par, opt6$par)
  
  best_idx <- which.min(objectives)
  
  list(
    tune_para = candidates[best_idx],
    obj_min   = objectives[best_idx]
  )
}


#' Compute Interval Forecast Error Metrics for Subnational MLFTS
interval_fore_subnational_cdf_MLFTS <- function(fdata_F, fdata_M, fdata_common, fore_method, 
                                                horizon, way_ncomp, level_sig, type_PI) {
  
  n_age <- ncol(fdata_F)
  n_val_steps <- 18 - horizon
  
  forecast_val_F <- matrix(NA, nrow = n_age, ncol = n_val_steps)
  forecast_val_M <- matrix(NA, nrow = n_age, ncol = n_val_steps)
  
  # Validation Forecast Phase
  if (fore_method == "CDF") {
    for (ij in 1:n_val_steps) {
      dum <- fore_national_cdf_MLFTS(
        data_set_F   = fdata_F[1:(15 + ij), ], 
        data_set_M   = fdata_M[1:(15 + ij), ],
        aux_variable = fdata_common[1:(15 + ij), ], 
        fh           = horizon, 
        fmethod      = "ets",
        method_ncomp = way_ncomp
      )
      forecast_val_F[, ij] <- dum$mlfts_fore_F
      forecast_val_M[, ij] <- dum$mlfts_fore_M
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:n_val_steps) {
      dum <- clr_MLFTS_fun(
        fdata_F         = fdata_F[1:(15 + ij), ], 
        fdata_M         = fdata_M[1:(15 + ij), ],
        ncomp_selection = way_ncomp, 
        fh              = horizon
      )
      forecast_val_F[, ij] <- dum$MLFTS_res_fore_F
      forecast_val_M[, ij] <- dum$MLFTS_res_fore_M
    }
  } else {
    stop("Forecasting method must be either 'CDF' or 'CLR'.")
  }
  
  # Holdout Validation Residuals
  holdout_val_F <- t(fdata_F[(16 + horizon):33, , drop = FALSE])
  holdout_val_M <- t(fdata_M[(16 + horizon):33, , drop = FALSE])
  
  resi_mat_F <- holdout_val_F - forecast_val_F
  resi_mat_M <- holdout_val_M - forecast_val_M
  
  sd_val_F <- apply(resi_mat_F, 1, sd)
  sd_val_M <- apply(resi_mat_M, 1, sd)
  
  # Find Optimal Tuning Parameters
  opt_F <- optimize_tuning_parameter(resi_mat_F, sd_val_F, level_sig, type_PI)
  opt_M <- optimize_tuning_parameter(resi_mat_M, sd_val_M, level_sig, type_PI)
  
  # Test Forecast Phase
  n_test_steps <- 17 - horizon
  forecast_test_F <- forecast_test_F_lb <- forecast_test_F_ub <- matrix(NA, nrow = n_age, ncol = n_test_steps)
  forecast_test_M <- forecast_test_M_lb <- forecast_test_M_ub <- matrix(NA, nrow = n_age, ncol = n_test_steps)
  
  if (fore_method == "CDF") {
    for (ij in 1:n_test_steps) {
      dum <- fore_national_cdf_MLFTS(
        data_set_F   = fdata_F[1:(32 + ij), ], 
        data_set_M   = fdata_M[1:(32 + ij), ],
        aux_variable = fdata_common[1:(32 + ij), ], 
        fh           = horizon, 
        fmethod      = "ets",
        method_ncomp = way_ncomp
      )
      forecast_test_F[, ij]    <- dum$mlfts_fore_F
      forecast_test_F_lb[, ij] <- forecast_test_F[, ij] - opt_F$tune_para * sd_val_F
      forecast_test_F_ub[, ij] <- forecast_test_F[, ij] + opt_F$tune_para * sd_val_F
      
      forecast_test_M[, ij]    <- dum$mlfts_fore_M
      forecast_test_M_lb[, ij] <- forecast_test_M[, ij] - opt_M$tune_para * sd_val_M
      forecast_test_M_ub[, ij] <- forecast_test_M[, ij] + opt_M$tune_para * sd_val_M
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:n_test_steps) {
      dum <- clr_MLFTS_fun(
        fdata_F         = fdata_F[1:(32 + ij), ], 
        fdata_M         = fdata_M[1:(32 + ij), ],
        ncomp_selection = way_ncomp, 
        fh              = horizon
      )
      forecast_test_F[, ij]    <- dum$MLFTS_res_fore_F
      forecast_test_F_lb[, ij] <- forecast_test_F[, ij] - opt_F$tune_para * sd_val_F
      forecast_test_F_ub[, ij] <- forecast_test_F[, ij] + opt_F$tune_para * sd_val_F
      
      forecast_test_M[, ij]    <- dum$MLFTS_res_fore_M
      forecast_test_M_lb[, ij] <- forecast_test_M[, ij] - opt_M$tune_para * sd_val_M
      forecast_test_M_ub[, ij] <- forecast_test_M[, ij] + opt_M$tune_para * sd_val_M
    }
  }
  
  # Final Evaluation on Holdout Testing Set
  holdout_test_F <- t(fdata_F[(33 + horizon):49, , drop = FALSE])
  holdout_test_M <- t(fdata_M[(33 + horizon):49, , drop = FALSE])
  
  alpha_param <- 1 - level_sig
  if (type_PI == "pointwise") {
    int_F_err <- interval_score(holdout = holdout_test_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = alpha_param)
    int_M_err <- interval_score(holdout = holdout_test_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = alpha_param)
  } else if (type_PI == "uniform") {
    int_F_err <- uniform_cpd(holdout = holdout_test_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = alpha_param)
    int_M_err <- uniform_cpd(holdout = holdout_test_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = alpha_param)
  }
  
  list(
    int_F_err            = int_F_err, 
    tune_para_find_F     = opt_F$tune_para,
    tune_para_find_F_obj = opt_F$obj_min,
    int_M_err            = int_M_err, 
    tune_para_find_M     = opt_M$tune_para,
    tune_para_find_M_obj = opt_M$obj_min
  )
}


# ==============================================================================
# 2. BATCH SIMULATION RUNNER
# ==============================================================================

run_simulation <- function(way_ncomp, level_sig, num_pref = 47, horizons = 16) {
  err_F <- err_M <- array(NA, dim = c(horizons, 3, num_pref))
  tune_F <- tune_obj_F <- tune_M <- tune_obj_M <- matrix(NA, horizons, num_pref)
  
  for (ij in 1:num_pref) {
    for (iw in 1:horizons) {
      dum <- interval_fore_subnational_cdf_MLFTS(
        fdata_F      = female_prefecture_dx[[ij]],
        fdata_M      = male_prefecture_dx[[ij]],
        fdata_common = NULL,
        fore_method  = "CDF",
        horizon      = iw,
        way_ncomp    = way_ncomp,
        level_sig    = level_sig,
        type_PI      = "pointwise"
      )
      
      err_F[iw, , ij]    <- dum$int_F_err
      tune_F[iw, ij]     <- dum$tune_para_find_F
      tune_obj_F[iw, ij] <- dum$tune_para_find_F_obj
      
      err_M[iw, , ij]    <- dum$int_M_err
      tune_M[iw, ij]     <- dum$tune_para_find_M
      tune_obj_M[iw, ij] <- dum$tune_para_find_M_obj
    }
    cat(sprintf("Completed Prefecture %d/%d for Config [%s | Sig: %.2f]\n", ij, num_pref, way_ncomp, level_sig))
  }
  
  err_F_mean <- apply(err_F, c(1, 2), mean)
  err_M_mean <- apply(err_M, c(1, 2), mean)
  
  colnames(err_F_mean) <- colnames(err_M_mean) <- c("ECP", "CPD", "score")
  rownames(err_F_mean) <- rownames(err_M_mean) <- 1:horizons
  
  list(
    err_F = err_F, err_M = err_M,
    tune_F = tune_F, tune_M = tune_M,
    tune_obj_F = tune_obj_F, tune_obj_M = tune_obj_M,
    err_F_mean = err_F_mean, err_M_mean = err_M_mean
  )
}

scenarios <- list(
  EVR_80 = list(way_ncomp = "EVR",     level_sig = 0.80, suffix = "_EVR_ETS"),
  K6_80  = list(way_ncomp = "provide", level_sig = 0.80, suffix = "_K6_ETS"),
  EVR_95 = list(way_ncomp = "EVR",     level_sig = 0.95, suffix = "_EVR_ETS_alpha_0.95"),
  K6_95  = list(way_ncomp = "provide", level_sig = 0.95, suffix = "_K6_ETS_alpha_0.95")
)

results <- lapply(scenarios, function(s) {
  run_simulation(way_ncomp = s$way_ncomp, level_sig = s$level_sig)
})


# ==============================================================================
# 3. SAVE RESULTS TO RDS
# ==============================================================================

for (key in names(scenarios)) {
  sfx <- scenarios[[key]]$suffix
  res <- results[[key]]
  
  saveRDS(res$err_F,      file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_F", sfx, ".rds")))
  saveRDS(res$err_M,      file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_M", sfx, ".rds")))
  saveRDS(res$tune_F,     file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_F", sfx, "_tune_para.rds")))
  saveRDS(res$tune_M,     file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_M", sfx, "_tune_para.rds")))
  saveRDS(res$tune_obj_F, file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_F", sfx, "_tune_para_obj.rds")))
  saveRDS(res$tune_obj_M, file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_M", sfx, "_tune_para_obj.rds")))
  saveRDS(res$err_F_mean, file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_F", sfx, "_mean.rds")))
  saveRDS(res$err_M_mean, file = file.path(dir_results, paste0("MLFTS_int_fore_subnational_err_M", sfx, "_mean.rds")))
}


# ==============================================================================
# 4. EXPORT TO SHINY APP FORMAT
# ==============================================================================

shiny_export_map <- list(
  F_EVR_80 = results$EVR_80$err_F,
  M_EVR_80 = results$EVR_80$err_M,
  F_K6_80  = results$K6_80$err_F,
  M_K6_80  = results$K6_80$err_M,
  F_EVR_95 = results$EVR_95$err_F,
  M_EVR_95 = results$EVR_95$err_M,
  F_K6_95  = results$K6_95$err_F,
  M_K6_95  = results$K6_95$err_M
)

for (name in names(shiny_export_map)) {
  arr <- shiny_export_map[[name]]
  arr_reordered <- aperm(arr, c(3, 1, 2))
  
  gender <- ifelse(startsWith(name, "F"), "female", "male")
  comp   <- ifelse(grepl("EVR", name), "EVR", "K")
  level  <- ifelse(grepl("80", name), "80", "95")
  
  file_name <- sprintf("MLFTS_%s_%s_%s.rds", gender, comp, level)
  saveRDS(arr_reordered, file = file.path(dir.shiny2, file_name))
}