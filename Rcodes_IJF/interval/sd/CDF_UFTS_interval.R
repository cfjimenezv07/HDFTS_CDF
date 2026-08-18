###############################################################################
# INTERVAL FORECASTING: SD / UFTS METHOD (CDF)
###############################################################################

# --- Source Auxiliary Functions ---
source("load_packages.R")
source(here("auxiliary_source", "auxiliary_interval.R"))

# --- Define Output Directories ---
dir.results <- here("results", "Results_Interval", "sd")
dir.shiny2  <- here("Shiny_App", "datasets_shiny_app", "IFE")

# Create directories if they do not exist
if (!dir.exists(dir.results)) dir.create(dir.results, recursive = TRUE)
if (!dir.exists(dir.shiny2))  dir.create(dir.shiny2, recursive = TRUE)

###############################################################################
# 1. FUNCTION DEFINITION
###############################################################################

# train_set: 1:16 | validation_set: 17:32 | test_set: 33:49
interval_fore_subnational_cdf <- function(fdata, method_ncomp, horizon, fore_method, 
                                          uni_fore_method, level_sig, type_PI) 
{
  n_age <- ncol(fdata)
  fore_validation <- matrix(NA, ncol(fdata), (18 - horizon))
  
  if (fore_method == "CDF") {
    for (ij in 1:(18 - horizon)) {
      fore_validation[, ij] <- fore_national_cdf(data_set = fdata[1:(15 + ij), ], 
                                                 ncomp_method = method_ncomp,
                                                 fh = horizon, 
                                                 fmethod = uni_fore_method)
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:(18 - horizon)) {
      fore_validation[, ij] <- as.numeric(clr_fun(fdata = fdata[1:(15 + ij), ], 
                                                  ncomp_selection = method_ncomp,
                                                  fh = horizon, 
                                                  fore_method = "ETS")$fore_count)
    }
  } else {
    stop("Forecasting method must either be CDF or CLR.")
  }
  
  # Holdout validation data
  holdout_validation_dum <- t(matrix(fdata[(16 + horizon):33, ], length((16 + horizon):33), ncol(fdata)))
  resi_mat <- holdout_validation_dum - fore_validation
  
  # Compute standard deviation of residuals across validation set
  sd_val_input <- apply(resi_mat, 1, sd)
  
  # Find optimal tuning parameter using multiple optimization routines
  tune_para_find_val_1 <- optimise(f = tune_para_find_function, interval = c(0, 1), 
                                   resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                   alpha_level = level_sig, PI_type = type_PI)
  
  tune_para_find_val_2 <- optimise(f = tune_para_find_function, interval = c(0, 5), 
                                   resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                   alpha_level = level_sig, PI_type = type_PI)
  
  tune_para_find_val_3 <- optimise(f = tune_para_find_function, interval = c(0, 10), 
                                   resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                   alpha_level = level_sig, PI_type = type_PI)
  
  tune_para_find_val_4 <- optimise(f = tune_para_find_function, interval = c(0, 20), 
                                   resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                   alpha_level = level_sig, PI_type = type_PI)
  
  tune_para_find_val_5 <- optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B", 
                                resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                alpha_level = level_sig, PI_type = type_PI)
  
  tune_para_find_val_6 <- optim(par = 1, fn = tune_para_find_function, method = "Nelder-Mead",
                                resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                alpha_level = level_sig, PI_type = type_PI)
  
  obj_val <- c(tune_para_find_val_1$objective, 
               tune_para_find_val_2$objective, 
               tune_para_find_val_3$objective,
               tune_para_find_val_4$objective,
               tune_para_find_val_5$value,
               tune_para_find_val_6$value)
  
  obj_val_min <- min(obj_val)
  
  tune_para_find <- c(tune_para_find_val_1$minimum, 
                      tune_para_find_val_2$minimum,
                      tune_para_find_val_3$minimum,
                      tune_para_find_val_4$minimum,
                      tune_para_find_val_5$par,
                      tune_para_find_val_6$par)[which.min(obj_val)]
  
  # Compute test set bounds
  fore_val <- fore_val_lb <- fore_val_ub <- matrix(NA, ncol(fdata), (17 - horizon))
  
  if (fore_method == "CDF") {
    for (ij in 1:(17 - horizon)) {
      fore_val[, ij] <- fore_national_cdf(data_set = fdata[1:(32 + ij), ], 
                                          ncomp_method = method_ncomp,
                                          fh = horizon, 
                                          fmethod = uni_fore_method)
      fore_val_lb[, ij] <- fore_val[, ij] - tune_para_find * sd_val_input
      fore_val_ub[, ij] <- fore_val[, ij] + tune_para_find * sd_val_input
    }
  } else if (fore_method == "CLR") {
    for (ij in 1:(17 - horizon)) {
      fore_val[, ij] <- as.numeric(clr_fun(fdata = fdata[1:(32 + ij), ], 
                                           ncomp_selection = method_ncomp,
                                           fh = horizon, 
                                           fore_method = "ETS")$fore_count)
      fore_val_lb[, ij] <- fore_val[, ij] - tune_para_find * sd_val_input
      fore_val_ub[, ij] <- fore_val[, ij] + tune_para_find * sd_val_input
    }
  }
  
  # Holdout testing data
  holdout_val_dum <- t(matrix(fdata[(33 + horizon):49, ], length((33 + horizon):49), ncol(fdata)))
  
  if (type_PI == "pointwise") {
    int_err <- interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, 
                              alpha = (1 - level_sig))
  } else if (type_PI == "uniform") {
    int_err <- uniform_cpd(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub,
                           alpha = (1 - level_sig))
  }
  
  return(list(int_err = int_err, tune_para_find = tune_para_find, tune_para_find_obj = obj_val_min))
}

###############################################################################
# 2. COMPUTATION LOOP FUNCTION
###############################################################################

run_interval_loop <- function(ncomp_method, alpha_sig) {
  err_F <- err_M <- array(NA, dim = c(16, 3, 47))
  tune_F <- tune_obj_F <- tune_M <- tune_obj_M <- matrix(NA, 16, 47)
  
  for (ij in 1:47) {
    for (iw in 1:16) {
      # Female (F)
      dum_F <- interval_fore_subnational_cdf(fdata = female_prefecture_dx[[ij]], 
                                             method_ncomp = ncomp_method,
                                             horizon = iw, fore_method = "CDF", 
                                             uni_fore_method = "ets",
                                             level_sig = alpha_sig, type_PI = "pointwise")
      err_F[iw, , ij]      <- dum_F$int_err
      tune_F[iw, ij]      <- dum_F$tune_para_find
      tune_obj_F[iw, ij]  <- dum_F$tune_para_find_obj
      
      # Male (M)
      dum_M <- interval_fore_subnational_cdf(fdata = male_prefecture_dx[[ij]], 
                                             method_ncomp = ncomp_method,
                                             horizon = iw, fore_method = "CDF", 
                                             uni_fore_method = "ets",
                                             level_sig = alpha_sig, type_PI = "pointwise")
      err_M[iw, , ij]      <- dum_M$int_err
      tune_M[iw, ij]      <- dum_M$tune_para_find
      tune_obj_M[iw, ij]  <- dum_M$tune_para_find_obj
    }
    cat(sprintf("Completed Prefecture %d/47 for Method: %s | Alpha: %.2f\n", ij, ncomp_method, alpha_sig))
  }
  
  # Calculate averages
  mean_F <- apply(err_F, c(1, 2), mean)
  mean_M <- apply(err_M, c(1, 2), mean)
  
  colnames(mean_F) <- colnames(mean_M) <- c("ECP", "CPD", "score")
  rownames(mean_F) <- rownames(mean_M) <- 1:16
  
  return(list(err_F = err_F, err_M = err_M, 
              tune_F = tune_F, tune_M = tune_M, 
              tune_obj_F = tune_obj_F, tune_obj_M = tune_obj_M,
              mean_F = mean_F, mean_M = mean_M))
}

###############################################################################
# 3. RUN EXPERIMENTS & SAVE RESULTS
###############################################################################

# --- Config 1: EVR (level_sig = 0.80) ---
res_EVR_80 <- run_interval_loop(ncomp_method = "EVR", alpha_sig = 0.80)

saveRDS(res_EVR_80$err_F,      file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS.rds"))
saveRDS(res_EVR_80$err_M,      file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS.rds"))
saveRDS(res_EVR_80$tune_F,     file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_tune_para.rds"))
saveRDS(res_EVR_80$tune_M,     file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_tune_para.rds"))
saveRDS(res_EVR_80$tune_obj_F, file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_tune_para_obj.rds"))
saveRDS(res_EVR_80$tune_obj_M, file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_tune_para_obj.rds"))
saveRDS(res_EVR_80$mean_F,     file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_mean.rds"))
saveRDS(res_EVR_80$mean_M,     file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_mean.rds"))

# --- Config 2: K = 6 (level_sig = 0.80) ---
res_K6_80 <- run_interval_loop(ncomp_method = "provide", alpha_sig = 0.80)

saveRDS(res_K6_80$err_F,      file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS.rds"))
saveRDS(res_K6_80$err_M,      file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS.rds"))
saveRDS(res_K6_80$tune_F,     file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_tune_para.rds"))
saveRDS(res_K6_80$tune_M,     file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_tune_para.rds"))
saveRDS(res_K6_80$tune_obj_F, file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_tune_para_obj.rds"))
saveRDS(res_K6_80$tune_obj_M, file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_tune_para_obj.rds"))
saveRDS(res_K6_80$mean_F,     file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_mean.rds"))
saveRDS(res_K6_80$mean_M,     file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_mean.rds"))

# --- Config 3: EVR (level_sig = 0.95) ---
res_EVR_95 <- run_interval_loop(ncomp_method = "EVR", alpha_sig = 0.95)

saveRDS(res_EVR_95$err_F,      file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_alpha_0.95.rds"))
saveRDS(res_EVR_95$err_M,      file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_alpha_0.95.rds"))
saveRDS(res_EVR_95$tune_F,     file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95.rds"))
saveRDS(res_EVR_95$tune_M,     file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95.rds"))
saveRDS(res_EVR_95$tune_obj_F, file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(res_EVR_95$tune_obj_M, file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(res_EVR_95$mean_F,     file = paste0(dir.results, "int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean.rds"))
saveRDS(res_EVR_95$mean_M,     file = paste0(dir.results, "int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean.rds"))

# --- Config 4: K = 6 (level_sig = 0.95) ---
res_K6_95 <- run_interval_loop(ncomp_method = "provide", alpha_sig = 0.95)

saveRDS(res_K6_95$err_F,      file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_alpha_0.95.rds"))
saveRDS(res_K6_95$err_M,      file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_alpha_0.95.rds"))
saveRDS(res_K6_95$tune_F,     file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95.rds"))
saveRDS(res_K6_95$tune_M,     file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95.rds"))
saveRDS(res_K6_95$tune_obj_F, file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(res_K6_95$tune_obj_M, file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(res_K6_95$mean_F,     file = paste0(dir.results, "int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean.rds"))
saveRDS(res_K6_95$mean_M,     file = paste0(dir.results, "int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean.rds"))

###############################################################################
# 4. EXPORT TO SHINY APP DIRECTORY (FORMAT: 47 x 16 x 3)
###############################################################################

arrays_list <- list(
  F_EVR_80 = res_EVR_80$err_F,
  M_EVR_80 = res_EVR_80$err_M,
  F_K6_80  = res_K6_80$err_F,
  M_K6_80  = res_K6_80$err_M,
  F_EVR_95 = res_EVR_95$err_F,
  M_EVR_95 = res_EVR_95$err_M,
  F_K6_95  = res_K6_95$err_F,
  M_K6_95  = res_K6_95$err_M
)

for (name in names(arrays_list)) {
  arr <- arrays_list[[name]]
  
  # Reorder dimensions: 16 x 3 x 47 -> 47 x 16 x 3
  arr_reordered <- aperm(arr, c(3, 1, 2))
  
  gender <- ifelse(grepl("^F", name), "female", "male")
  comp   <- ifelse(grepl("EVR", name), "EVR", "K")
  level  <- ifelse(grepl("80", name), "80", "95")
  
  file_name <- paste0("UFTS_", gender, "_", comp, "_", level, ".rds")
  saveRDS(arr_reordered, file = file.path(dir.shiny2, file_name))
}

