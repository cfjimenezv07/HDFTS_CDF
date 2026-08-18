# ==============================================================================
# File Directory Setup
# ==============================================================================
source("load_packages.R")

dir_source  <- here("auxiliary_source")
dir_results <- here("results", "Results_Interval", "sd")
dir.shiny2  <- here("results", "Shiny_App", "datasets_shiny_app", "IFE")

# Load dependencies using relative pathing
source(file.path(dir_source, "auxiliary_interval.R"))

# Ensure destination export directories exist
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)
if (!dir.exists(dir.shiny2))  dir.create(dir.shiny2, recursive = TRUE)

# fdata_F: female data
# fdata_M: male data
# horizon: forecast horizon
# first_order: number of components
# second_order: number of components
# transformation: type of transformation
# level_sig: level of significance

int_hdfpca_fun <- function(fdata_F, fdata_M, horizon, first_order, second_order, transformation, level_sig)
{
  n_age = ncol(fdata_F)
  forecast_validation_F = forecast_validation_M = matrix(NA, ncol(fdata_F), (18 - horizon))
  if(transformation == "direct")
  {
    for(ij in 1:(18 - horizon))
    {
      data_F = fdata_F[1:(15+ij),]
      data_M = fdata_M[1:(15+ij),]
      data_comb = list()
      data_comb[[1]] = t(data_F)
      data_comb[[2]] = t(data_M)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      forecast_validation_F[,ij] = (fore_val[[1]])[,horizon]
      forecast_validation_M[,ij] = (fore_val[[2]])[,horizon]
      rm(ij); rm(fore_val)
    }
  }
  else if(transformation == "CLR")
  {
    for(ij in 1:(18 - horizon))
    {
      data_F = as.matrix(clr(fdata_F[1:(15+ij),]))
      data_M = as.matrix(clr(fdata_M[1:(15+ij),]))
      data_comb = list()
      data_comb[[1]] = t(data_F)
      data_comb[[2]] = t(data_M)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      forecast_validation_F[,ij] = as.numeric(clrInv((fore_val[[1]])[,horizon])) * 10^5
      forecast_validation_M[,ij] = as.numeric(clrInv((fore_val[[2]])[,horizon])) * 10^5
      rm(ij); rm(fore_val)
    }	
  }
  else if(transformation == "CDF")
  {
    for(ijk in 1:(18 - horizon))
    {
      data_F = fdata_F[1:(15+ijk),]/10^5
      data_M = fdata_M[1:(15+ijk),]/10^5
      
      data_F_cumsum_dum = data_M_cumsum_dum = matrix(NA, nrow(data_F), ncol(data_F))
      for(iw in 1:nrow(data_F))
      {
        data_F_cumsum_dum[iw,] = cumsum(data_F[iw,])
        data_M_cumsum_dum[iw,] = cumsum(data_M[iw,])
        rm(iw)
      }
      
      # check if any cumsum values equal to 0
      if(any(data_F_cumsum_dum == 0))
      {
        data_F_cumsum = replace(data_F_cumsum_dum, which(data_F_cumsum_dum == 0), 10^-5)
      }
      else
      {
        data_F_cumsum = data_F_cumsum_dum
      }
      
      if(any(data_M_cumsum_dum == 0))
      {
        data_M_cumsum = replace(data_M_cumsum_dum, which(data_M_cumsum_dum == 0), 10^-5)
      }
      else
      {
        data_M_cumsum = data_M_cumsum_dum
      }
      rm(data_F_cumsum_dum); rm(data_M_cumsum_dum)
      
      # logit transformation
      
      data_F_cumsum_logit = data_M_cumsum_logit = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
      for(ij in 1:nrow(data_F))
      {
        data_F_cumsum_logit[ij,] = logit(data_F_cumsum[ij, 1:(ncol(data_F) - 1)])
        data_M_cumsum_logit[ij,] = logit(data_M_cumsum[ij, 1:(ncol(data_M) - 1)])
        rm(ij)
      }
      
      data_comb = list()
      data_comb[[1]] = t(data_F_cumsum_logit)
      data_comb[[2]] = t(data_M_cumsum_logit)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      fore_val_F = (fore_val[[1]])[,horizon]
      fore_val_M = (fore_val[[2]])[,horizon]
      
      data_cumsum_logit_fore_add_F = c(invlogit(fore_val_F), 1)
      data_cumsum_logit_fore_add_M = c(invlogit(fore_val_M), 1)
      
      data_cumsum_logit_fore_add_diff_F = c(data_cumsum_logit_fore_add_F[1], diff(data_cumsum_logit_fore_add_F))
      data_cumsum_logit_fore_add_diff_M = c(data_cumsum_logit_fore_add_M[1], diff(data_cumsum_logit_fore_add_M))
      
      forecast_validation_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
      forecast_validation_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
      rm(ijk); rm(data_F); rm(data_M); rm(fore_val)
    }
  }
  else
  {
    warning("none, CLR, CDF transformation allowed only.")
  }
  
  # holdout validation data
  
  holdout_validation_F = t(matrix(fdata_F[(16 + horizon):33,], length((16 + horizon):33), ncol(fdata_F)))
  holdout_validation_M = t(matrix(fdata_M[(16 + horizon):33,], length((16 + horizon):33), ncol(fdata_M)))
  resi_mat_F = holdout_validation_F - forecast_validation_F
  resi_mat_M = holdout_validation_M - forecast_validation_M
  
  # compute standard deviation of residuals
  
  sd_val_F = apply(resi_mat_F, 1, sd)
  sd_val_M = apply(resi_mat_M, 1, sd)
  
  # find the optimal tuning parameter
  
  tune_para_find_val_F_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                    resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_F_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                    resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_F_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                    resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_F_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                    resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_F_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                 resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                 alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_F_6 = optim(par = 1, fn = tune_para_find_function, 
                                 resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                 alpha_level = level_sig, PI_type = "pointwise")
  
  obj_val = c(tune_para_find_val_F_1$objective, 
              tune_para_find_val_F_2$objective, 
              tune_para_find_val_F_3$objective, 
              tune_para_find_val_F_4$objective, 
              tune_para_find_val_F_5$value, 
              tune_para_find_val_F_6$value)
  obj_val_min_F = min(obj_val)
  
  tune_para_find_F = c(tune_para_find_val_F_1$minimum, 
                       tune_para_find_val_F_2$minimum, 
                       tune_para_find_val_F_3$minimum, 
                       tune_para_find_val_F_4$minimum, 
                       tune_para_find_val_F_5$par,
                       tune_para_find_val_F_6$par)[which.min(obj_val)]
  rm(obj_val)
  
  # male
  
  tune_para_find_val_M_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                    resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_M_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                    resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_M_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                    resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_M_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                    resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                    alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_M_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                 resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                 alpha_level = level_sig, PI_type = "pointwise")
  
  tune_para_find_val_M_6 = optim(par = 1, fn = tune_para_find_function, 
                                 resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                 alpha_level = level_sig, PI_type = "pointwise")
  
  obj_val = c(tune_para_find_val_M_1$objective, 
              tune_para_find_val_M_2$objective, 
              tune_para_find_val_M_3$objective, 
              tune_para_find_val_M_4$objective, 
              tune_para_find_val_M_5$value, 
              tune_para_find_val_M_6$value)
  obj_val_min_M = min(obj_val)
  
  tune_para_find_M = c(tune_para_find_val_M_1$minimum, 
                       tune_para_find_val_M_2$minimum, 
                       tune_para_find_val_M_3$minimum, 
                       tune_para_find_val_M_4$minimum, 
                       tune_para_find_val_M_5$par,
                       tune_para_find_val_M_6$par)[which.min(obj_val)]
  rm(obj_val)
  
  forecast_test_F = forecast_test_F_lb = forecast_test_F_ub = matrix(NA, ncol(fdata_F), (17 - horizon))
  forecast_test_M = forecast_test_M_lb = forecast_test_M_ub = matrix(NA, ncol(fdata_M), (17 - horizon))
  if(transformation == "direct")
  {
    for(ij in 1:(17 - horizon))
    {
      data_F = fdata_F[1:(32+ij),]
      data_M = fdata_M[1:(32+ij),]
      data_comb = list()
      data_comb[[1]] = t(data_F)
      data_comb[[2]] = t(data_M)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      forecast_test_F[,ij] = (fore_val[[1]])[,horizon]
      forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
      forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
      
      forecast_test_M[,ij] = (fore_val[[2]])[,horizon]
      forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
      forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
      rm(ij)
    }		
  }
  else if(transformation == "CLR")
  {
    for(ij in 1:(17 - horizon))
    {
      data_F = as.matrix(clr(fdata_F[1:(32+ij),]))
      data_M = as.matrix(clr(fdata_M[1:(32+ij),]))
      data_comb = list()
      data_comb[[1]] = t(data_F)
      data_comb[[2]] = t(data_M)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      forecast_test_F[,ij] = as.numeric(clrInv((fore_val[[1]])[,horizon])) * 10^5
      forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
      forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
      
      forecast_test_M[,ij] = as.numeric(clrInv((fore_val[[2]])[,horizon])) * 10^5
      forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
      forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
      rm(ij); rm(fore_val)
    }	
  }
  else if(transformation == "CDF")
  {
    for(ijk in 1:(17 - horizon))
    {
      data_F = fdata_F[1:(32+ijk),]/10^5
      data_M = fdata_M[1:(32+ijk),]/10^5
      
      data_F_cumsum_dum = data_M_cumsum_dum = matrix(NA, nrow(data_F), ncol(data_F))
      for(iw in 1:nrow(data_F))
      {
        data_F_cumsum_dum[iw,] = cumsum(data_F[iw,])
        data_M_cumsum_dum[iw,] = cumsum(data_M[iw,])
        rm(iw)
      }
      
      # check if any cumsum values equal to 0
      if(any(data_F_cumsum_dum == 0))
      {
        data_F_cumsum = replace(data_F_cumsum_dum, which(data_F_cumsum_dum == 0), 10^-5)
      }
      else
      {
        data_F_cumsum = data_F_cumsum_dum
      }
      
      if(any(data_M_cumsum_dum == 0))
      {
        data_M_cumsum = replace(data_M_cumsum_dum, which(data_M_cumsum_dum == 0), 10^-5)
      }
      else
      {
        data_M_cumsum = data_M_cumsum_dum
      }
      rm(data_F_cumsum_dum); rm(data_M_cumsum_dum)
      
      # logit transformation
      
      data_F_cumsum_logit = data_M_cumsum_logit = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
      for(ij in 1:nrow(data_F))
      {
        data_F_cumsum_logit[ij,] = logit(data_F_cumsum[ij, 1:(ncol(data_F) - 1)])
        data_M_cumsum_logit[ij,] = logit(data_M_cumsum[ij, 1:(ncol(data_M) - 1)])
        rm(ij)
      }
      
      data_comb = list()
      data_comb[[1]] = t(data_F_cumsum_logit)
      data_comb[[2]] = t(data_M_cumsum_logit)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      fore_val_F = (fore_val[[1]])[,horizon]
      fore_val_M = (fore_val[[2]])[,horizon]
      
      data_cumsum_logit_fore_add_F = c(invlogit(fore_val_F), 1)
      data_cumsum_logit_fore_add_M = c(invlogit(fore_val_M), 1)
      
      data_cumsum_logit_fore_add_diff_F = c(data_cumsum_logit_fore_add_F[1], diff(data_cumsum_logit_fore_add_F))
      data_cumsum_logit_fore_add_diff_M = c(data_cumsum_logit_fore_add_M[1], diff(data_cumsum_logit_fore_add_M))
      
      forecast_test_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
      forecast_test_F_lb[,ijk] = forecast_test_F[,ijk] - tune_para_find_F * sd_val_F
      forecast_test_F_ub[,ijk] = forecast_test_F[,ijk] + tune_para_find_F * sd_val_F
      
      forecast_test_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
      forecast_test_M_lb[,ijk] = forecast_test_M[,ijk] - tune_para_find_M * sd_val_M
      forecast_test_M_ub[,ijk] = forecast_test_M[,ijk] + tune_para_find_M * sd_val_M
      rm(ijk); rm(data_F); rm(data_M)
    }
  }	
  
  # holdout testing data
  
  holdout_val_F = t(matrix(fdata_F[(33 + horizon):49,], length((33 + horizon):49), ncol(fdata_F)))
  holdout_val_M = t(matrix(fdata_M[(33 + horizon):49,], length((33 + horizon):49), ncol(fdata_M)))
  
  int_F_err = interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
  int_M_err = interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
  return(list(int_F_err = int_F_err, tune_para_find_F = tune_para_find_F, 
              tune_para_find_F_obj = obj_val_min_F,
              
              int_M_err = int_M_err, tune_para_find_M = tune_para_find_M,
              tune_para_find_M_obj = obj_val_min_M))    
}

# Helper function to encapsulate the prefecture loops
run_eval_loop <- function(level_sig, first_order, second_order) {
  err_F <- err_M <- array(NA, dim = c(47, 16, 3), dimnames = list(state, 1:16, c("ECP", "CPD", "score")))
  tune_F <- tune_F_obj <- tune_M <- tune_M_obj <- matrix(NA, 47, 16)
  
  for(ij in 1:47) {
    for(iw in 1:16) {
      dum <- int_hdfpca_fun(
        fdata_F = female_prefecture_dx[[ij]], 
        fdata_M = male_prefecture_dx[[ij]], 
        horizon = iw, 
        first_order = first_order, 
        second_order = second_order, 
        transformation = "CDF", 
        level_sig = level_sig
      )
      
      err_F[ij, iw, ]      <- dum$int_F_err
      tune_F[ij, iw]       <- dum$tune_para_find_F
      tune_F_obj[ij, iw]   <- dum$tune_para_find_F_obj
      
      err_M[ij, iw, ]      <- dum$int_M_err
      tune_M[ij, iw]       <- dum$tune_para_find_M
      tune_M_obj[ij, iw]   <- dum$tune_para_find_M_obj
      rm(dum)
    }
  }
  
  horizon_F <- apply(err_F, c(2, 3), mean)
  horizon_M <- apply(err_M, c(2, 3), mean)
  
  return(list(
    err_F = err_F, tune_F = tune_F, tune_F_obj = tune_F_obj, horizon_F = horizon_F,
    err_M = err_M, tune_M = tune_M, tune_M_obj = tune_M_obj, horizon_M = horizon_M
  ))
}

# ==============================================================================
# COMPUTATIONS & SAVING (EVR vs. K6)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. EVR Model (first_order = 6, second_order = 2)
# ------------------------------------------------------------------------------

# alpha = 0.8
res_EVR_80 <- run_eval_loop(level_sig = 0.8, first_order = 6, second_order = 2)

saveRDS(res_EVR_80$err_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_CDF.rds"))
saveRDS(res_EVR_80$err_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_CDF.rds"))
saveRDS(res_EVR_80$tune_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_CDF.rds"))
saveRDS(res_EVR_80$tune_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_CDF.rds"))
saveRDS(res_EVR_80$tune_F_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_CDF.rds"))
saveRDS(res_EVR_80$tune_M_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_CDF.rds"))
saveRDS(res_EVR_80$horizon_F, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_CDF.rds"))
saveRDS(res_EVR_80$horizon_M, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_CDF.rds"))

# alpha = 0.95
res_EVR_95 <- run_eval_loop(level_sig = 0.95, first_order = 6, second_order = 2)

saveRDS(res_EVR_95$err_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$err_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$tune_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$tune_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$tune_F_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$tune_M_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$horizon_F, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_CDF.rds"))
saveRDS(res_EVR_95$horizon_M, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_CDF.rds"))


# ------------------------------------------------------------------------------
# 2. K6 Model (Fixed-K model specifications)
# Note: Adjust parameters if K6 requires different order/r parameters than EVR
# ------------------------------------------------------------------------------

# alpha = 0.8
res_K6_80 <- run_eval_loop(level_sig = 0.8, first_order = 6, second_order = 2)

saveRDS(res_K6_80$err_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_CDF.rds"))
saveRDS(res_K6_80$err_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_CDF.rds"))
saveRDS(res_K6_80$tune_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_tune_para_CDF.rds"))
saveRDS(res_K6_80$tune_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_tune_para_CDF.rds"))
saveRDS(res_K6_80$tune_F_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_tune_para_obj_CDF.rds"))
saveRDS(res_K6_80$tune_M_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_tune_para_obj_CDF.rds"))
saveRDS(res_K6_80$horizon_F, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_F_K6_ETS_CDF.rds"))
saveRDS(res_K6_80$horizon_M, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_M_K6_ETS_CDF.rds"))

# alpha = 0.95
res_K6_95 <- run_eval_loop(level_sig = 0.95, first_order = 6, second_order = 2)

saveRDS(res_K6_95$err_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$err_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$tune_F, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$tune_M, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$tune_F_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$tune_M_obj, file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$horizon_F, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_F_K6_ETS_alpha_0.95_CDF.rds"))
saveRDS(res_K6_95$horizon_M, file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_M_K6_ETS_alpha_0.95_CDF.rds"))


# ==============================================================================
# READ RESULTS AND EXPORT TO SHINY APP
# ==============================================================================

# ----------------
# Level of significance = 0.8
# ----------------
hdfpca_int_fore_subnational_err_F_EVR_ETS <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_CDF.rds"))
hdfpca_int_fore_subnational_err_M_EVR_ETS <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_CDF.rds"))

# FIXED: Points to the distinct K6 RDS files instead of duplicating EVR
hdfpca_int_fore_subnational_err_F_K6_ETS  <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_CDF.rds"))
hdfpca_int_fore_subnational_err_M_K6_ETS  <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_CDF.rds"))

# ----------------
# Level of significance = 0.95
# ----------------
hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95 <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_CDF.rds"))
hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95 <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_CDF.rds"))

# FIXED: Points to the distinct K6 RDS files instead of duplicating EVR
hdfpca_int_fore_subnational_err_F_K6_ETS_alpha_0.95  <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_F_K6_ETS_alpha_0.95_CDF.rds"))
hdfpca_int_fore_subnational_err_M_K6_ETS_alpha_0.95  <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_M_K6_ETS_alpha_0.95_CDF.rds"))

# List of arrays to process (47 × 16 × 3)
arrays_list <- list(
  F_EVR_80 = hdfpca_int_fore_subnational_err_F_EVR_ETS,
  M_EVR_80 = hdfpca_int_fore_subnational_err_M_EVR_ETS,
  F_K6_80  = hdfpca_int_fore_subnational_err_F_K6_ETS,
  M_K6_80  = hdfpca_int_fore_subnational_err_M_K6_ETS,
  F_EVR_95 = hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95,
  M_EVR_95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95,
  F_K6_95  = hdfpca_int_fore_subnational_err_F_K6_ETS_alpha_0.95,
  M_K6_95  = hdfpca_int_fore_subnational_err_M_K6_ETS_alpha_0.95
)

# Loop and save with clean Shiny names
for(name in names(arrays_list)) {
  
  arr <- arrays_list[[name]]
  
  # Determine gender
  gender <- ifelse(grepl("^F", name), "female", "male")
  
  # Determine component (EVR or K)
  comp <- ifelse(grepl("EVR", name), "EVR", "K")
  
  # Determine level
  level <- ifelse(grepl("80", name), "80", "95")
  
  # Construct Shiny-friendly file name
  file_name <- paste0("HDFPCA_", gender, "_", comp, "_", level, ".rds")
  
  # Save to Shiny directory
  saveRDS(arr, file = file.path(dir.shiny2, file_name))
}

