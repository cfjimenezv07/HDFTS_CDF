# fdata_F: female data
# fdata_M: male data
# fore_method: transformation
# horizon: forecast horizon
# way_ncomp: way of selecting number of components
# uni_fore_method: univariate time-series forecasting method
# level_sig: level of significance
# type_PI: pointwise or uniform

interval_fore_subnational_cdf_MFTS <- function(fdata_F, fdata_M, fore_method, horizon, way_ncomp, 
                                            uni_fore_method, level_sig, type_PI)
{
    n_age = ncol(fdata_F)
    forecast_validation_F = forecast_validation_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- fore_national_cdf_MFTS(data_set_F = fdata_F[1:(15+ij),], data_set_M = fdata_M[1:(15+ij),],
                                          fh = horizon, fmethod = uni_fore_method, method_ncomp = way_ncomp)
            forecast_validation_F[,ij] = dum$mfts_fore_F
            forecast_validation_M[,ij] = dum$mfts_fore_M
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- clr_MFTS_fun(fdata_F = fdata_F[1:(15+ij),], fdata_M = fdata_M[1:(15+ij),],
                                ncomp_selection = way_ncomp, fh = horizon, fore_method = "ets")
            forecast_validation_F[,ij] = dum$MFTS_res_fore_F
            forecast_validation_M[,ij] = dum$MFTS_res_fore_M
            rm(ij)
        }
    }
    else
    {
        warning("Forecasting method must either be CDF or CLR.")
    }
    
    # holdout validation data
    
    holdout_validation_F = t(matrix(fdata_F[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_F)))
    holdout_validation_M = t(matrix(fdata_M[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_M)))
    resi_mat_F = holdout_validation_F - forecast_validation_F
    resi_mat_M = holdout_validation_M - forecast_validation_M
    
    # compute standard deviation of residuals
    
    sd_val_F = apply(resi_mat_F, 1, sd)
    sd_val_M = apply(resi_mat_M, 1, sd)
    
    # find the optimal tuning parameter
    
    tune_para_find_val_F_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                      resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_F_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                      resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_F_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                      resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_F_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                      resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_F_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                   resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                   alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_F_6 = optim(par = 1, fn = tune_para_find_function, method = "Nelder-Mead",
                                   resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                   alpha_level = level_sig, PI_type = type_PI)
    
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
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_M_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                      resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_M_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                      resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_M_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                      resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                      alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_M_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                   resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                   alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_M_6 = optim(par = 1, fn = tune_para_find_function, method = "Nelder-Mead",
                                   resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                   alpha_level = level_sig, PI_type = type_PI)
    
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
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- fore_national_cdf_MFTS(data_set_F = fdata_F[1:(31+ij),], data_set_M = fdata_M[1:(31+ij),],
                                          fh = horizon, fmethod = uni_fore_method, method_ncomp = way_ncomp)
            forecast_test_F[,ij] = dum$mfts_fore_F
            forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
            forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
            
            forecast_test_M[,ij] = dum$mfts_fore_M
            forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
            forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
            rm(dum); rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- clr_MFTS_fun(fdata_F = fdata_F[1:(31+ij),], fdata_M = fdata_M[1:(31+ij),],
                                ncomp_selection = way_ncomp, fh = horizon, fore_method = "ets")
            forecast_test_F[,ij] = dum$MFTS_res_fore_F
            forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
            forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
            
            forecast_test_M[,ij] = dum$MFTS_res_fore_M
            forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
            forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
            rm(dum); rm(ij)
        }
    }
    
    # holdout testing data
    
    holdout_val_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_M)))
    
    if(type_PI == "pointwise")
    {
      int_F_err = interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
      int_M_err = interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
    }
    else if(type_PI == "uniform")
    {
      int_F_err = uniform_cpd(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
      int_M_err = uniform_cpd(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
    }
    
    return(list(int_F_err = int_F_err, tune_para_find_F = tune_para_find_F, 
                tune_para_find_F_obj = obj_val_min_F,
                int_M_err = int_M_err, tune_para_find_M = tune_para_find_M,
                tune_para_find_M_obj = obj_val_min_M))
}

###############################
## level of significance = 0.8
###############################

# EVR

MFTS_int_fore_subnational_err_F_EVR_ETS = MFTS_int_fore_subnational_err_M_EVR_ETS = array(NA, dim = c(15, 3, 47))
MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para = MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para_obj = 
MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para = MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para_obj = matrix(NA, 15, 47)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = interval_fore_subnational_cdf_MFTS(fdata_F = female_prefecture_dx[[ij]], 
                                                 fdata_M = male_prefecture_dx[[ij]], 
                                                 fore_method = "CDF", 
                                                 horizon = iw, way_ncomp = "EVR", uni_fore_method = "ets",
                                                 level_sig = 0.8, type_PI = "pointwise")
        MFTS_int_fore_subnational_err_F_EVR_ETS[iw,,ij] = dum$int_F_err
        MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para[iw,ij] = dum$tune_para_find_F
        MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_F_obj
        
        MFTS_int_fore_subnational_err_M_EVR_ETS[iw,,ij] = dum$int_M_err
        MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para[iw,ij] = dum$tune_para_find_M
        MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_M_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

MFTS_int_fore_subnational_err_F_EVR_ETS_mean = apply(MFTS_int_fore_subnational_err_F_EVR_ETS, c(1, 2), mean)
MFTS_int_fore_subnational_err_M_EVR_ETS_mean = apply(MFTS_int_fore_subnational_err_M_EVR_ETS, c(1, 2), mean)

colnames(MFTS_int_fore_subnational_err_F_EVR_ETS_mean) = colnames(MFTS_int_fore_subnational_err_M_EVR_ETS_mean) = c("ECP", "CPD", "score")
rownames(MFTS_int_fore_subnational_err_F_EVR_ETS_mean) = rownames(MFTS_int_fore_subnational_err_M_EVR_ETS_mean) = 1:15

# K = 6

MFTS_int_fore_subnational_err_F_K6_ETS = MFTS_int_fore_subnational_err_M_K6_ETS = array(NA, dim = c(15, 3, 47))
MFTS_int_fore_subnational_err_F_K6_ETS_tune_para = MFTS_int_fore_subnational_err_F_K6_ETS_tune_para_obj = 
MFTS_int_fore_subnational_err_M_K6_ETS_tune_para = MFTS_int_fore_subnational_err_M_K6_ETS_tune_para_obj = matrix(NA, 15, 47)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = interval_fore_subnational_cdf_MFTS(fdata_F = female_prefecture_dx[[ij]], 
                                              fdata_M = male_prefecture_dx[[ij]], 
                                              fore_method = "CDF", 
                                              horizon = iw, way_ncomp = "provide", uni_fore_method = "ets",
                                              level_sig = 0.8, type_PI = "pointwise")
        MFTS_int_fore_subnational_err_F_K6_ETS[iw,,ij] = dum$int_F_err
        MFTS_int_fore_subnational_err_F_K6_ETS_tune_para[iw,ij] = dum$tune_para_find_F
        MFTS_int_fore_subnational_err_F_K6_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_F_obj
        
        MFTS_int_fore_subnational_err_M_K6_ETS[iw,,ij] = dum$int_M_err
        MFTS_int_fore_subnational_err_M_K6_ETS_tune_para[iw,ij] = dum$tune_para_find_M
        MFTS_int_fore_subnational_err_M_K6_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_M_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

MFTS_int_fore_subnational_err_F_K6_ETS_mean = apply(MFTS_int_fore_subnational_err_F_K6_ETS, c(1, 2), mean)
MFTS_int_fore_subnational_err_M_K6_ETS_mean = apply(MFTS_int_fore_subnational_err_M_K6_ETS, c(1, 2), mean)

colnames(MFTS_int_fore_subnational_err_F_K6_ETS_mean) = colnames(MFTS_int_fore_subnational_err_M_K6_ETS_mean) = c("ECP", "CPD", "score")
rownames(MFTS_int_fore_subnational_err_F_K6_ETS_mean) = rownames(MFTS_int_fore_subnational_err_M_K6_ETS_mean) = 1:15

################################
## level of significance = 0.95
################################

# EVR

MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95 = MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95 = array(NA, dim = c(15, 3, 47))
MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95 = MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95 = 
MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95 = MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95 = matrix(NA, 15, 47)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = interval_fore_subnational_cdf_MFTS(fdata_F = female_prefecture_dx[[ij]], 
                                                 fdata_M = male_prefecture_dx[[ij]], 
                                                 fore_method = "CDF", 
                                                 horizon = iw, way_ncomp = "EVR", uni_fore_method = "ets",
                                                 level_sig = 0.95, type_PI = "pointwise")
        MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95[iw,,ij] = dum$int_F_err
        MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find_F
        MFTS_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_F_obj
        
        MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95[iw,,ij] = dum$int_M_err
        MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find_M
        MFTS_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_M_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean = apply(MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95, c(1, 2), mean)
MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean = apply(MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95, c(1, 2), mean)

colnames(MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean) = colnames(MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean) = rownames(MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean) = 1:15

# K = 6

MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95 = MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95 = array(NA, dim = c(15, 3, 47))
MFTS_int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95 = MFTS_int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95 = 
MFTS_int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95 = MFTS_int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95 = matrix(NA, 15, 47)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = interval_fore_subnational_cdf_MFTS(fdata_F = female_prefecture_dx[[ij]], 
                                                 fdata_M = male_prefecture_dx[[ij]], 
                                                 fore_method = "CDF", 
                                                 horizon = iw, way_ncomp = "provide", uni_fore_method = "ets",
                                                 level_sig = 0.95, type_PI = "pointwise")
        MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95[iw,,ij] = dum$int_F_err
        MFTS_int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find_F
        MFTS_int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_F_obj
        
        MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95[iw,,ij] = dum$int_M_err
        MFTS_int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find_M
        MFTS_int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_M_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean = apply(MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95, c(1, 2), mean)
MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean = apply(MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95, c(1, 2), mean)

colnames(MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean) = colnames(MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean) = rownames(MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean) = 1:15

