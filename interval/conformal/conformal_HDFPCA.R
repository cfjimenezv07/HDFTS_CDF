###########################################
# conformal pointwise prediction intervals
###########################################

# fdata_F_data: female data series
# fdata_M_data: male data series
# fore_method: forecasting method
# horizon: forecast horizon
# level_sig: level of significance

hdfpca_interval_fore_subnational_cdf_conformal <- function(fdata_F_data, fdata_M_data, fore_method, horizon, level_sig)
{
    n_age = ncol(fdata_F_data)
    if(fore_method == "CDF")
    {
        dum <- hdfpca_fun_fore(fdata_F = fdata_F_data[1:32,], fdata_M = fdata_M_data[1:32,],
                               horizon = horizon, first_order = 6, second_order = 2, transformation = "CDF")
        forecast_validation_F = dum$forecast_val_F
        forecast_validation_M = dum$forecast_val_M
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- clr_MFTS_fun(fdata_F = fdata_F_data[1:(15+ij),], fdata_M = fdata_M_data[1:(15+ij),],
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
    
    holdout_validation_F = t(matrix(fdata_F_data[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_F_data)))
    holdout_validation_M = t(matrix(fdata_M_data[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_M_data)))
    resi_mat_F = holdout_validation_F - forecast_validation_F
    resi_mat_M = holdout_validation_M - forecast_validation_M
    
    quantile_resid_F <- apply(resi_mat_F, 1, function(x) quantile(abs(x), probs = level_sig))
    quantile_resid_M <- apply(resi_mat_M, 1, function(x) quantile(abs(x), probs = level_sig))
    
    if(fore_method == "CDF")
    {
        dum <- hdfpca_fun_fore(fdata_F = fdata_F_data, fdata_M = fdata_M_data,
                          horizon = horizon, first_order = 6, second_order = 2, transformation = "CDF")
        forecast_test_F = dum$forecast_val_F
        forecast_test_F_lb = forecast_test_F - quantile_resid_F
        forecast_test_F_ub = forecast_test_F + quantile_resid_F
            
        forecast_test_M = dum$forecast_val_M
        forecast_test_M_lb = forecast_test_M - quantile_resid_M
        forecast_test_M_ub = forecast_test_M + quantile_resid_M
        rm(dum)
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- clr_MFTS_fun(fdata_F = fdata_F_data[1:(31+ij),], fdata_M = fdata_M_data[1:(31+ij),],
                                ncomp_selection = way_ncomp, fh = horizon, fore_method = "ets")
            forecast_test_F[,ij] = dum$MFTS_res_fore_F
            forecast_test_F_lb[,ij] = forecast_test_F[,ij] - quantile_resid_F
            forecast_test_F_ub[,ij] = forecast_test_F[,ij] + quantile_resid_F
            
            forecast_test_M[,ij] = dum$MFTS_res_fore_M
            forecast_test_M_lb[,ij] = forecast_test_M[,ij] - quantile_resid_M
            forecast_test_M_ub[,ij] = forecast_test_M[,ij] + quantile_resid_M
            rm(dum); rm(ij)
        }
    }
    
    # holdout testing data
    
    holdout_val_F = t(matrix(fdata_F_data[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F_data)))
    holdout_val_M = t(matrix(fdata_M_data[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_M_data)))
    
    int_F_err = interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
    int_M_err = interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
    return(list(int_F_err = int_F_err, int_M_err = int_M_err))
}

###############################
## level of significance = 0.8
###############################

hdfpca_int_fore_subnational_err_F_EVR_conformal = hdfpca_int_fore_subnational_err_M_EVR_conformal = array(NA, dim = c(15, 3, 47))
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = hdfpca_interval_fore_subnational_cdf_conformal(fdata_F_data = female_prefecture_dx[[ij]], 
                                                             fdata_M_data = male_prefecture_dx[[ij]], 
                                                             fore_method = "CDF", 
                                                             horizon = iw, level_sig = 0.8)
        hdfpca_int_fore_subnational_err_F_EVR_conformal[iw,,ij] = dum$int_F_err
        hdfpca_int_fore_subnational_err_M_EVR_conformal[iw,,ij] = dum$int_M_err
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

hdfpca_int_fore_subnational_err_F_EVR_conformal_mean = apply(hdfpca_int_fore_subnational_err_F_EVR_conformal, c(1, 2), mean)
hdfpca_int_fore_subnational_err_M_EVR_conformal_mean = apply(hdfpca_int_fore_subnational_err_M_EVR_conformal, c(1, 2), mean)

colnames(hdfpca_int_fore_subnational_err_F_EVR_conformal_mean) = colnames(hdfpca_int_fore_subnational_err_M_EVR_conformal_mean) = c("ECP", "CPD", "score")
rownames(hdfpca_int_fore_subnational_err_F_EVR_conformal_mean) = rownames(hdfpca_int_fore_subnational_err_M_EVR_conformal_mean) = 1:15

################################
## level of significance = 0.95
################################

hdfpca_int_fore_subnational_err_F_EVR_conformal_alpha_0.95 =
hdfpca_int_fore_subnational_err_M_EVR_conformal_alpha_0.95 = array(NA, dim = c(15, 3, 47))
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = hdfpca_interval_fore_subnational_cdf_conformal(fdata_F_data = female_prefecture_dx[[ij]], 
                                                             fdata_M_data = male_prefecture_dx[[ij]], 
                                                             fore_method = "CDF", 
                                                             horizon = iw, level_sig = 0.95)
        hdfpca_int_fore_subnational_err_F_EVR_conformal_alpha_0.95[iw,,ij] = dum$int_F_err
        hdfpca_int_fore_subnational_err_M_EVR_conformal_alpha_0.95[iw,,ij] = dum$int_M_err
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

hdfpca_int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean = apply(hdfpca_int_fore_subnational_err_F_EVR_conformal_alpha_0.95, c(1, 2), mean)
hdfpca_int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean = apply(hdfpca_int_fore_subnational_err_M_EVR_conformal_alpha_0.95, c(1, 2), mean)

colnames(hdfpca_int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = colnames(hdfpca_int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(hdfpca_int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = rownames(hdfpca_int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = 1:15

