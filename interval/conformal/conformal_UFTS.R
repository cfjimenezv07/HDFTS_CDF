###########################
# split conform prediction
###########################

source("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/auxiliary_source/auxiliary_interval.R")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/HLS/"


# fdata: functional data of dimension n by p
# method_ncomp: EVR or K = 6
# horizon: forecast horizon
# fore_method: transformation
# uni_fore_method: forecasting method
# level_sig: level of significance

interval_fore_subnational_cdf_conformal <- function(fdata, method_ncomp, horizon, fore_method, uni_fore_method,
                                                 level_sig)
{
    n_age = ncol(fdata)
    fore_validation = matrix(NA, ncol(fdata), (18 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(18 - horizon))
        {
            fore_validation[,ij] = fore_national_cdf(data_set = fdata[1:(15 + ij),],
                                                     ncomp_method = method_ncomp,
                                                     fh = horizon, fmethod = uni_fore_method)
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(18 - horizon))
        {
            fore_validation[,ij] = as.numeric(clr_fun(fdata = fdata[1:(15 + ij),],
                                                      ncomp_selection = method_ncomp,
                                                      fh = horizon, fore_method = "ETS")$fore_count)
            rm(ij)
        }
    }
    else
    {
        warning("forecasting method must either be CDF or CLR.")
    }
    
    # holdout validation data
    
    holdout_validation_dum = t(matrix(fdata[(16 + horizon):33,], length((16 + horizon):33), ncol(fdata)))
    resi_mat = holdout_validation_dum - fore_validation
    
    quantile_resid <- apply(resi_mat, 1, function(x) quantile(abs(x), probs = level_sig))
    
    fore_val = fore_val_lb = fore_val_ub = matrix(NA, ncol(fdata), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] <- fore_national_cdf(data_set = fdata[1:(32 + ij),], ncomp_method = method_ncomp,
                                               fh = horizon, fmethod = uni_fore_method)
            fore_val_lb[,ij] = fore_val[,ij] - quantile_resid
            fore_val_ub[,ij] = fore_val[,ij] + quantile_resid
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] = as.numeric(clr_fun(fdata = fdata[1:(32 + ij),], ncomp_selection = method_ncomp,
                                               fh = horizon, fore_method = "ETS")$fore_count)
            fore_val_lb[,ij] = fore_val[,ij] - quantile_resid
            fore_val_ub[,ij] = fore_val[,ij] + quantile_resid
            rm(ij)
        }
    }
    
    # holdout testing data
    
    holdout_val_dum = t(matrix(fdata[(33 + horizon):49,], length((33 + horizon):49), ncol(fdata)))
    
    int_err = interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, alpha = (1 - level_sig))
    return(list(int_err = int_err, dimn = ncol(resi_mat)))
}

######
# CPD
######

### CDF

## level_sig = 0.8

# EVR

int_fore_subnational_err_F_EVR_conformal = int_fore_subnational_err_M_EVR_conformal = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
        
        dum = interval_fore_subnational_cdf_conformal(fdata = female_prefecture_dx[[ij]], method_ncomp = "EVR",
                                                   horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                   level_sig = 0.8)
        int_fore_subnational_err_F_EVR_conformal[iw,,ij] = dum$int_err
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf_conformal(fdata = male_prefecture_dx[[ij]], method_ncomp = "EVR",
                                                   horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                   level_sig = 0.8)
        int_fore_subnational_err_M_EVR_conformal[iw,,ij] = dum$int_err
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_EVR_conformal_mean = apply(int_fore_subnational_err_F_EVR_conformal, c(1, 2), mean)
int_fore_subnational_err_M_EVR_conformal_mean = apply(int_fore_subnational_err_M_EVR_conformal, c(1, 2), mean)

colnames(int_fore_subnational_err_F_EVR_conformal_mean) = colnames(int_fore_subnational_err_M_EVR_conformal_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_EVR_conformal_mean) = rownames(int_fore_subnational_err_M_EVR_conformal_mean) = 1:16

# K = 6

int_fore_subnational_err_F_K6_conformal = int_fore_subnational_err_M_K6_conformal = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
        
        dum = interval_fore_subnational_cdf_conformal(fdata = female_prefecture_dx[[ij]], method_ncomp = "provide",
                                                   horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                   level_sig = 0.8)
        int_fore_subnational_err_F_K6_conformal[iw,,ij] = dum$int_err
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf_conformal(fdata = male_prefecture_dx[[ij]], method_ncomp = "provide",
                                                   horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                   level_sig = 0.8)
        int_fore_subnational_err_M_K6_conformal[iw,,ij] = dum$int_err
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_K6_conformal_mean = apply(int_fore_subnational_err_F_K6_conformal, c(1, 2), mean)
int_fore_subnational_err_M_K6_conformal_mean = apply(int_fore_subnational_err_M_K6_conformal, c(1, 2), mean)

colnames(int_fore_subnational_err_F_K6_conformal_mean) = colnames(int_fore_subnational_err_M_K6_conformal_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_K6_conformal_mean) = rownames(int_fore_subnational_err_M_K6_conformal_mean) = 1:16

## level_sig = 0.95

# EVR

int_fore_subnational_err_F_EVR_conformal_alpha_0.95 = int_fore_subnational_err_M_EVR_conformal_alpha_0.95 = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
    
          dum = interval_fore_subnational_cdf_conformal(fdata = female_prefecture_dx[[ij]], method_ncomp = "EVR",
                                                     horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                     level_sig = 0.95)
          int_fore_subnational_err_F_EVR_conformal_alpha_0.95[iw,,ij] = dum$int_err
          rm(dum)
          
          ## (M)
          
          dum = interval_fore_subnational_cdf_conformal(fdata = male_prefecture_dx[[ij]], method_ncomp = "EVR",
                                                     horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                     level_sig = 0.95)
          int_fore_subnational_err_M_EVR_conformal_alpha_0.95[iw,,ij] = dum$int_err
          rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_F_EVR_conformal_alpha_0.95, c(1, 2), mean)
int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_M_EVR_conformal_alpha_0.95, c(1, 2), mean)

colnames(int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = colnames(int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = 1:16

# K = 6

int_fore_subnational_err_F_K6_conformal_alpha_0.95 = int_fore_subnational_err_M_K6_conformal_alpha_0.95 = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
        
        dum = interval_fore_subnational_cdf_conformal(fdata = female_prefecture_dx[[ij]], method_ncomp = "provide",
                                                   horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                   level_sig = 0.95)
        int_fore_subnational_err_F_K6_conformal_alpha_0.95[iw,,ij] = dum$int_err
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf_conformal(fdata = male_prefecture_dx[[ij]], method_ncomp = "provide",
                                                   horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                   level_sig = 0.95)
        int_fore_subnational_err_M_K6_conformal_alpha_0.95[iw,,ij] = dum$int_err
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_F_K6_conformal_alpha_0.95, c(1, 2), mean)
int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_M_K6_conformal_alpha_0.95, c(1, 2), mean)

colnames(int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean) = colnames(int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean) = 1:16

# ===============================
# level_sig = 0.8
# ===============================

# --- EVR ---
saveRDS(int_fore_subnational_err_F_EVR_conformal,
        file = paste0(dir.r, "int_fore_subnational_err_F_EVR_conformal.rds"))
saveRDS(int_fore_subnational_err_M_EVR_conformal,
        file = paste0(dir.r, "int_fore_subnational_err_M_EVR_conformal.rds"))

saveRDS(int_fore_subnational_err_F_EVR_conformal_mean,
        file = paste0(dir.r, "int_fore_subnational_err_F_EVR_conformal_mean.rds"))
saveRDS(int_fore_subnational_err_M_EVR_conformal_mean,
        file = paste0(dir.r, "int_fore_subnational_err_M_EVR_conformal_mean.rds"))


# --- K = 6 ---
saveRDS(int_fore_subnational_err_F_K6_conformal,
        file = paste0(dir.r, "int_fore_subnational_err_F_K6_conformal.rds"))
saveRDS(int_fore_subnational_err_M_K6_conformal,
        file = paste0(dir.r, "int_fore_subnational_err_M_K6_conformal.rds"))

saveRDS(int_fore_subnational_err_F_K6_conformal_mean,
        file = paste0(dir.r, "int_fore_subnational_err_F_K6_conformal_mean.rds"))
saveRDS(int_fore_subnational_err_M_K6_conformal_mean,
        file = paste0(dir.r, "int_fore_subnational_err_M_K6_conformal_mean.rds"))


# ===============================
# level_sig = 0.95
# ===============================

# --- EVR ---
saveRDS(int_fore_subnational_err_F_EVR_conformal_alpha_0.95,
        file = paste0(dir.r, "int_fore_subnational_err_F_EVR_conformal_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_conformal_alpha_0.95,
        file = paste0(dir.r, "int_fore_subnational_err_M_EVR_conformal_alpha_0.95.rds"))

saveRDS(int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean,
        file = paste0(dir.r, "int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean.rds"))
saveRDS(int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean,
        file = paste0(dir.r, "int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean.rds"))


# --- K = 6 ---
saveRDS(int_fore_subnational_err_F_K6_conformal_alpha_0.95,
        file = paste0(dir.r, "int_fore_subnational_err_F_K6_conformal_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_K6_conformal_alpha_0.95,
        file = paste0(dir.r, "int_fore_subnational_err_M_K6_conformal_alpha_0.95.rds"))

saveRDS(int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean,
        file = paste0(dir.r, "int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean.rds"))
saveRDS(int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean,
        file = paste0(dir.r, "int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean.rds"))
