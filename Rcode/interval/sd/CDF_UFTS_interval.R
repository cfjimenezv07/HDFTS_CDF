
source("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/auxiliary_source/auxiliary_interval.R")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/HLS/"


#########################
# train_set: 1:16
# validation_set: 17:32
# test_set: 33:49
#########################

# fdata: functional data n by p
# method_ncomp: EVR or K = 6
# horizon: forecast horizon
# fore_method: transformation
# uni_fore_method: forecasting method
# level_sig: level of significance
# type_PI: prediction interval or prediction band

interval_fore_subnational_cdf <- function(fdata, method_ncomp, horizon, fore_method, uni_fore_method,
                                       level_sig, type_PI)
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
           print(ij) 
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
    
    # compute standard deviation of residuals
    
    sd_val_input = apply(resi_mat, 1, sd)
    
    # find the optimal tuning parameter
    
    tune_para_find_val_1 = optimise(f = tune_para_find_function, interval = c(0, 1), 
                                    resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                    alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_2 = optimise(f = tune_para_find_function, interval = c(0, 5), 
                                    resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                    alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_3 = optimise(f = tune_para_find_function, interval = c(0, 10), 
                                    resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                    alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_4 = optimise(f = tune_para_find_function, interval = c(0, 20), 
                                    resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                    alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B", 
                                 resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                 alpha_level = level_sig, PI_type = type_PI)
    
    tune_para_find_val_6 = optim(par = 1, fn = tune_para_find_function, method = "Nelder-Mead",
                                 resi_mat = resi_mat, sd_val_input = sd_val_input, 
                                 alpha_level = level_sig, PI_type = type_PI)
    
    obj_val = c(tune_para_find_val_1$objective, 
                tune_para_find_val_2$objective, 
                tune_para_find_val_3$objective,
                tune_para_find_val_4$objective,
                tune_para_find_val_5$value,
                tune_para_find_val_6$value)
    obj_val_min = min(obj_val)
    
    tune_para_find = c(tune_para_find_val_1$minimum, 
                       tune_para_find_val_2$minimum,
                       tune_para_find_val_3$minimum,
                       tune_para_find_val_4$minimum,
                       tune_para_find_val_5$par,
                       tune_para_find_val_6$par)[which.min(obj_val)]
    
    fore_val = fore_val_lb = fore_val_ub = matrix(NA, ncol(fdata), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] <- fore_national_cdf(data_set = fdata[1:(32 + ij),], ncomp_method = method_ncomp,
                                               fh = horizon, fmethod = uni_fore_method)
            fore_val_lb[,ij] = fore_val[,ij] - tune_para_find * sd_val_input
            fore_val_ub[,ij] = fore_val[,ij] + tune_para_find * sd_val_input
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] = as.numeric(clr_fun(fdata = fdata[1:(32 + ij),], ncomp_selection = method_ncomp,
                                               fh = horizon, fore_method = "ETS")$fore_count)
            fore_val_lb[,ij] = fore_val[,ij] - tune_para_find * sd_val_input
            fore_val_ub[,ij] = fore_val[,ij] + tune_para_find * sd_val_input
            rm(ij)
        }
    }
    
    # holdout testing data
    
    holdout_val_dum = t(matrix(fdata[(33 + horizon):49,], length((33 + horizon):49), ncol(fdata)))
    
    if(type_PI == "pointwise")
    {
        int_err = interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, 
                                 alpha = (1 - level_sig))
    }
    else if(type_PI == "uniform")
    {
        int_err = uniform_cpd(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub,
                              alpha = (1 - level_sig))
    }
    return(list(int_err = int_err, tune_para_find = tune_para_find, tune_para_find_obj = obj_val_min))
}

########
### CDF
########

## level_sig = 0.8

# EVR

int_fore_subnational_err_F_EVR_ETS = int_fore_subnational_err_M_EVR_ETS = array(NA, dim = c(16, 3, 47))
int_fore_subnational_err_F_EVR_ETS_tune_para = int_fore_subnational_err_F_EVR_ETS_tune_para_obj = 
int_fore_subnational_err_M_EVR_ETS_tune_para = int_fore_subnational_err_M_EVR_ETS_tune_para_obj = matrix(NA, 16, 47)
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
    
        dum = interval_fore_subnational_cdf(fdata = female_prefecture_dx[[ij]], method_ncomp = "EVR",
                                         horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                         level_sig = 0.8, type_PI = "pointwise")
        int_fore_subnational_err_F_EVR_ETS[iw,,ij] = dum$int_err
        int_fore_subnational_err_F_EVR_ETS_tune_para[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_F_EVR_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_obj
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf(fdata = male_prefecture_dx[[ij]], method_ncomp = "EVR",
                                         horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                         level_sig = 0.8, type_PI = "pointwise")
        int_fore_subnational_err_M_EVR_ETS[iw,,ij] = dum$int_err
        int_fore_subnational_err_M_EVR_ETS_tune_para[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_M_EVR_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_EVR_ETS_mean = apply(int_fore_subnational_err_F_EVR_ETS, c(1, 2), mean)
int_fore_subnational_err_M_EVR_ETS_mean = apply(int_fore_subnational_err_M_EVR_ETS, c(1, 2), mean)

colnames(int_fore_subnational_err_F_EVR_ETS_mean) = colnames(int_fore_subnational_err_M_EVR_ETS_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_EVR_ETS_mean) = rownames(int_fore_subnational_err_M_EVR_ETS_mean) = 1:16

# K = 6

int_fore_subnational_err_F_K6_ETS = int_fore_subnational_err_M_K6_ETS = array(NA, dim = c(16, 3, 47))
int_fore_subnational_err_F_K6_ETS_tune_para = int_fore_subnational_err_F_K6_ETS_tune_para_obj = 
int_fore_subnational_err_M_K6_ETS_tune_para = int_fore_subnational_err_M_K6_ETS_tune_para_obj = matrix(NA, 16, 47)
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
        
        dum = interval_fore_subnational_cdf(fdata = female_prefecture_dx[[ij]], method_ncomp = "provide",
                                         horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                         level_sig = 0.8, type_PI = "pointwise")
        int_fore_subnational_err_F_K6_ETS[iw,,ij] = dum$int_err
        int_fore_subnational_err_F_K6_ETS_tune_para[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_F_K6_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_obj
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf(fdata = male_prefecture_dx[[ij]], method_ncomp = "provide",
                                         horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                         level_sig = 0.8, type_PI = "pointwise")
        int_fore_subnational_err_M_K6_ETS[iw,,ij] = dum$int_err
        int_fore_subnational_err_M_K6_ETS_tune_para[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_M_K6_ETS_tune_para_obj[iw,ij] = dum$tune_para_find_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_K6_ETS_mean = apply(int_fore_subnational_err_F_K6_ETS, c(1, 2), mean)
int_fore_subnational_err_M_K6_ETS_mean = apply(int_fore_subnational_err_M_K6_ETS, c(1, 2), mean)

colnames(int_fore_subnational_err_F_K6_ETS_mean) = colnames(int_fore_subnational_err_M_K6_ETS_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_K6_ETS_mean) = rownames(int_fore_subnational_err_M_K6_ETS_mean) = 1:16

## level_sig = 0.95

# EVR

int_fore_subnational_err_F_EVR_ETS_alpha_0.95 = int_fore_subnational_err_M_EVR_ETS_alpha_0.95 = array(NA, dim = c(16, 3, 47))
int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95 = int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95 = 
int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95 = int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95 = matrix(NA, 16, 47)
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
        
        dum = interval_fore_subnational_cdf(fdata = female_prefecture_dx[[ij]], method_ncomp = "EVR",
                                            horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                            level_sig = 0.95, type_PI = "pointwise")
        int_fore_subnational_err_F_EVR_ETS_alpha_0.95[iw,,ij] = dum$int_err
        int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_obj
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf(fdata = male_prefecture_dx[[ij]], method_ncomp = "EVR",
                                            horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                            level_sig = 0.95, type_PI = "pointwise")
        int_fore_subnational_err_M_EVR_ETS_alpha_0.95[iw,,ij] = dum$int_err
        int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean = apply(int_fore_subnational_err_F_EVR_ETS_alpha_0.95, c(1, 2), mean)
int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean = apply(int_fore_subnational_err_M_EVR_ETS_alpha_0.95, c(1, 2), mean)

colnames(int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean) = colnames(int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean) = 1:16

# K = 6

int_fore_subnational_err_F_K6_ETS_alpha_0.95 = int_fore_subnational_err_M_K6_ETS_alpha_0.95 = array(NA, dim = c(16, 3, 47))
int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95 = int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95 = 
int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95 = int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95 = matrix(NA, 16, 47)
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## (F)
        
        dum = interval_fore_subnational_cdf(fdata = female_prefecture_dx[[ij]], method_ncomp = "provide",
                                         horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                         level_sig = 0.95, type_PI = "pointwise")
        int_fore_subnational_err_F_K6_ETS_alpha_0.95[iw,,ij] = dum$int_err
        int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_obj
        rm(dum)
        
        ## (M)
        
        dum = interval_fore_subnational_cdf(fdata = male_prefecture_dx[[ij]], method_ncomp = "provide",
                                         horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                         level_sig = 0.95, type_PI = "pointwise")
        int_fore_subnational_err_M_K6_ETS_alpha_0.95[iw,,ij] = dum$int_err
        int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95[iw,ij] = dum$tune_para_find
        int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95[iw,ij] = dum$tune_para_find_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean = apply(int_fore_subnational_err_F_K6_ETS_alpha_0.95, c(1, 2), mean)
int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean = apply(int_fore_subnational_err_M_K6_ETS_alpha_0.95, c(1, 2), mean)

colnames(int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean) = colnames(int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean) = 1:16

# ----------------------------
# SAVE ALL RESULTS
# ----------------------------

# EVR (level_sig = 0.8)
saveRDS(int_fore_subnational_err_F_EVR_ETS, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_tune_para.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_tune_para.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para_obj, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_tune_para_obj.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para_obj, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_tune_para_obj.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_mean, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_mean.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_mean, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_mean.rds"))

# K = 6 (level_sig = 0.8)
saveRDS(int_fore_subnational_err_F_K6_ETS, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS.rds"))
saveRDS(int_fore_subnational_err_F_K6_ETS_tune_para, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_tune_para.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_tune_para, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_tune_para.rds"))
saveRDS(int_fore_subnational_err_F_K6_ETS_tune_para_obj, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_tune_para_obj.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_tune_para_obj, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_tune_para_obj.rds"))
saveRDS(int_fore_subnational_err_F_K6_ETS_mean, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_mean.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_mean, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_mean.rds"))

# EVR (level_sig = 0.95)
saveRDS(int_fore_subnational_err_F_EVR_ETS_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean, file = paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean.rds"))
saveRDS(int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean, file = paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean.rds"))

# K = 6 (level_sig = 0.95)
saveRDS(int_fore_subnational_err_F_K6_ETS_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_tune_para_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_tune_para_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_tune_para_obj_alpha_0.95.rds"))
saveRDS(int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean, file = paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean.rds"))
saveRDS(int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean, file = paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean.rds"))



#save for the shiny app them
dir.shiny2 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Shiny_App_CDF/datasets_shiny_app2/IFE/Japan/"


# EVR (level_sig = 0.8)
int_fore_subnational_err_F_EVR_ETS <- readRDS(paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS.rds"))
int_fore_subnational_err_M_EVR_ETS <- readRDS(paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS.rds"))

# K = 6 (level_sig = 0.8)
int_fore_subnational_err_F_K6_ETS <- readRDS(paste0(dir.r, "int_fore_subnational_err_F_K6_ETS.rds"))
int_fore_subnational_err_M_K6_ETS <- readRDS(paste0(dir.r, "int_fore_subnational_err_M_K6_ETS.rds"))

# EVR (level_sig = 0.95)
int_fore_subnational_err_F_EVR_ETS_alpha_0.95 <- readRDS(paste0(dir.r, "int_fore_subnational_err_F_EVR_ETS_alpha_0.95.rds"))
int_fore_subnational_err_M_EVR_ETS_alpha_0.95 <- readRDS(paste0(dir.r, "int_fore_subnational_err_M_EVR_ETS_alpha_0.95.rds"))

# K = 6 (level_sig = 0.95)
int_fore_subnational_err_F_K6_ETS_alpha_0.95 <- readRDS(paste0(dir.r, "int_fore_subnational_err_F_K6_ETS_alpha_0.95.rds"))
int_fore_subnational_err_M_K6_ETS_alpha_0.95 <- readRDS(paste0(dir.r, "int_fore_subnational_err_M_K6_ETS_alpha_0.95.rds"))

# Directory to save
dir.shiny2 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Shiny_App_CDF/datasets_shiny_app2/IFE/Japan/"

# List of arrays to process
arrays_list <- list(
  F_EVR_80      = int_fore_subnational_err_F_EVR_ETS,
  M_EVR_80      = int_fore_subnational_err_M_EVR_ETS,
  F_K6_80       = int_fore_subnational_err_F_K6_ETS,
  M_K6_80       = int_fore_subnational_err_M_K6_ETS,
  F_EVR_95      = int_fore_subnational_err_F_EVR_ETS_alpha_0.95,
  M_EVR_95      = int_fore_subnational_err_M_EVR_ETS_alpha_0.95,
  F_K6_95       = int_fore_subnational_err_F_K6_ETS_alpha_0.95,
  M_K6_95       = int_fore_subnational_err_M_K6_ETS_alpha_0.95
)

# Loop and save
for(name in names(arrays_list)) {
  
  arr <- arrays_list[[name]]
  
  # Reorder dimensions: 16 × 3 × 47 -> 47 × 16 × 3
  arr_new <- aperm(arr, c(3, 1, 2))
  
  # Extract gender, component, and level
  gender <- ifelse(grepl("^F", name), "female", "male")
  comp   <- ifelse(grepl("EVR", name), "EVR", "K")
  level  <- ifelse(grepl("80", name), "80", "95")
  
  # Construct file name
  file_name <- paste0("UFTS_", gender, "_", comp, "_", level, ".rds")
  
  # Save
  saveRDS(arr_new, file = file.path(dir.shiny2, file_name))
}

