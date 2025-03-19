###########################
# split conform prediction
###########################

##################
# Dataset entries
##################

year = 1975:2022
n_year = length(year)
age = 0:110
n_age = length(age)
n_prefectures=47

# Row partition
part_list = list()
for(ik in 1:n_prefectures)
{
  part_list[[ik]] = (n_year*ik-(n_year-1)):(n_year*ik)
}

#Column partition
n_populations=2
part_list_c = list()
for(ik in 1:n_populations)
{
  part_list_c[[ik]] = (n_age*ik-(n_age-1)):(n_age*ik)
}

######################
# model HDFTS of CDFs
######################
# read Japanese Subnational Human Mortality Data

state = c("Hokkaido", 
          "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
          "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
          "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
          "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

# change file directory dir.d

female_prefecture_qx = male_prefecture_qx = total_prefecture_qx = list()
for(ik in 1:length(state))
{
  female_prefecture_qx[[ik]] = t(matrix(read.table(paste(dir.d,paste("female_prefecture_", ik, sep = ""), ".txt", sep = ""), header = TRUE)$qx, n_age, n_year))
  male_prefecture_qx[[ik]]   = t(matrix(read.table(paste(dir.d,paste("male_prefecture_", ik, sep = ""), ".txt", sep = ""),   header = TRUE)$qx, n_age, n_year))
  total_prefecture_qx[[ik]]  = t(matrix(read.table(paste(dir.d,paste("total_prefecture_", ik, sep = ""), ".txt", sep = ""),  header = TRUE)$qx, n_age, n_year))
  print(ik); rm(ik)
}

female_prefecture_dx = male_prefecture_dx = total_prefecture_dx = list()
for(iw in 1:length(state))
{
  female_prefecture_dum = male_prefecture_dum = total_prefecture_dum = matrix(NA, n_year, n_age)
  for(ij in 1:n_year)
  {
    # set radix (normalising to 1)
    start_pop_female = start_pop_male = start_pop_total = 10^5
    for(ik in 1:n_age)
    {
      female_prefecture_dum[ij,ik] = (female_prefecture_qx[[iw]])[ij,ik] * start_pop_female
      start_pop_female = start_pop_female - female_prefecture_dum[ij,ik]
      
      male_prefecture_dum[ij,ik] = (male_prefecture_qx[[iw]])[ij,ik] * start_pop_male
      start_pop_male = start_pop_male - male_prefecture_dum[ij,ik]
      
      total_prefecture_dum[ij,ik] = (total_prefecture_qx[[iw]])[ij,ik] * start_pop_total
      start_pop_total = start_pop_total - total_prefecture_dum[ij,ik]
    }
  }
  female_prefecture_dx[[iw]] = t(female_prefecture_dum)
  male_prefecture_dx[[iw]]   = t(male_prefecture_dum)
  total_prefecture_dx[[iw]]  = t(total_prefecture_dum)
  rm(female_prefecture_dum); rm(male_prefecture_dum); rm(total_prefecture_dum)
  print(iw); rm(iw)
}

All_Japan_female_qx <- female_prefecture_dx
All_Japan_male_qx   <- male_prefecture_dx
n_states <- length(state)

############################################
# Functional mean ANOVA approach (FM-ANOVA)
############################################

new_age <- 0:109
# This function computes the functional mean ANOVA decomposition based on means
FANOVA_means <- hdftsa::FANOVA(data_pop1 = t(all_unconstrained_male), data_pop2 = t(all_unconstrained_female), year, new_age, n_states, n_populations)

#This function computes the functional residuals after removing the deterministic components
# obtained from the FANOVA function.
Residuals_means <- hdftsa::Two_way_Residuals_means(data_pop1 = t(all_unconstrained_male), data_pop2=t(all_unconstrained_female), year, new_age, n_states, n_populations)

Res1_means = Residuals_means$residuals1_mean
Res2_means = Residuals_means$residuals2_mean
Residuals_mean <- cbind(Res1_means,Res2_means)

# Reconstructed data
RR <- Residuals_means$rd #Matrix with the original data reconstructed from the FMP decomposition

# It's the proof of the reconstruction of the residuals.
Residuals_means$R #The result should be a vector with two entries TRUE, TRUE.

#Indicating that after adding both deterministic and time-varying components the FTS are recovered.
Fixed_part_means<-Residuals_means$Fixed_comp_mean # deterministic components to be added up after forecasting
nn_age <- length(new_age)
Fixed_part_means_1 <- Fixed_part_means[,1:nn_age]
Fixed_part_means_2 <- Fixed_part_means[,(nn_age+1):(2*nn_age)]

######################################################################################################
# point forecast (FANOVA)
######################################################################################################
# This function takes fixed and residuals matrices of dimensions n_age by n_year for each prefecture.
# and do point forecast for the residuals later add the Fixed components and transformed them back to the original space.
fore_FANOVA_cdf <- function(Fixed, Residuals, ncomp_method, fh, fmethod,est_method )
{
  
  fmethod="ets"
  data_cumsum_logit_fore =  Pref_forecast_curves(fixed_com = Fixed,
                                                 Residuals_f = Residuals,
                                                 est_method = est_method,
                                                 fh = fh, B = 1000,
                                                 prediction_method=fmethod,select_K=ncomp_method, K=6)$med_polish_curve_forecast
  
  # h-step-ahead mean forecast
  
  data_cumsum_logit_fore_add = c(invlogit(data_cumsum_logit_fore[,fh]), 1)
  data_cumsum_logit_fore_add_diff = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
  
  return(data_cumsum_logit_fore_add_diff * 10^5)
}
###############################
# prediction interval (FANOVA)
###############################

# This function requires inputs:
# Fixed: a n_age by n_year matrix coming from the FANOVA decomposition for a particular prefecture
# Residuals: a n_age by n_year matrix coming from the FANOVA decomposition for a particular prefecture
# Holdout: for a particular gender. A matrix of dimension n_age by n_year
# method_ncomp: EVR or K = 6
# est_method: cov or long-run covariance
# horizon: forecast horizon
# uni_fore_method: univariate time-series forecasting method
# level_sig: level of significance

interval_fore_FANOVA_cdf_conformal <- function(Fixed, Residuals, Holdout, method_ncomp, est_method, horizon, fore_method = "CDF",
                                               uni_fore_method, level_sig)
{
  fdata = t(Holdout)
  n_age = ncol(fdata)
  n_age = ncol(fdata)
  fore_validation = matrix(NA, ncol(fdata), (17 - horizon))
  if(fore_method == "CDF")
  {
    for(ij in 1:(17 - horizon))
    {
      fore_validation[,ij] = fore_FANOVA_cdf(Fixed = Fixed[1:(15 + ij),], Residuals = Residuals[1:(15 + ij),],
                                             ncomp_method = method_ncomp, fh = horizon, fmethod = uni_fore_method, est_method = est_method)
      rm(ij)
    }
  }
  else if(fore_method == "CLR")
  {
    for(ij in 1:(17 - horizon))
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
  
  holdout_validation_dum = t(matrix(fdata[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata)))
  resi_mat = holdout_validation_dum - fore_validation
  
  quantile_resid <- apply(resi_mat, 1, function(x) quantile(abs(x), probs = level_sig))
  
  fore_val = fore_val_lb = fore_val_ub = matrix(NA, ncol(fdata), (17 - horizon))
  if(fore_method == "CDF")
  {
    for(ij in 1:(17 - horizon))
    {
      fore_val[,ij] <- fore_FANOVA_cdf(Fixed=Fixed[1:(31 + ij),], Residuals = Residuals[1:(31 + ij),],
                                       ncomp_method=method_ncomp, fh=horizon, fmethod=uni_fore_method,est_method = est_method)
      fore_val_lb[,ij] = fore_val[,ij] - quantile_resid
      fore_val_ub[,ij] = fore_val[,ij] + quantile_resid
      rm(ij)
    }
  }
  else if(fore_method == "CLR")
  {
    for(ij in 1:(17 - horizon))
    {
      fore_val[,ij] = as.numeric(clr_fun(fdata = fdata[1:(31 + ij),], ncomp_selection = method_ncomp,
                                         fh = horizon, fore_method = "ETS")$fore_count)
      fore_val_lb[,ij] = fore_val[,ij] - quantile_resid
      fore_val_ub[,ij] = fore_val[,ij] + quantile_resid
      rm(ij)
    }
  }
  
  # holdout testing data
  
  holdout_val_dum = t(matrix(fdata[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata)))
  
  int_err = interval_score(holdout = holdout_val_dum, lb = fore_val_lb, ub = fore_val_ub, alpha = (1 - level_sig))
  return(list(int_err = int_err, dimn = ncol(resi_mat)))
}

######
# CPD
######

### CDF

## level_sig = 0.8

# EVR

int_fore_subnational_err_F_EVR_conformal = int_fore_subnational_err_M_EVR_conformal = array(NA, dim = c(15, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:15)
  {
    ## (F)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=female_prefecture_fixed_means[[ij]],Residuals=female_prefecture_res_means[[ij]],
                                   Holdout=female_prefecture_dx[[ij]], method_ncomp="EVR", est_method="cov",
                                   horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.8)
    int_fore_subnational_err_F_EVR_conformal[iw,,ij] = dum$int_err
    rm(dum)
    
    ## (M)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=male_prefecture_fixed_means[[ij]],Residuals=male_prefecture_res_means[[ij]],
                                             Holdout=male_prefecture_dx[[ij]], method_ncomp="EVR", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.8)
    int_fore_subnational_err_M_EVR_conformal[iw,,ij] = dum$int_err
    rm(dum); print(iw); rm(iw)
  }
  print(ij); rm(ij)
}

int_fore_subnational_err_F_EVR_conformal_mean = apply(int_fore_subnational_err_F_EVR_conformal, c(1, 2), mean)
int_fore_subnational_err_M_EVR_conformal_mean = apply(int_fore_subnational_err_M_EVR_conformal, c(1, 2), mean)

colnames(int_fore_subnational_err_F_EVR_conformal_mean) = colnames(int_fore_subnational_err_M_EVR_conformal_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_EVR_conformal_mean) = rownames(int_fore_subnational_err_M_EVR_conformal_mean) = 1:15

# K = 6

int_fore_subnational_err_F_K6_conformal = int_fore_subnational_err_M_K6_conformal = array(NA, dim = c(15, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:15)
  {
    ## (F)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=female_prefecture_fixed_means[[ij]],Residuals=female_prefecture_res_means[[ij]],
                                             Holdout=female_prefecture_dx[[ij]], method_ncomp="Fixed", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.8)
    int_fore_subnational_err_F_K6_conformal[iw,,ij] = dum$int_err
    rm(dum)
    
    ## (M)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=male_prefecture_fixed_means[[ij]],Residuals=male_prefecture_res_means[[ij]],
                                             Holdout=male_prefecture_dx[[ij]], method_ncomp="Fixed", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.8)
    int_fore_subnational_err_M_K6_conformal[iw,,ij] = dum$int_err
    rm(dum); print(iw); rm(iw)
  }
  print(ij); rm(ij)
}

int_fore_subnational_err_F_K6_conformal_mean = apply(int_fore_subnational_err_F_K6_conformal, c(1, 2), mean)
int_fore_subnational_err_M_K6_conformal_mean = apply(int_fore_subnational_err_M_K6_conformal, c(1, 2), mean)

colnames(int_fore_subnational_err_F_K6_conformal_mean) = colnames(int_fore_subnational_err_M_K6_conformal_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_K6_conformal_mean) = rownames(int_fore_subnational_err_M_K6_conformal_mean) = 1:15

####################
## level_sig = 0.95
####################

# EVR

int_fore_subnational_err_F_EVR_conformal_alpha_0.95 = int_fore_subnational_err_M_EVR_conformal_alpha_0.95 = array(NA, dim = c(15, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:15)
  {
    ## (F)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=female_prefecture_fixed_means[[ij]],Residuals=female_prefecture_res_means[[ij]],
                                             Holdout=female_prefecture_dx[[ij]], method_ncomp="EVR", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.95)
    int_fore_subnational_err_F_EVR_conformal_alpha_0.95[iw,,ij] = dum$int_err
    rm(dum)
    
    ## (M)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=male_prefecture_fixed_means[[ij]],Residuals=male_prefecture_res_means[[ij]],
                                             Holdout=male_prefecture_dx[[ij]], method_ncomp="EVR", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.95)
    int_fore_subnational_err_M_EVR_conformal_alpha_0.95[iw,,ij] = dum$int_err
    rm(dum); print(iw); rm(iw)
  }
  print(ij); rm(ij)
}

int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_F_EVR_conformal_alpha_0.95, c(1, 2), mean)
int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_M_EVR_conformal_alpha_0.95, c(1, 2), mean)

colnames(int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = colnames(int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = 1:15

# K = 6

int_fore_subnational_err_F_K6_conformal_alpha_0.95 = int_fore_subnational_err_M_K6_conformal_alpha_0.95 = array(NA, dim = c(15, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:15)
  {
    ## (F)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=female_prefecture_fixed_means[[ij]],Residuals=female_prefecture_res_means[[ij]],
                                             Holdout=female_prefecture_dx[[ij]], method_ncomp="Fixed", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.95)
    int_fore_subnational_err_F_K6_conformal_alpha_0.95[iw,,ij] = dum$int_err
    rm(dum)
    
    ## (M)
    
    dum = interval_fore_FANOVA_cdf_conformal(Fixed=male_prefecture_fixed_means[[ij]],Residuals=male_prefecture_res_means[[ij]],
                                             Holdout=male_prefecture_dx[[ij]], method_ncomp="Fixed", est_method="cov",
                                             horizon=iw, fore_method="CDF", uni_fore_method="ets", level_sig=0.95)
    int_fore_subnational_err_M_K6_conformal_alpha_0.95[iw,,ij] = dum$int_err
    rm(dum); print(iw); rm(iw)
  }
  print(ij); rm(ij)
}

int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_F_K6_conformal_alpha_0.95, c(1, 2), mean)
int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean = apply(int_fore_subnational_err_M_K6_conformal_alpha_0.95, c(1, 2), mean)

colnames(int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean) = colnames(int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean) = c("ECP", "CPD", "score")
rownames(int_fore_subnational_err_F_K6_conformal_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_K6_conformal_alpha_0.95_mean) = 1:15

