###########################
# split conform prediction
###########################
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/FANOVA/"

##################
# Dataset entries
##################

year = 1975:2023
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

# dir.d is the working directory where data are stored
dir.d <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/updated_data/"
female_prefecture_qx <- male_prefecture_qx <- total_prefecture_qx <- list()

for (ik in seq_along(state)) {
  
  # -------------------- FEMALE --------------------
  female_data <- read.table(
    paste0(dir.d, "female_prefecture_", ik, ".txt"),
    header = TRUE, skip = 2
  )
  female_data <- female_data[female_data$Year >= 1975, ]
  female_prefecture_qx[[ik]] <- t(matrix(female_data$qx, n_age, n_year))
  
  
  # -------------------- MALE --------------------
  male_data <- read.table(
    paste0(dir.d, "male_prefecture_", ik, ".txt"),
    header = TRUE, skip = 2
  )
  male_data <- male_data[male_data$Year >= 1975, ]
  male_prefecture_qx[[ik]] <- t(matrix(male_data$qx, n_age, n_year))
  
  
  # -------------------- TOTAL --------------------
  total_data <- read.table(
    paste0(dir.d, "total_prefecture_", ik, ".txt"),
    header = TRUE, skip = 2
  )
  total_data <- total_data[total_data$Year >= 1975, ]
  total_prefecture_qx[[ik]] <- t(matrix(total_data$qx, n_age, n_year))
  
  # Print progress
  print(paste("Finished prefecture:", ik))
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

###############################
# Apply the CDF transformation
###############################

source("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/auxiliary_source/CDF_transformation.R")

Transformed_female <- list()
Transformed_male <- list()
for(i in 1:n_states)
{
  Transformed_female[[i]]   <- t(cdf_transformation(t(All_Japan_female_qx[[i]])/10^5,year))
  Transformed_male[[i]]     <- t(cdf_transformation(t(All_Japan_male_qx[[i]])/10^5,year))
}

all_unconstrained_female<-t(list.cbind(Transformed_female))
all_unconstrained_male<-t(list.cbind(Transformed_male))

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
#########################
# train_set: 1:16
# validation_set: 17:33
# test_set: 34:49
#########################

select_k <- function(tau, eigenvalue)
{
  
  k_max = length(eigenvalue)
  
  k_all = rep(0, k_max-1)
  
  for(k in 1:(k_max-1))
    
  {
    
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
    
  }
  
  K_hat = which.min(k_all)
  
  return(K_hat)
  
}

# interval score

interval_score <- function(holdout, lb, ub, alpha)
{
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
  cpd = abs(cover - (1 - alpha))
  return(c(cover, cpd, mean(score)))
}

# Since for the FANOVA decomposition you used all data for all prefectures and gender at once. 
# You need now to create a list of 47 prefectures of matrices n_age by n_year for both residuals and fixed effects. 


# Split residuals per prefecture
male_prefecture_res_means <- lapply(1:length(part_list), 
                                    function(k){Res1_means[part_list[[k]], ]})
female_prefecture_res_means <-lapply(1:length(part_list), 
                                     function(k){Res2_means[part_list[[k]], ]})

#split the fixed components by prefecture
male_prefecture_fixed_means <-lapply(1:length(part_list), 
                                     function(k){Fixed_part_means_1[part_list[[k]], ]})

female_prefecture_fixed_means <-lapply(1:length(part_list), 
                                       function(k){Fixed_part_means_2[part_list[[k]], ]})

######################################################################################################
# point forecast (FANOVA)
######################################################################################################
Pref_forecast_curves<-function(fixed_com,Residuals_f,
                               est_method = c("lrc", "cov"),
                               fh = 30, 
                               B = 1000, 
                               prediction_method=c("ARIMA","VAR","ets"),select_K=c("Fixed","EVR"), K=6){
  med_polish_resi=t(Residuals_f)
  if(est_method == "lrc"){
    # estimate long-run covariance by kernel sandwich estimator
    med_polish_resi_lrc = long_run_covariance_estimation(med_polish_resi)
  }else if(est_method == "cov"){
    # estimate empirical covariance function
    med_polish_resi_lrc = cov(t(med_polish_resi))
  }
  # perform eigen-decomposition
  med_polish_resi_eigen = eigen(med_polish_resi_lrc)
  
  if(select_K=="Fixed"){
    
    retain_component = K
    
  }else if(select_K=="EVR"){
    # determine retained number of components via eigenvalue ratio
    lambda_val = med_polish_resi_eigen$values
    retain_component = select_k(tau = 10^-2, eigenvalue = lambda_val)
  }
  
  var_total_variations = (sum(med_polish_resi_eigen$values[1:retain_component])/sum(med_polish_resi_eigen$values))*100
  
  
  # determine 1st set of basis function and its scores
  med_polish_resi_basis = as.matrix(med_polish_resi_eigen$vectors[,1:retain_component])
  med_polish_resi_score = crossprod(med_polish_resi, med_polish_resi_basis)
  
  # obtain forecasts of PC scores via auto.arima
  med_polish_resi_score_forecast = matrix(NA, retain_component, fh)
  med_polish_resi_score_forecast_boot = array(NA, dim = c(retain_component, fh, B))
  if(prediction_method=="ARIMA"){
    for(ip in 1:retain_component){
      dum = forecast_Arima(object=auto.arima(med_polish_resi_score[,ip]), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip,] = dum$mean
      med_polish_resi_score_forecast_boot[ip,,] = t(dum$sim)
      rm(ip); rm(dum)
    }
  }else if(prediction_method=="VAR"){
    object=med_polish_resi_score
    colnames(object)<-1:dim(object)[2]
    lag=VARselect(y=object,type = "const")$selection[1]
    model_VAR <- VAR(y=object,type = "const",ic="AIC",p=lag)
    pred=predict(model_VAR,n.ahead=fh)$fcst
    for (ip in 1:retain_component) {
      pred1=pred[[ip]]
      med_polish_resi_score_forecast[ip,]=pred1[,1]
    }
    
  }else if(prediction_method=="ets"){
    
    for(ip in 1:retain_component){
      dum = forecast:::forecast.ets(object=ets(y=as.vector(med_polish_resi_score[,ip])), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ip,] = as.vector(dum$mean)
      # med_polish_resi_score_forecast_boot[ip,,] = t(dum$sim)
      rm(ip); rm(dum)
    }
  }
  
  med_polish_resi_forecast = med_polish_resi_basis %*% med_polish_resi_score_forecast
  
  # add the fixed parts
  
  Fixed=t(fixed_com)[,1:fh]
  med_polish_curve_forecast = med_polish_resi_forecast + Fixed
  
  return(list(med_polish_curve_forecast=med_polish_curve_forecast, 
              med_polish_resi_forecast=  med_polish_resi_forecast,TV = var_total_variations))
  
}

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
  fore_validation = matrix(NA, ncol(fdata), (18 - horizon))
  if(fore_method == "CDF")
  {
    for(ij in 1:(18 - horizon))
    {
      fore_validation[,ij] = fore_FANOVA_cdf(Fixed = Fixed[1:(15 + ij),], Residuals = Residuals[1:(15 + ij),],
                                             ncomp_method = method_ncomp, fh = horizon, fmethod = uni_fore_method, est_method = est_method)
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
      fore_val[,ij] <- fore_FANOVA_cdf(Fixed=Fixed[1:(32 + ij),], Residuals = Residuals[1:(32 + ij),],
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
rownames(int_fore_subnational_err_F_EVR_conformal_mean) = rownames(int_fore_subnational_err_M_EVR_conformal_mean) = 1:16

# K = 6

int_fore_subnational_err_F_K6_conformal = int_fore_subnational_err_M_K6_conformal = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:16)
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
rownames(int_fore_subnational_err_F_K6_conformal_mean) = rownames(int_fore_subnational_err_M_K6_conformal_mean) = 1:16

####################
## level_sig = 0.95
####################

# EVR

int_fore_subnational_err_F_EVR_conformal_alpha_0.95 = int_fore_subnational_err_M_EVR_conformal_alpha_0.95 = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:16)
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
rownames(int_fore_subnational_err_F_EVR_conformal_alpha_0.95_mean) = rownames(int_fore_subnational_err_M_EVR_conformal_alpha_0.95_mean) = 1:16

# K = 6

int_fore_subnational_err_F_K6_conformal_alpha_0.95 = int_fore_subnational_err_M_K6_conformal_alpha_0.95 = array(NA, dim = c(16, 3, 47))
for(ij in 1:47)
{
  for(iw in 1:16)
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
