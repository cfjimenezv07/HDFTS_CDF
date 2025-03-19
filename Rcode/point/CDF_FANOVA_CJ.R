################################
# CDf transformation + FM-ANOVA
################################

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
n_populations = 2
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

# dir.d is the working directory where data are stored

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

###############################
# Apply the CDF transformation
###############################

source("CDF_transformation.R")

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
FANOVA_means <- hdftsa::FANOVA(data_pop1=t(all_unconstrained_male),data_pop2=t(all_unconstrained_female),year,new_age,n_states,n_populations)

#This function computes the functional residuals after removing the deterministic components
# obtained from the FANOVA function.

Residuals_means<-hdftsa::Two_way_Residuals_means(data_pop1=t(all_unconstrained_male), data_pop2=t(all_unconstrained_female), year, new_age, n_states, n_populations)

Res1_means=Residuals_means$residuals1_mean
Res2_means=Residuals_means$residuals2_mean
Residuals_mean<-cbind(Res1_means,Res2_means)

# Reconstructed data
RR<-Residuals_means$rd #Matrix with the original data reconstructed from the FMP decomposition

#  It's the proof of the reconstruction of the residuals.
Residuals_means$R #The result should be a vector with two entries TRUE, TRUE.

#Indicating that after adding both deterministic and time-varying components the FTS are recovered.
Fixed_part_means<-Residuals_means$Fixed_comp_mean # deterministic components to be added up after forecasting
nn_age <- length(new_age)
Fixed_part_means_1 <- Fixed_part_means[,1:nn_age]
Fixed_part_means_2 <- Fixed_part_means[,(nn_age+1):(2*nn_age)]

################################################################################
# Computation of the point forecasts based on functional mean ANOVA (FM-ANOVA)
################################################################################

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

# do the forecasting and compute the errors

Boot_female_EVR <- list()
Boot_female_K <- list()
Boot_male_EVR <- list()
Boot_male_K <- list()
max_h=16
no_core=detectCores()-2
for (i in 1:n_states)
{
  pref=(n_year * i - (n_year - 1)):(n_year * i)
  n_training_ini =length(pref)-max_h #32
  fixed_com1 = Fixed_part_means[pref,]
  Residuals_f1 = Residuals_mean[pref,]
  registerDoMC(no_core)
  # iwk is the expanding window. Therefore you iterate max_h times. Boot is a list of matrices 2*n_age by K=16,...,1.
  boot_EVR = foreach(iwk = 1:max_h ) %dopar% Pref_forecast_curves(fixed_com = fixed_com1[1:(n_training_ini - 1 + iwk),],
                                                              Residuals_f = Residuals_f1[1:(n_training_ini - 1 + iwk),],
                                                              est_method = "cov",
                                                              fh = (max_h -iwk+1), B = 1000,
                                                              prediction_method="ets",select_K="EVR", K=6)$med_polish_curve_forecast
  boot_K = foreach(iwk = 1:max_h ) %dopar% Pref_forecast_curves(fixed_com = fixed_com1[1:(n_training_ini - 1 + iwk),],
                                                                  Residuals_f = Residuals_f1[1:(n_training_ini - 1 + iwk),],
                                                                  est_method = "cov",
                                                                  fh = (max_h -iwk+1), B = 1000,
                                                                  prediction_method="ets",select_K="Fixed", K=6)$med_polish_curve_forecast
  # Create empty lists
  Forecast_female_EVR <- list()
  Forecast_male_EVR <- list()
  Forecast_female_K <- list()
  Forecast_male_K <- list()


  # Extract Forecast_female and Forecast_male from boot
  for (h in seq_along(boot_K)) {
    Forecast_female_EVR[[h]] <- boot_EVR[[h]][111:220, ]  # Extract rows 111 to 220
    Forecast_male_EVR[[h]] <- boot_EVR[[h]][1:110, ]      # Extract rows 1 to 110

    Forecast_female_K[[h]] <- boot_K[[h]][111:220, ]  # Extract rows 111 to 220
    Forecast_male_K[[h]] <- boot_K[[h]][1:110, ]
  }

  
  Forecast_female_transformed_EVR <- list()
  Forecast_male_transformed_EVR <- list()
  Forecast_female_transformed_K <- list()
  Forecast_male_transformed_K <- list()
  # Apply transformations
  for (h in seq_along(Forecast_female_EVR)) {
    fore_val_EVR_1 <- as.matrix(Forecast_female_EVR[[h]])
    fore_val_EVR_2 <- as.matrix(Forecast_male_EVR[[h]])

    fore_val_K_1 <- as.matrix(Forecast_female_K[[h]])
    fore_val_K_2 <- as.matrix(Forecast_male_K[[h]])

    n_age <- nrow(fore_val_EVR_1)  # Get number of ages (should be 110)
    n_cols <- ncol(fore_val_EVR_1) # Number of columns (varies from 16 to 1)

    # Initialize matrices for EVR
    f_x_t_star_fore_EVR_1 <- f_x_t_star_fore_EVR_2 <-
      d_x_t_star_fore_EVR_1 <- d_x_t_star_fore_EVR_2 <- matrix(NA, (n_age+1), n_cols)

    for (ik in 1:n_cols) {
      f_x_t_star_fore_EVR_1[, ik] <- c(invlogit(fore_val_EVR_1[, ik]), 1)
      f_x_t_star_fore_EVR_2[, ik] <- c(invlogit(fore_val_EVR_2[, ik]), 1)

      d_x_t_star_fore_EVR_1[, ik] <- c(f_x_t_star_fore_EVR_1[1, ik], diff(f_x_t_star_fore_EVR_1[, ik]))
      d_x_t_star_fore_EVR_2[, ik] <- c(f_x_t_star_fore_EVR_2[1, ik], diff(f_x_t_star_fore_EVR_2[, ik]))
    }

    # Initialize matrices for K
    f_x_t_star_fore_K_1 <- f_x_t_star_fore_K_2 <-
      d_x_t_star_fore_K_1 <- d_x_t_star_fore_K_2 <- matrix(NA, (n_age+1), n_cols)

    for (ik in 1:n_cols) {
      f_x_t_star_fore_K_1[, ik] <- c(invlogit(fore_val_K_1[, ik]), 1)
      f_x_t_star_fore_K_2[, ik] <- c(invlogit(fore_val_K_2[, ik]), 1)

      d_x_t_star_fore_K_1[, ik] <- c(f_x_t_star_fore_K_1[1, ik], diff(f_x_t_star_fore_K_1[, ik]))
      d_x_t_star_fore_K_2[, ik] <- c(f_x_t_star_fore_K_2[1, ik], diff(f_x_t_star_fore_K_2[, ik]))
    }


    # Store transformed matrices in lists
    Forecast_female_transformed_EVR[[h]] <-  d_x_t_star_fore_EVR_1* 10^5
    Forecast_male_transformed_EVR[[h]] <-    d_x_t_star_fore_EVR_2* 10^5
    Forecast_female_transformed_K[[h]] <-  d_x_t_star_fore_K_1* 10^5
    Forecast_male_transformed_K[[h]] <-    d_x_t_star_fore_K_2* 10^5
  }

  Boot_female_EVR[[i]] <- Forecast_female_transformed_EVR
  Boot_female_K[[i]]   <- Forecast_female_transformed_K
  Boot_male_EVR[[i]]   <- Forecast_male_transformed_EVR
  Boot_male_K[[i]]     <- Forecast_male_transformed_K

}

##############################################################################
#Compute the forecast errors
##############################################################################
source("forecast_errors.R")

compute_error <- function(forecast, true)
{
  KLD    <- KLD(forecast, true)
  JS2    <- JSD_geom(forecast, true)
  return(c(KLD,JS2))
}

pairwise_errors <- function(j, forecast, true) {
  forecast <- forecast[, j]
  if (is.null(dim(true))) true <- matrix(true, ncol = 1)
  true <- true[, j]
  compute_error(forecast, true)
}

# function to pass from errors at each EW to errors per forecast horizon
create_new_list <- function(A) {
  # Number of columns in the largest matrix (A[[1]]) is 16
  num_columns <- 16
  
  # Initialize an empty list for B
  B <- vector("list", num_columns)
  
  # Loop over each column index (1 to 16)
  for (col in 1:num_columns) {
    # For each column in B, collect the col-th column from each matrix in A
    B[[col]] <- sapply(A, function(mat) {
      if (ncol(mat) >= col) {
        # If the matrix has enough columns, extract the column as a vector
        mat[, col]
      } else {
        # If the matrix does not have enough columns, return NA
        rep(NA, nrow(mat))
      }
    })
  }
  
  return(B)
}

# Compute the errors per gender and error type
errors_KLD_female_EVR <- matrix(NA,16,47)
errors_KLD_male_EVR <- matrix(NA,16,47)

errors_JSD_female_EVR <- matrix(NA,16,47)
errors_JSD_male_EVR <- matrix(NA,16,47)

errors_KLD_female_K <- matrix(NA,16,47)
errors_KLD_male_K <- matrix(NA,16,47)

errors_JSD_female_K <- matrix(NA,16,47)
errors_JSD_male_K <- matrix(NA,16,47)
for (i in 1:n_states) {
  
 forecast_female_transformed_EVR <- Boot_female_EVR[[i]]
 forecast_female_transformed_K   <- Boot_female_K[[i]]
 forecast_male_transformed_EVR   <- Boot_male_EVR[[i]]
 forecast_male_transformed_K     <- Boot_male_K[[i]]
  
  
  # Holdout data
  Holdout_female_1   <- female_prefecture_dx[[i]]
  Holdout_male_1     <- male_prefecture_dx[[i]]
  
  
  if(any(Holdout_female_1 == 0 | Holdout_male_1 == 0)){
    Holdout_female = replace(x = Holdout_female_1, list = which(Holdout_female_1 == 0), values = 10^-5)
    Holdout_male   = replace(x = Holdout_male_1, list = which(Holdout_male_1 == 0), values = 10^-5)
  }else{
    Holdout_female = Holdout_female_1
    Holdout_male = Holdout_male_1
  }
  
  # This are the list of errors for each expanding window (EW). That is for EW=1, I compute the pairwise error between forecasts h=1,...,h=16
  # with holdout data for years 33 to 48.
  errors_female_EVR_list <- vector("list", 16)  # Each entry will store 2 errors for each horizon
  errors_male_EVR_list   <- vector("list", 16)
  errors_female_K_list <- vector("list", 16)  # Each entry will store 2 errors for each horizon
  errors_male_K_list   <- vector("list", 16)
  
  # Iterate through the Expanding windows (h = 1 to 16)
  for (h in seq_along(forecast_female_transformed_EVR)) {
    # Get forecast matrices
    forecast_female_EVR <- as.matrix(forecast_female_transformed_EVR[[h]])
    forecast_male_EVR   <- as.matrix(forecast_male_transformed_EVR[[h]])
    forecast_female_K <- as.matrix(forecast_female_transformed_K[[h]])
    forecast_male_K   <- as.matrix(forecast_male_transformed_K[[h]])
    
    
    # Define corresponding holdout years (column indices)
    holdout_start <- 33 + (h - 1)
    holdout_end   <- 48
    holdout_cols  <- holdout_start:holdout_end
    
    # Extract the appropriate part of the holdout dataset
    true_female <- as.matrix(Holdout_female[, holdout_cols])
    true_male   <- as.matrix(Holdout_male[, holdout_cols])
    
    # Compute errors for each column of forecasts
    error_list_female_EVR <- sapply(seq_len(ncol(forecast_female_EVR)), function(j) {
      pairwise_errors(j, forecast_female_EVR, true_female)
    })
    
    error_list_male_EVR <- sapply(seq_len(ncol(forecast_male_EVR)), function(j) {
      pairwise_errors(j, forecast_male_EVR, true_male)
    })
    
    error_list_female_K <- sapply(seq_len(ncol(forecast_female_K)), function(j) {
      pairwise_errors(j, forecast_female_K, true_female)
    })
    
    error_list_male_K <- sapply(seq_len(ncol(forecast_male_K)), function(j) {
      pairwise_errors(j, forecast_male_K, true_male)
    })
    
    # Store computed errors (2 × k matrix for this h)
    errors_female_EVR_list[[h]] <- error_list_female_EVR
    errors_male_EVR_list[[h]]   <- error_list_male_EVR
    
    errors_female_K_list[[h]] <- error_list_female_K
    errors_male_K_list[[h]]   <- error_list_male_K
  }
  
  # Now we need to aggregate errors over all h. That is we need to cbind for h=1 (in the 16 EW),h=2(in the 15 EW),...,h=16(in 1 EW). 
  
  # what this does is to get the h=1 from each expanding window (16 times). Then the h=2 (15 times)...h=16 (one time). Then create matrices with these errors. 
  
  new_female_EVR <-create_new_list(errors_female_EVR_list)
  new_female_K <-create_new_list(errors_female_K_list)
  new_male_EVR <-create_new_list(errors_male_EVR_list)
  new_male_K <-create_new_list(errors_male_K_list)
  
  # Initialize final matrices
  errors_female_EVR_final <- matrix(NA, 2, 16)  # 2 error metrics × 16 forecast horizons
  errors_male_EVR_final   <- matrix(NA, 2, 16)
 
  errors_female_K_final <- matrix(NA, 2, 16)  # 2 error metrics × 16 forecast horizons
  errors_male_K_final   <- matrix(NA, 2, 16)
   
  # Compute mean errors for each horizon h
  for (h in 1:16) {
    # Compute column-wise mean over available computations
    errors_female_EVR_final[, h] <- rowMeans( new_female_EVR[[h]], na.rm = TRUE)
    errors_male_EVR_final[, h]   <- rowMeans( new_male_EVR[[h]], na.rm = TRUE)
    errors_female_K_final[, h] <- rowMeans( new_female_K[[h]], na.rm = TRUE)
    errors_male_K_final[, h]   <- rowMeans( new_male_K[[h]], na.rm = TRUE)
  }
  
  errors_KLD_female_EVR[,i]  = errors_female_EVR_final[1, ]
  errors_KLD_male_EVR[,i]  = errors_male_EVR_final[1, ]
  errors_JSD_female_EVR[,i]  = errors_female_EVR_final[2, ]
  errors_JSD_male_EVR[,i]  = errors_male_EVR_final[2, ]
  
  errors_KLD_female_K[,i]  = errors_female_K_final[1, ]
  errors_KLD_male_K[,i]  = errors_male_K_final[1, ]
  errors_JSD_female_K[,i]  = errors_female_K_final[2, ]
  errors_JSD_male_K[,i]  = errors_male_K_final[2, ]
}

#errors 16 by 47

errors_pref_EVR <-  list(errors_KLD_female_EVR,errors_JSD_female_EVR,errors_KLD_male_EVR,errors_JSD_male_EVR)
errors_pref_K <-  list(errors_KLD_female_K,errors_JSD_female_K,errors_KLD_male_K,errors_JSD_male_K)

# see the errors averaged by prefecture per error type

errors_EVR <- cbind(apply(errors_KLD_female_EVR,1,mean),apply(errors_JSD_female_EVR,1,mean),apply(errors_KLD_male_EVR,1,mean),apply(errors_JSD_male_EVR,1,mean))
colnames(errors_EVR) <- c("KLD_Female","JSD_Female","KLD_Male","JSD_Male")
errors_EVR
colMeans(errors_EVR)
apply(errors_EVR, 2, median)

errors_K <- cbind(apply(errors_KLD_female_K,1,mean),apply(errors_JSD_female_K,1,mean),apply(errors_KLD_male_K,1,mean),apply(errors_JSD_male_K,1,mean))
colnames(errors_K) <- c("KLD_Female","JSD_Female","KLD_Male","JSD_Male")
errors_K
colMeans(errors_K)
apply(errors_K, 2, median)

