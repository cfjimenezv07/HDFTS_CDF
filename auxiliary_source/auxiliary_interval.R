##################################################################
# determining the optimal tuning parameter via standard deviation
##################################################################
# Datasets
year = 1975:2023
n_year = length(year)
age = 0:110
n_age = length(age)
n_prefectures=47
state = c("Hokkaido", 
          "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
          "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
          "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
          "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

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
  female_prefecture_dx[[iw]] = (female_prefecture_dum)
  male_prefecture_dx[[iw]]   = (male_prefecture_dum)
  total_prefecture_dx[[iw]]  = (total_prefecture_dum)
  rm(female_prefecture_dum); rm(male_prefecture_dum); rm(total_prefecture_dum)
  print(iw); rm(iw)
}

# selecting the number of components
# tau: a tuning parameter
# eigenvalue: estimated eigenvalues

select_K <- function(tau, eigenvalue)
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

# tune_para: tuning parameter
# resi_mat: residual functions
# sd_val_input: functional standard deviation (pointwise)
# alpha_level: level of significance
# PI_type: pointwise or uniform

years= 1975:2023
ages = 1:111

tune_para_find_function <- function(tune_para, resi_mat, sd_val_input, alpha_level, PI_type)
{
  n_age = nrow(resi_mat)
  if(PI_type == "pointwise")
  {
    ind = matrix(NA, n_age, ncol(resi_mat))
    for(iw in 1:ncol(resi_mat))
    {
      ind[,iw] = ifelse(between(resi_mat[,iw], -tune_para * sd_val_input, tune_para * sd_val_input), 1, 0)
      rm(iw)
    }
    ecp = sum(ind)/(n_age * ncol(resi_mat))
  }
  else if(PI_type == "uniform")
  {
    ind = vector("numeric", ncol(resi_mat))
    for(iw in 1:ncol(resi_mat))
    {
      ind[iw] = ifelse(all(between(resi_mat[,iw], -tune_para * sd_val_input, 
                                   tune_para * sd_val_input)), 1, 0)
      rm(iw)
    }
    ecp = sum(ind)/ncol(resi_mat)
  }
  else
  {
    warning("PI type must either be pointwise or uniform.")
  }
  rm(ind)
  return(abs(ecp - alpha_level))
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

# uniform cpd

uniform_cpd <- function(holdout, lb, ub, alpha)
{
  lb_ind = ifelse(any(holdout < lb), 1, 0)
  ub_ind = ifelse(any(holdout > ub), 1, 0)
  cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/ncol(holdout)
  cpd = abs(cover - (1 - alpha))
  return(c(cover, cpd))
}

####################################
# univariate functional time series 
####################################

# data_set: life-table death counts
# ncomp_method: way of selecting the number of components
# fh: forecast horizon
# fmethod: ets

fore_national_cdf <- function(data_set, ncomp_method, fh, fmethod)
{
  data = data_set/10^5
  data_cumsum_dum = matrix(NA, nrow(data), ncol(data))
  for(ij in 1:nrow(data))
  {
    data_cumsum_dum[ij,] = cumsum(data[ij,])
    rm(ij)
  }
  
  # check if any cumsum values equal to 0
  
  if(any(data_cumsum_dum == 0))
  {
    data_cumsum = replace(data_cumsum_dum, which(data_cumsum_dum == 0), 10^-5)
  }
  else
  {
    data_cumsum = data_cumsum_dum
  }
  rm(data_cumsum_dum)
  
  # logit transformation
  
  data_cumsum_logit = matrix(NA, nrow(data), (ncol(data) - 1))
  for(ij in 1:nrow(data))
  {
    data_cumsum_logit[ij,] = logit(data_cumsum[ij, 1:(ncol(data) - 1)])
    rm(ij)
  }
  rm(data_cumsum)
  rownames(data_cumsum_logit) = years[1:nrow(data)]
  
  # fitting a functional time series forecasting method
  
  if(ncomp_method == "EVR")
  {
    ncomp = select_K(tau = 10^-3, eigenvalue = (svd(data_cumsum_logit)$d)^2)
  }
  else if(ncomp_method == "provide")
  {
    ncomp = 6
  }
  else
  {
    warning("The number of components is required.")
  }
  data_cumsum_logit_fore = forecast(object = ftsm(fts(ages[1:110], t(data_cumsum_logit)), order = ncomp), 
                                    h = fh, method = fmethod, B = 399)
  
  # h-step-ahead mean forecast
  
  data_cumsum_logit_fore_add = c(invlogit(data_cumsum_logit_fore$mean$y[,fh]), 1)
  data_cumsum_logit_fore_add_diff = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
  
  return(data_cumsum_logit_fore_add_diff * 10^5)
}

####################################
# multivariate functional time series
####################################

# data_input: array
# ncomp_method: EVR or K = 6
# fh: forecast horizon
# fore_method: forecasting method

MFTS_model <- function(data_input, ncomp_method, fh, fore_method)
{
    n_age  = dim(data_input)[1]
    n_year = dim(data_input)[2]
    n_pop  = dim(data_input)[3]
    
    # clean the data
    
    data_set_array = array(NA, dim = c(n_age, n_year, n_pop))
    if(any(!is.finite(data_input)))
    {
      for(iw in 1:n_pop)
      {
        for(ij in 1:n_age)
        {
          data_set_array[ij,,iw] = na.interp(data_input[ij,,iw])
        }
      }
    }
    else
    {
      data_set_array = data_input
    }
    
    rowmeans_object = sd_object = decenter_object = list()
    for(ik in 1:n_pop)
    {
      # compute mean and sd function
      rowmeans_object[[ik]] = rowMeans(data_set_array[,,ik], na.rm = TRUE)
      sd_object[[ik]] = apply(data_set_array[,,ik], 1, sd, na.rm = TRUE)
      
      # de-center functional data
      decenter_object[[ik]] = t(scale(t(data_set_array[,,ik]), center = TRUE, scale = TRUE))
    }
    
    comb_object = do.call(rbind, decenter_object)
    # comb_object = do.call(rbind, comb_object_raw)
    colnames(comb_object) = 1:ncol(comb_object)
    
    eigen_value = eigen(cov(t(comb_object)))$values
    if(ncomp_method == "EVR")
    {
      ncomp = select_K(tau = 10^-2, eigenvalue = eigen_value)
    }
    else if(ncomp_method == "provide")
    {
      ncomp = 6
    }
    else
    {
      warning("The number of components is required.")
    }
    fore_ftsm = forecast(ftsm(fts(1:nrow(comb_object), comb_object), order = ncomp), h = fh, 
                         method = fore_method)
    res_fore = as.matrix(fore_ftsm$mean$y[,fh] * do.call(c, sd_object) + do.call(c, rowmeans_object))
    return(res_fore)
}
  
####################################################
# Multivariate functional time series decomposition
####################################################

# data_F: female data
# data_M: male data
# fh: forecast horizon
# fmethod: forecasting method
# object_interest: point or interval forecasts
# alpha: level of significance
# method_ncomp: K = 6 or fixed

fore_national_cdf_MFTS <- function(data_set_F, data_set_M, fh, fmethod, method_ncomp)
{
    data_F = data_set_F/10^5
    data_M = data_set_M/10^5
    
    data_cumsum_dum_F = data_cumsum_dum_M = matrix(NA, nrow(data_F), ncol(data_F))
    for(ij in 1:nrow(data_F))
    {
      data_cumsum_dum_F[ij,] = cumsum(data_F[ij,])
      data_cumsum_dum_M[ij,] = cumsum(data_M[ij,])
      rm(ij)
    }
    
    if(any(data_cumsum_dum_F == 0))
    {
      data_cumsum_F = replace(data_cumsum_dum_F, which(data_cumsum_dum_F == 0), 10^-5)
    }
    else
    {
      data_cumsum_F = data_cumsum_dum_F
    }
    if(any(data_cumsum_dum_M == 0))
    {
      data_cumsum_M = replace(data_cumsum_dum_M, which(data_cumsum_dum_M == 0), 10^-5)
    }
    else
    {
      data_cumsum_M = data_cumsum_dum_M
    }
    rm(data_cumsum_dum_F); rm(data_cumsum_dum_M)
    
    data_cumsum_logit_F = data_cumsum_logit_M = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
    for(ij in 1:nrow(data_F))
    {
      data_cumsum_logit_F[ij,] = logit(data_cumsum_F[ij, 1:(ncol(data_F) - 1)])
      data_cumsum_logit_M[ij,] = logit(data_cumsum_M[ij, 1:(ncol(data_M) - 1)])
      rm(ij)
    }
    rownames(data_cumsum_logit_F) = rownames(data_cumsum_logit_M) = years[1:nrow(data_F)]
    
    data_comb = array(NA, dim = c((ncol(data_F) - 1), nrow(data_F), 2))
    data_comb[,,1] = t(data_cumsum_logit_F)
    data_comb[,,2] = t(data_cumsum_logit_M)
    
    dum = MFTS_model(data_input = data_comb, ncomp_method = method_ncomp, fh = fh, fore_method = fmethod)
    
    data_cumsum_logit_F_fore = dum[1:(ncol(data_F) - 1),]
    data_cumsum_logit_M_fore = dum[ncol(data_F):(2 * (ncol(data_F) - 1)),]
    rm(dum)
    
    data_cumsum_logit_F_fore_add = c(invlogit(data_cumsum_logit_F_fore), 1)
    data_cumsum_logit_M_fore_add = c(invlogit(data_cumsum_logit_M_fore), 1)
    
    data_cumsum_logit_F_fore_add_diff = c(data_cumsum_logit_F_fore_add[1], diff(data_cumsum_logit_F_fore_add))
    data_cumsum_logit_M_fore_add_diff = c(data_cumsum_logit_M_fore_add[1], diff(data_cumsum_logit_M_fore_add))
    return(list(mfts_fore_F = data_cumsum_logit_F_fore_add_diff * 10^5, 
                mfts_fore_M = data_cumsum_logit_M_fore_add_diff * 10^5))
}

#################################
# High-dimensional FPCA (HDFPCA)
#################################

sco.resamp = ftsa:::sco.resamp

# fdata_F: female data
# fdata_M: male data
# horizon: forecast horizon
# first_order: number of retained components
# second_order: number of retained components
# transformation: direct, clr or CDF

hdfpca_fun <- function(fdata_F, fdata_M, horizon, first_order, second_order, transformation)
{
    forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    if(transformation == "direct")
    {
        for(ij in 1:(17 - horizon))
        {
            data_F = fdata_F[1:(nrow(fdata_F)-1+ij),]
            data_M = fdata_M[1:(nrow(fdata_F)-1+ij),]
            data_comb = list()
            data_comb[[1]] = t(data_F)
            data_comb[[2]] = t(data_M)
            
            fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
            forecast_val_F[,ij] = (fore_val[[1]])[,horizon]
            forecast_val_M[,ij] = (fore_val[[2]])[,horizon]
            rm(ij)
        }
    }
    else if(transformation == "clr")
    {
        for(ij in 1:(17 - horizon))
        {
            data_F = as.matrix(clr(fdata_F[1:(nrow(fdata_F)-1+ij),]))
            data_M = as.matrix(clr(fdata_M[1:(nrow(fdata_F)-1+ij),]))
            data_comb = list()
            data_comb[[1]] = t(data_F)
            data_comb[[2]] = t(data_M)
            
            fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
            forecast_val_F[,ij] = as.numeric(clrInv((fore_val[[1]])[,horizon])) * 10^5
            forecast_val_M[,ij] = as.numeric(clrInv((fore_val[[2]])[,horizon])) * 10^5
            rm(ij)
        }
    }
    else if(transformation == "CDF")
    {
        for(ijk in 1:(17 - horizon))
        {
            data_F = fdata_F[1:(nrow(fdata_F)-1+ijk),]/10^5
            data_M = fdata_M[1:(nrow(fdata_F)-1+ijk),]/10^5
            
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
            
            forecast_val_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
            forecast_val_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
            rm(ijk); rm(data_F); rm(data_M)
        }
    }
    holdout_val_dum_F = t(matrix(fdata_F[(nrow(fdata_F) + horizon):(nrow(fdata_F)+16),], length((nrow(fdata_F) + horizon):(nrow(fdata_F)+16)), ncol(fdata_F)))
    holdout_val_dum_M = t(matrix(fdata_M[(nrow(fdata_F) + horizon):(nrow(fdata_F)+16),], length((nrow(fdata_F) + horizon):(nrow(fdata_F)+16)), ncol(fdata_F)))
    if(any(holdout_val_dum_F == 0))
    {
        holdout_val_F = replace(holdout_val_dum_F, which(holdout_val_dum_F == 0), 10^-5)
    }
    else
    {
        holdout_val_F = holdout_val_dum_F
    }
    if(any(holdout_val_dum_M == 0))
    {
        holdout_val_M = replace(holdout_val_dum_M, which(holdout_val_dum_M == 0), 10^-5)
    }
    else
    {
        holdout_val_M = holdout_val_dum_M
    }
    rm(holdout_val_dum_F); rm(holdout_val_dum_M)
    
    KL_div_val_F = JS_div_val_F = L1_dist_F = L2_dist_F = 
    KL_div_val_M = JS_div_val_M = L1_dist_M = L2_dist_M = vector("numeric", (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        # symmetric KL dist
        
        KL_div_val_F[ij] = mean(KLdiv(cbind(forecast_val_F[,ij], holdout_val_F[,ij]))[2:3])
        
        # Jensen-Shannon dist
        
        JS_div_val_F[ij] = sqrt(mean(KLdiv(cbind(forecast_val_F[,ij], apply(cbind(forecast_val_F[,ij], holdout_val_F[,ij]), 1, geometric.mean)))[2:3]))
        
        # L1 dist
        
        L1_dist_F[ij] = wasserstein1d(a = forecast_val_F[,ij], b = holdout_val_F[,ij], p = 1)
        
        # L2 dist
        
        L2_dist_F[ij] = wasserstein1d(a = forecast_val_F[,ij], b = holdout_val_F[,ij], p = 2)
        
        # symmetric KL dist
        
        KL_div_val_M[ij] = mean(KLdiv(cbind(forecast_val_M[,ij], holdout_val_M[,ij]))[2:3])
        
        # Jensen-Shannon dist
        
        JS_div_val_M[ij] = sqrt(mean(KLdiv(cbind(forecast_val_M[,ij], apply(cbind(forecast_val_M[,ij], holdout_val_M[,ij]), 1, geometric.mean)))[2:3]))
        
        # L1 dist
        
        L1_dist_M[ij] = wasserstein1d(a = forecast_val_M[,ij], b = holdout_val_M[,ij], p = 1)
        
        # L2 dist
        
        L2_dist_M[ij] = wasserstein1d(a = forecast_val_M[,ij], b = holdout_val_M[,ij], p = 2)
        rm(ij)
    }
    err_F = c(mean(KL_div_val_F), mean(JS_div_val_F), sqrt(mean(L1_dist_F^2)), sqrt(mean(L2_dist_F^2)))
    err_M = c(mean(KL_div_val_M), mean(JS_div_val_M), sqrt(mean(L1_dist_M^2)), sqrt(mean(L2_dist_M^2)))
    
    return(list(forecast_pdf_F = forecast_val_F, forecast_pdf_M = forecast_val_M,
                holdout_pdf_F = holdout_val_F, holdout_pdf_M = holdout_val_M,
                err_F = err_F, err_M = err_M))
}

#########
# HDFPCA
#########

# fdata_F: female data
# fdata_M: male data
# horizon: forecast horizon
# first_order: number of components in the first stage
# second_order: number of components in the second stage
# transformation: CLR, CDF, direct

###########################
# forecast HDFPCA function
###########################


# object: functional data object
# h: forecast horizon
# level: level of prediction
# B: number of replicates

forecast.hdfpca <- function(object, h, level, B) 
{
  order = object$order
  r = object$r
  m = object$m
  p = object$p
  n = dim(object$y[[1]])[2]
  resid <- list()
  for(im in 1:m) 
  {
    resid[[im]] <- object$y[[im]] - object$fitted[[im]]
  }
  mod.fore = load.fore <- list()
  for(io in 1:order) 
  {
    load.fore[[io]] <- array(NA, dim = c(h, r))
    mod.fore[[io]] <- list()
    for(ir in 1:r) 
    {
      mod <- auto.arima(object$model$model2[[io]]$coef[,ir])
      mod.fore[[io]][[ir]] <- forecast(mod, h)
      load.fore[[io]][, ir] <- mod.fore[[io]][[ir]]$mean
    }
  }
  score.fore <- list()
  for(io in 1:order) 
  {
    score.fore[[io]] <- object$model$model2[[io]]$basis[, 1] + object$model$model2[[io]]$basis[, 2:(1 + r)] %*% t(load.fore[[io]])
  }
  score <- list()
  for(im in 1:m) 
  {
    score[[im]] <- sapply(score.fore, "[", ((im - 1) * h + 1):(im * h))
  }
  fun.fore <- list()
  if(h == 1) 
  {
    for(im in 1:m) 
    {
      fun.fore[[im]] <- object$model$model1[[im]]$basis[, 1] + object$model$model1[[im]]$basis[, 2:(1 + order)] %*% score[[im]]
    }
  }
  #else if (h > (n/2)) 
  #{
  #  warning("forecast horizon is too big considering the sample size")
  #}
  else
  {
    for(im in 1:m) 
    {
      fun.fore[[im]] <- object$model$model1[[im]]$basis[, 1] + object$model$model1[[im]]$basis[, 2:(1 + order)] %*% t(score[[im]])
    }
  }
  return(structure(list(forecast = fun.fore), class = "forecast.hdfpca"))
}



hdfpca_fun_fore <- function(fdata_F, fdata_M, horizon, first_order, second_order, transformation)
{
  forecast_val_F_test = forecast_val_M_test = matrix(NA, ncol(fdata_F), (17 - horizon))
  forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (18 - horizon))
  if(transformation == "direct")
  {
    for(ij in 1:(18 - horizon))
    {
      data_F = fdata_F[1:(nrow(fdata_F)-16+ij),]
      data_M = fdata_M[1:(nrow(fdata_F)-16+ij),]
      data_comb = list()
      data_comb[[1]] = t(data_F)
      data_comb[[2]] = t(data_M)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      forecast_val_F[,ij] = (fore_val[[1]])[,horizon]
      forecast_val_M[,ij] = (fore_val[[2]])[,horizon]
      rm(ij)
    }
  }
  else if(transformation == "clr")
  {
    for(ij in 1:(18 - horizon))
    {
      data_F = as.matrix(clr(fdata_F[1:(nrow(fdata_F)-18+ij),]))
      data_M = as.matrix(clr(fdata_M[1:(nrow(fdata_F)-18+ij),]))
      data_comb = list()
      data_comb[[1]] = t(data_F)
      data_comb[[2]] = t(data_M)
      
      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
      forecast_val_F[,ij] = as.numeric(clrInv((fore_val[[1]])[,horizon])) * 10^5
      forecast_val_M[,ij] = as.numeric(clrInv((fore_val[[2]])[,horizon])) * 10^5
      rm(ij)
    }
  }
  else if(transformation == "CDF")
  {
    for(ijk in 1:(18 - horizon))
    {
      data_F = fdata_F[1:(nrow(fdata_F) - 18 + ijk),]/10^5
      data_M = fdata_M[1:(nrow(fdata_F) - 18 + ijk),]/10^5
      
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
      
      forecast_val_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
      forecast_val_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
      rm(ijk); rm(data_F); rm(data_M)
    }
  }else if(transformation == "CDF_test")
  {
    
    for(ijk in 1:(17 - horizon))
    {
      data_F = fdata_F[1:(nrow(fdata_F) - 17 + ijk),]/10^5
      data_M = fdata_M[1:(nrow(fdata_F) - 17 + ijk),]/10^5
      
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
      
      forecast_val_F_test[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
      forecast_val_M_test[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
      rm(ijk); rm(data_F); rm(data_M)
    }
  }
  return(list(forecast_val_F = forecast_val_F, forecast_val_M = forecast_val_M,
              forecast_val_F_test = forecast_val_F_test, forecast_val_M_test = forecast_val_M_test))
}











MLFTS_model<- function(data_input, aux_var, ncomp_method, fh, fore_method)
{
  n_age  = dim(data_input)[1]
  n_year = dim(data_input)[2]
  n_pop  = dim(data_input)[3]
  
  # clean the data
  
  data_set_array = array(NA, dim = c(n_age, n_year, n_pop))
  if(any(!is.finite(data_input)))
  {
    for(iw in 1:n_pop)
    {
      for(ij in 1:n_age)
      {
        data_set_array[ij,,iw] = na.interp(data_input[ij,,iw])
      }
    }
  }
  else
  {
    data_set_array = data_input
  }
  
  # compute the mean function
  
  mean_function_list = list()
  for(ik in 1:n_pop)
  {
    mean_function_list[[ik]] = rowMeans(data_set_array[,,ik], na.rm = TRUE)
    rm(ik)
  }
  
  data_set = array(NA, dim = c(n_age, n_year, n_pop))
  for(ik in 1:n_pop)
  {
    data_set[,,ik] = t(scale(t(data_set_array[,,ik]), center = TRUE, scale = FALSE))
    rm(ik)
  }
  
  if(missing(aux_var)|is.null(aux_var))
  {
    aggregate_data = apply(data_set, c(1, 2), mean)
  }
  else
  {
    aggregate_data = t(aux_var)
  }
  colnames(aggregate_data) = 1:n_year
  rownames(aggregate_data) = 1:n_age
  
  # 1st FPCA
  
  eigen_value_aggregate = eigen(cov(t(aggregate_data)))$values
  if(ncomp_method == "EVR")
  {
    ncomp_aggregate = select_K(tau = 10^-3, eigenvalue = eigen_value_aggregate)
  }
  else if(ncomp_method == "provide")
  {
    ncomp_aggregate = 6
  }
  ftsm_aggregate = ftsm(fts(1:n_age, aggregate_data), order = ncomp_aggregate)
  
  # calculate sum of lambda_k
  sum_lambda_k = sum(eigen_value_aggregate[1:ncomp_aggregate])
  
  # compute the residual trend
  data_residual = array(NA, dim = c(n_age, n_year, n_pop))
  for(iw in 1:n_pop)
  {
    data_residual[,,iw] = data_set[,,iw] - ftsm_aggregate$fitted$y
    colnames(data_residual[,,iw]) = 1:n_year
    rownames(data_residual[,,iw]) = 1:n_age
    rm(iw)
  }
  
  # 2nd FPCA
  
  if(ncomp_method == "EVR")
  {
    ncomp_resi = vector("numeric", n_pop)
    for(iw in 1:n_pop)
    {
      eigen_value_resi = eigen(cov(t(data_residual[,,iw])))$values
      ncomp_resi[iw] = select_K(tau = 10^-3, eigenvalue = eigen_value_resi)
    }
  }
  else if(ncomp_method == "provide")
  {
    ncomp_resi = rep(6, n_pop)
  }
  
  sum_lambda_l = vector("numeric", n_pop)
  for(iw in 1:n_pop)
  {
    eigen_value_resi = eigen(cov(t(data_residual[,,iw])))$values
    
    # calculate sum of lambda_l
    sum_lambda_l[iw] = sum(eigen_value_resi[1:(ncomp_resi[iw])])
  }
  
  ftsm_resi = list()
  for(iw in 1:n_pop)
  {
    ftsm_resi[[iw]] = ftsm(fts(1:n_age, data_residual[,,iw]), order = ncomp_resi[iw])
    rm(iw)
  }
  
  # within-cluster variability
  
  within_cluster_variability = vector("numeric", n_pop)
  for(iw in 1:n_pop)
  {
    within_cluster_variability[iw] = sum_lambda_k/(sum_lambda_k + sum_lambda_l[iw])
  }
  
  # reconstruction
  
  coef_fore = matrix(NA, ncomp_aggregate, fh)
  if(fore_method == "arima")
  {
    for(ik in 1:ncomp_aggregate)
    {
      coef_fore[ik,] = forecast(auto.arima(ftsm_aggregate$coeff[,ik+1]), h = fh)$mean
    }
  }
  else if(fore_method == "ets")
  {
    for(ik in 1:ncomp_aggregate)
    {
      coef_fore[ik,] = forecast(ets(ftsm_aggregate$coeff[,ik+1]), h = fh)$mean
    }
  }
  else
  {
    warning("Forecasting method can either be ARIMA or ETS.")
  }
  rownames(coef_fore) = 1:ncomp_aggregate
  colnames(coef_fore) = 1:fh
  
  if(ncomp_aggregate == 1)
  {
    aggregate_fore = as.matrix(ftsm_aggregate$basis[,2]) %*% matrix(coef_fore, nrow = 1)
  }
  else
  {
    aggregate_fore = ftsm_aggregate$basis[,2:(ncomp_aggregate+1)] %*% coef_fore
  }
  
  # residual forecasts
  
  coef_fore_resi_list = list()
  for(iw in 1:n_pop)
  {
    coef_fore_resi = matrix(NA, ncomp_resi[iw], fh)
    if(fore_method == "arima")
    {
      for(ik in 1:ncomp_resi[iw])
      {
        coef_fore_resi[ik,] = forecast(auto.arima(ftsm_resi[[iw]]$coeff[,ik+1]), h = fh)$mean
      }
    }
    else if(fore_method == "ets")
    {
      for(ik in 1:ncomp_resi[iw])
      {
        coef_fore_resi[ik,] = forecast(ets(ftsm_resi[[iw]]$coeff[,ik+1]), h = fh)$mean
      }
    }
    else
    {
      warning("Forecasting method can either be ARIMA or ETS.")
    }
    coef_fore_resi_list[[iw]] = coef_fore_resi
    rm(iw)
  }
  
  resi_fore = list()
  for(iw in 1:n_pop)
  {
    resi_fore[[iw]] = ftsm_resi[[iw]]$basis[,2:(ncomp_resi[iw] + 1)] %*% coef_fore_resi_list[[iw]]
    rm(iw)
  }
  
  final_fore = list()
  for(iw in 1:n_pop)
  {
    final_fore[[iw]] = mean_function_list[[iw]] + (aggregate_fore + resi_fore[[iw]])[,fh]
    rm(iw)
  }
  return(final_fore)
}

fore_national_cdf_MLFTS <-function(data_set_F, data_set_M, aux_variable, fh, fmethod, method_ncomp)
{
  data_F = data_dum_F = data_set_F/10^5
  data_M = data_dum_M = data_set_M/10^5
  
  if(any(data_F[,1] == 0)|any(data_M[,1] == 0))
  {
    data_F[,1] = replace(data_dum_F[,1], which(data_dum_F[,1] == 0), 10^-5)
    data_M[,1] = replace(data_dum_M[,1], which(data_dum_M[,1] == 0), 10^-5)
  }
  else
  {
    data_F = data_dum_F
    data_M = data_dum_M
  }
  rm(data_dum_F); rm(data_dum_M)
  
  if(missing(aux_variable)|is.null(aux_variable))
  {
    data_cumsum_F = data_cumsum_M = matrix(NA, nrow(data_F), ncol(data_F))
    for(ij in 1:nrow(data_F))
    {
      data_cumsum_F[ij,] = cumsum(data_F[ij,])
      data_cumsum_M[ij,] = cumsum(data_M[ij,])
      rm(ij)
    }
    
    data_cumsum_logit_F = data_cumsum_logit_M = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
    for(ij in 1:nrow(data_F))
    {
      data_cumsum_logit_F[ij,] = logit(data_cumsum_F[ij, 1:(ncol(data_F) - 1)])
      data_cumsum_logit_M[ij,] = logit(data_cumsum_M[ij, 1:(ncol(data_M) - 1)])
      rm(ij)
    }
    rownames(data_cumsum_logit_F) = rownames(data_cumsum_logit_M) = years[1:nrow(data_F)]
    colnames(data_cumsum_logit_F) = colnames(data_cumsum_logit_M) = 1:(ncol(data_F) - 1)
    data_common = NULL
  }
  else
  {
    data_cumsum_F = data_cumsum_M = data_cumsum_T = matrix(NA, nrow(data_F), ncol(data_F))
    for(ij in 1:nrow(data_F))
    {
      data_cumsum_F[ij,] = cumsum(data_F[ij,])
      data_cumsum_M[ij,] = cumsum(data_M[ij,])
      data_cumsum_T[ij,] = cumsum(aux_variable[ij,])
      rm(ij)
    }
    
    data_cumsum_logit_F = data_cumsum_logit_M = data_cumsum_logit_T = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
    for(ij in 1:nrow(data_F))
    {
      data_cumsum_logit_F[ij,] = logit(data_cumsum_F[ij, 1:(ncol(data_F) - 1)])
      data_cumsum_logit_M[ij,] = logit(data_cumsum_M[ij, 1:(ncol(data_M) - 1)])
      data_cumsum_logit_T[ij,] = logit(data_cumsum_T[ij, 1:(ncol(data_M) - 1)])
      rm(ij)
    }
    rownames(data_cumsum_logit_F) = rownames(data_cumsum_logit_M) = rownames(data_cumsum_logit_T) = years[1:nrow(data_F)]
    colnames(data_cumsum_logit_F) = colnames(data_cumsum_logit_M) = colnames(data_cumsum_logit_T) = 1:(ncol(data_F) - 1)
    data_common = data_cumsum_logit_T
  }
  data_comb = array(NA, dim = c((ncol(data_F) - 1), nrow(data_F), 2))
  data_comb[,,1] = t(data_cumsum_logit_F)
  data_comb[,,2] = t(data_cumsum_logit_M)
  
  # implementing the multilevel functional data model
  
  dum = MLFTS_model(data_input = data_comb, aux_var = data_common, ncomp_method = method_ncomp,
                    fh = fh, fore_method = fmethod)
  data_cumsum_logit_F_fore = dum[[1]]
  data_cumsum_logit_M_fore = dum[[2]]
  rm(dum)
  
  data_cumsum_logit_F_fore_add = c(invlogit(data_cumsum_logit_F_fore), 1)
  data_cumsum_logit_M_fore_add = c(invlogit(data_cumsum_logit_M_fore), 1)
  
  data_cumsum_logit_F_fore_add_diff = c(data_cumsum_logit_F_fore_add[1], diff(data_cumsum_logit_F_fore_add))
  data_cumsum_logit_M_fore_add_diff = c(data_cumsum_logit_M_fore_add[1], diff(data_cumsum_logit_M_fore_add))
  return(list(mlfts_fore_F = data_cumsum_logit_F_fore_add_diff * 10^5,
              mlfts_fore_M = data_cumsum_logit_M_fore_add_diff * 10^5))
}
