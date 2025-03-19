##################################################################
# determining the optimal tuning parameter via standard deviation
##################################################################

# tune_para: tuning parameter
# resi_mat: residual functions
# sd_val_input: functional standard deviation (pointwise)
# alpha_level: level of significance
# PI_type: pointwise or uniform

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
                                    h = fh, method = fmethod, pimethod = "nonparametric", B = 399)
  
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

hdfpca_fun_fore <- function(fdata_F, fdata_M, horizon, first_order, second_order, transformation)
{
    forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    if(transformation == "direct")
    {
        for(ij in 1:(17 - horizon))
        {
            data_F = fdata_F[1:(nrow(fdata_F)-17+ij),]
            data_M = fdata_M[1:(nrow(fdata_F)-17+ij),]
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
            
            forecast_val_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
            forecast_val_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
            rm(ijk); rm(data_F); rm(data_M)
        }
    }
    return(list(forecast_val_F = forecast_val_F, forecast_val_M = forecast_val_M))
}

