####################################
# univariate functional time series 
####################################

# data_set: life-table death counts
# ncomp_method: way of selecting the number of components
# fh: forecast horizon
# fmethod: ets (we also compare auto.arima, but prefer ets)

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
    data_cumsum_logit_fore = forecast(ftsm(fts(ages[1:110], t(data_cumsum_logit)), order = ncomp), h = fh,
                                      method = fmethod)
    
    # h-step-ahead forecast
        
    data_cumsum_logit_fore_add = c(invlogit(data_cumsum_logit_fore$mean$y[,fh]), 1)
    data_cumsum_logit_fore_add_diff = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
    return(data_cumsum_logit_fore_add_diff * 10^5)
}


###########################################
# Univariate functional time series method
###########################################

# fdata: a data matrix of dimension (n by p)
# sex: female or male series
# method_ncomp: way of selecting the number of components
# horizon: forecast horizon 1 to 16
# fore_method: forecasting method
# variable_interest: point or interval forecasts (not needed)
# CLR_ncomp_selection: when the fore_method = "CLR", it requires a way for selecting number of components

point_fore_national_cdf <- function(fdata, sex, method_ncomp, horizon, fore_method, uni_fore_method,
                                    CLR_ncomp_selection)
{
    fore_val = matrix(NA, ncol(fdata), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] <- fore_national_cdf(data_set = fdata[1:(31 + ij),], ncomp_method = method_ncomp,
                                               fh = horizon, fmethod = uni_fore_method)
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] <- as.numeric(clr_fun(fdata = fdata[1:(31 + ij),], ncomp_selection = CLR_ncomp_selection,
                                                fh = horizon, fore_method = "ETS")$fore_count)
            rm(ij)
        }
    }
    else
    {
        warning("forecasting method must either be CDF or CLR.")
    }
    
    holdout_val_dum = t(matrix(fdata[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata)))
    if(any(holdout_val_dum == 0))
    {
        holdout_val = replace(x = holdout_val_dum, list = which(holdout_val_dum == 0), values = 10^-5)
    }
    else
    {
        holdout_val = holdout_val_dum
    }
    rm(holdout_val_dum)
    
    # compute the KL divergence and JS divergence
    
    KL_div_val = JS_div_val = L1_dist = L2_dist = vector("numeric", (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        # symmetric KL dist
        
        KL_div_val[ij] = mean(KLdiv(cbind(fore_val[,ij], holdout_val[,ij]))[2:3])
        
        # Jensen-Shannon dist
        
        JS_div_val[ij] = sqrt(mean(KLdiv(cbind(fore_val[,ij], 
                                    apply(cbind(fore_val[,ij], holdout_val[,ij]), 1, geometric.mean)))[2:3]))
        
        # L1 Wasserstein (not report)
         
        L1_dist[ij] = wasserstein1d(a = fore_val[,ij], b = holdout_val[,ij], p = 1)
        
        # L2 Wasserstein (not report)
        
        L2_dist[ij] = wasserstein1d(a = fore_val[,ij], b = holdout_val[,ij], p = 2)      
    }
    err = c(mean(KL_div_val), mean(JS_div_val), sqrt(mean(L1_dist^2)), sqrt(mean(L2_dist^2)))
    return(list(err = err, forecast_pdf = fore_val, holdout_pdf = holdout_val))
}

######################################
# multivariate functional time series
######################################

# data_input: multivariate functional time series
# ncomp_method: way of selecting the number of components
# fh: forecast horizon
# fore_method: forecasting method
# object_interest: point or interval forecast accuracy
# boot_number: number of bootstrap samples
# PI_level: nominal coverage probability

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
    #res_fore = as.matrix(fore_ftsm$mean$y[,fh])
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

######################################################################
# Multivariate functional time series method (evaluation of accuracy)
######################################################################

# fdata_F: a data matrix of dimension (n by p)
# fdata_M: a data matrix of dimension (n by p)
# fore_method: type of transformation
# horizon: forecast horizon 1 to 16
# way_ncomp: K = 6 or EVR
# uni_fore_method: univariate time-series forecasting method

point_fore_national_cdf_MFTS <- function(fdata_F, fdata_M, fore_method, horizon, way_ncomp, uni_fore_method)
{
    forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- fore_national_cdf_MFTS(data_set_F = fdata_F[1:(31+ij),], data_set_M = fdata_M[1:(31+ij),],
                                          fh = horizon, fmethod = uni_fore_method, object_interest = "point",
                                          alpha = 0.2, method_ncomp = way_ncomp)
            forecast_val_F[,ij] = dum$mfts_fore_F
            forecast_val_M[,ij] = dum$mfts_fore_M
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- clr_MFTS_fun(fdata_F = fdata_F[1:(31+ij),], fdata_M = fdata_M[1:(31+ij),], 
                                ncomp_selection = way_ncomp, fh = horizon, fore_method = "ets")
            forecast_val_F[,ij] = dum$MFTS_res_fore_F
            forecast_val_M[,ij] = dum$MFTS_res_fore_M
            rm(ij)
        }
    }
    
    holdout_val_dum_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_dum_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
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
        # symmetric KL dist (Female)
      
        KL_div_val_F[ij] = mean(KLdiv(cbind(forecast_val_F[,ij], holdout_val_F[,ij]))[2:3])
        
        # Jensen-Shannon dist
        
        JS_div_val_F[ij] = sqrt(mean(KLdiv(cbind(forecast_val_F[,ij], apply(cbind(forecast_val_F[,ij], holdout_val_F[,ij]), 1, geometric.mean)))[2:3]))
        
        # L1 dist
        
        L1_dist_F[ij] = wasserstein1d(a = forecast_val_F[,ij], b = holdout_val_F[,ij], p = 1)
        
        # L2 dist
        
        L2_dist_F[ij] = wasserstein1d(a = forecast_val_F[,ij], b = holdout_val_F[,ij], p = 2)
        
        # symmetric KL dist (Male)
        
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

##########################################
# Multilevel functional time series model
##########################################

# data_set: a list of p by n data matrix
# aux_var: an aggregated p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, such as "ETS"

MLFTS_model <- function(data_input, aux_var, ncomp_method, fh, fore_method)
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

## producing MLFTS point forecasts
# data_set_F: female data
# data_set_M: male data
# aux_variable: common data
# fh: forecast horizon
# fmethod: forecasting method
# method_ncomp: way of selecting the number of components

fore_national_cdf_MLFTS <- function(data_set_F, data_set_M, aux_variable, fh, fmethod, method_ncomp)
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

# fdata_F: a data matrix of dimension (n by p)
# fdata_M: a data matrix of dimension (n by p)
# fdata_common: a common data matrix
# fore_method: CDF or CLR
# horizon: forecast horizon 1 to 16
# way_ncomp: way of selecting the number of components

point_fore_national_cdf_MLFTS <- function(fdata_F, fdata_M, fdata_common, fore_method, horizon, way_ncomp)
{
    forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- fore_national_cdf_MLFTS(data_set_F = fdata_F[1:(31+ij),], data_set_M = fdata_M[1:(31+ij),],
                                           aux_variable = fdata_common[1:(31+ij),], fh = horizon, fmethod = "ets",
                                           method_ncomp = way_ncomp)
            forecast_val_F[,ij] = dum$mlfts_fore_F
            forecast_val_M[,ij] = dum$mlfts_fore_M
            rm(ij); rm(dum)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum = clr_MLFTS_fun(fdata_F = fdata_F[1:(31+ij),], fdata_M = fdata_M[1:(31+ij),],
                                ncomp_selection = way_ncomp, fh = horizon)
            forecast_val_F[,ij] = dum$MLFTS_res_fore_F
            forecast_val_M[,ij] = dum$MLFTS_res_fore_M
            rm(ij); rm(dum)
        }
    }
    rownames(forecast_val_F) = rownames(forecast_val_M) = 1:ncol(fdata_F)
    colnames(forecast_val_F) = colnames(forecast_val_M) = 1:(17 - horizon)
    
    holdout_val_dum_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_dum_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
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
      
        # L1_dist
        
        L1_dist_F[ij] = wasserstein1d(a = forecast_val_F[,ij], b = holdout_val_F[,ij], p = 1)
        
        # L2_dist
        
        L2_dist_F[ij] = wasserstein1d(a = forecast_val_F[,ij], b = holdout_val_F[,ij], p = 2)
        
        # symmetric KL dist
        
        KL_div_val_M[ij] = mean(KLdiv(cbind(forecast_val_M[,ij], holdout_val_M[,ij]))[2:3])
        
        # Jensen-Shannon dist
        
        JS_div_val_M[ij] = sqrt(mean(KLdiv(cbind(forecast_val_M[,ij], apply(cbind(forecast_val_M[,ij], holdout_val_M[,ij]), 1, geometric.mean)))[2:3]))
        
        # L1_dist
        
        L1_dist_M[ij] = wasserstein1d(a = forecast_val_M[,ij], b = holdout_val_M[,ij], p = 1)
        
        # L2_dist
        
        L2_dist_M[ij] = wasserstein1d(a = forecast_val_M[,ij], b = holdout_val_M[,ij], p = 2)
    }
    
    err_F = c(mean(KL_div_val_F), mean(JS_div_val_F), sqrt(mean(L1_dist_F^2)), sqrt(mean(L2_dist_F^2)))
    err_M = c(mean(KL_div_val_M), mean(JS_div_val_M), sqrt(mean(L1_dist_M^2)), sqrt(mean(L2_dist_M^2)))
    return(list(forecast_pdf_F = forecast_val_F, forecast_pdf_M = forecast_val_M,
                holdout_pdf_F = holdout_val_F, holdout_pdf_M = holdout_val_M,
                err_F = err_F, err_M = err_M))
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

###########################
# log-ratio approach (clr)
###########################

# fdata: n by p data matrix
# ncomp_selection: method for selecting the number of retained components
# fh: forecast horizon
# fore_method: forecasting method

clr_fun <- function(fdata, ncomp_selection, fh, fore_method)
{
    n_age = ncol(fdata)
    n_year = nrow(fdata)
    h_x_t = CLR(fdata)$LR

    SVD_decomp = svd(h_x_t)
    if(ncomp_selection == "EVR")
    {
        ncomp = select_K(tau = 0.001, eigenvalue = SVD_decomp$d^2)
    }
    else if(ncomp_selection == "fixed")
    {
        ncomp = 6
    }
    else
    {
        warning("The number of retained component must be chosen by EVR or fixed at 6.")
    }
    basis = SVD_decomp$v[,1:ncomp]
    score = t(basis) %*% t(h_x_t)
    recon = basis %*% score
    resi = t(h_x_t) - recon

    # reconstruction (model in-sample fitting)

    recon = invCLR(t(recon))

    # forecasts of principal component scores

    score_fore = matrix(NA, ncomp, 1)
    for(ik in 1:ncomp)
    {
        if(fore_method == "RWF_no_drift")
        {
            score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = FALSE)$mean[fh]
        }
        else if(fore_method == "RWF_drift")
        {
            score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = TRUE)$mean[fh]
        }
        else if(fore_method == "ETS")
        {
            score_fore[ik,] = forecast(ets(as.numeric(score[ik,])), h = fh)$mean[fh]
        }
        else if(fore_method == "ARIMA")
        {
            score_fore[ik,] = forecast(auto.arima(as.numeric(score[ik,])), h = fh)$mean[fh]
        }
        else
        {
            warning("Univariate time series forecasting method is not on the list.")
        }
    }

    # obtain forecasts in real-valued space

    fore_val = basis %*% score_fore
    fore_count = invCLR(t(fore_val)) * 10^5
    return(list(ncomp = ncomp, fore_count = fore_count))
}

## CLR multivariate functional time series
# fdata_F: female data
# fdata_M: male data
# ncomp_selection: way of selecting number of components
# fh: forecast horizon
# fore_method: forecasting method

clr_MFTS_fun <- function(fdata_F, fdata_M, ncomp_selection, fh, fore_method)
{
    n_age = ncol(fdata_F)
    n_year = nrow(fdata_F)
    
    h_x_t_F = CLR(fdata_F)$LR
    h_x_t_M = CLR(fdata_M)$LR
    
    h_x_t_comb = array(NA, dim = c(n_age, n_year, 2))
    h_x_t_comb[,,1] = t(h_x_t_F)
    h_x_t_comb[,,2] = t(h_x_t_M)
    
    n_pop = dim(h_x_t_comb)[3]
    MFTS_res = MFTS_model(data_input = h_x_t_comb, ncomp_method = ncomp_selection, fh = fh, 
                          fore_method = fore_method, object_interest = "point", PI_level = 80)
    MFTS_res_F = MFTS_res[1:n_age,]
    MFTS_res_M = MFTS_res[(n_age + 1):(n_age * 2),]
    MFTS_res_fore_F = as.numeric(invCLR(t(MFTS_res_F))) * 10^5
    MFTS_res_fore_M = as.numeric(invCLR(t(MFTS_res_M))) * 10^5
    return(list(MFTS_res_fore_F = MFTS_res_fore_F, MFTS_res_fore_M = MFTS_res_fore_M))    
}

## CLR multilevel functional time series
# fdata_F: female data
# fdata_M: male data
# ncomp_selection: way of selecting number of components
# fh: forecast horizon

clr_MLFTS_fun <- function(fdata_F, fdata_M, ncomp_selection, fh)
{
    n_age = ncol(fdata_F)
    n_year = nrow(fdata_F)
    
    h_x_t_F = CLR(fdata_F)$LR
    h_x_t_M = CLR(fdata_M)$LR
    
    h_x_t_comb = array(NA, dim = c(n_age, n_year, 2))
    h_x_t_comb[,,1] = t(h_x_t_F)
    h_x_t_comb[,,2] = t(h_x_t_M)
    
    n_pop = dim(h_x_t_comb)[3]
    MLFTS_res = MLFTS_model(data_input = h_x_t_comb, aux_var = NULL, ncomp_method = ncomp_selection, fh = fh, 
                            fore_method = "ets")
    MLFTS_res_F = MLFTS_res[[1]]
    MLFTS_res_M = MLFTS_res[[2]]
    MLFTS_res_fore_F = as.numeric(invCLR(t(MLFTS_res_F))) * 10^5
    MLFTS_res_fore_M = as.numeric(invCLR(t(MLFTS_res_M))) * 10^5
    return(list(MLFTS_res_fore_F = MLFTS_res_fore_F, MLFTS_res_fore_M = MLFTS_res_fore_M))    
}

## forecasting via high-dimensional functional principal component analysis
# object: high-dimensional functional time series

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

############################################
## multilevel functional time-series method
############################################

# fdata_F: female data
# fdata_M: male data
# fdata_common: common data
# fore_method: forecasting method
# horizon: forecast horizon
# way_ncomp: way of selecting the number of components
# level_sig: level of significance

interval_fore_national_cdf_MLFTS <- function(fdata_F, fdata_M, fdata_common, fore_method, horizon, way_ncomp, 
                                             level_sig)
{
    n_age = ncol(fdata_F)
    forecast_validation_F = forecast_validation_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- fore_national_cdf_MLFTS(data_set_F = fdata_F[1:(15+ij),], data_set_M = fdata_M[1:(15+ij),],
                                           aux_variable = fdata_common[1:(15+ij),], fh = horizon, fmethod = "ets",
                                           method_ncomp = way_ncomp)
            forecast_validation_F[,ij] = dum$mlfts_fore_F
            forecast_validation_M[,ij] = dum$mlfts_fore_M
            rm(ij); rm(dum)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum = clr_MLFTS_fun(fdata_F = fdata_F[1:(15+ij),], fdata_M = fdata_M[1:(15+ij),],
                                ncomp_selection = way_ncomp, fh = horizon)
            forecast_validation_F[,ij] = dum$MLFTS_res_fore_F
            forecast_validation_M[,ij] = dum$MLFTS_res_fore_M
            rm(ij); rm(dum)
        }
    }
    else
    {
        warning("Forecasting method must either be CDF or CLR.")
    }
    rownames(forecast_validation_F) = rownames(forecast_validation_M) = 1:ncol(fdata_F)
    colnames(forecast_validation_F) = colnames(forecast_validation_M) = 1:(17 - horizon)
    
    # holdout validation data
    
    holdout_validation_F = t(matrix(fdata_F[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_F)))
    holdout_validation_M = t(matrix(fdata_M[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_M)))
    resi_mat_F = holdout_validation_F - forecast_validation_F
    resi_mat_M = holdout_validation_M - forecast_validation_M
    
    # compute standard deviation of residuals
    
    sd_val_F = apply(resi_mat_F, 1, sd)
    sd_val_M = apply(resi_mat_M, 1, sd)
    
    # find the optimal tuning parameter
    
    tune_para_find_val_F_1 = optimise(f = tune_para_find, interval = c(0, 10),
                                      resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                      alpha_level = level_sig, PI_type = "pointwise")
    
    tune_para_find_val_F_2 = optim(par = 1, fn = tune_para_find, lower = 0, method = "L-BFGS-B",
                                   resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                   alpha_level = level_sig, PI_type = "pointwise")
    
    tune_para_find_F = ifelse(which.min(c(tune_para_find_val_F_1$objective, tune_para_find_val_F_2$value)) == 1,
                              tune_para_find_val_F_1$minimum,   tune_para_find_val_F_2$par)
    
    tune_para_find_val_M_1 = optimise(f = tune_para_find, interval = c(0, 10),
                                      resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                      alpha_level = level_sig, PI_type = "pointwise")
    
    tune_para_find_val_M_2 = optim(par = 1, fn = tune_para_find, lower = 0, method = "L-BFGS-B",
                                   resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                   alpha_level = level_sig, PI_type = "pointwise")
    
    tune_para_find_M = ifelse(which.min(c(tune_para_find_val_M_1$objective, tune_para_find_val_M_2$value)) == 1,
                              tune_para_find_val_M_1$minimum,   tune_para_find_val_M_2$par)
    
    forecast_test_F = forecast_test_F_lb = forecast_test_F_ub = matrix(NA, ncol(fdata_F), (17 - horizon))
    forecast_test_M = forecast_test_M_lb = forecast_test_M_ub = matrix(NA, ncol(fdata_M), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            dum <- fore_national_cdf_MLFTS(data_set_F = fdata_F[1:(31+ij),], data_set_M = fdata_M[1:(31+ij),],
                                           aux_variable = fdata_common[1:(31+ij),], fh = horizon, fmethod = "ets",
                                           method_ncomp = way_ncomp)
            forecast_test_F[,ij] = dum$mlfts_fore_F
            forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
            forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
            
            forecast_test_M[,ij] = dum$mlfts_fore_M
            forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
            forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
            rm(ij); rm(dum)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            dum = clr_MLFTS_fun(fdata_F = fdata_F[1:(31+ij),], fdata_M = fdata_M[1:(31+ij),],
                                ncomp_selection = way_ncomp, fh = horizon)
            forecast_test_F[,ij] = dum$MLFTS_res_fore_F
            forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
            forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
            
            forecast_test_M[,ij] = dum$MLFTS_res_fore_M
            forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
            forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
            rm(ij); rm(dum)
        }
    }
    
    # holdout testing data
    
    holdout_val_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_M)))
    
    int_F_err = interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
    int_M_err = interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
    return(list(int_F_err = int_F_err, tune_para_find_F = tune_para_find_F,
                tune_para_find_F_obj = min(c(tune_para_find_val_F_1$objective, tune_para_find_val_F_2$value)),
                int_M_err = int_M_err, tune_para_find_M = tune_para_find_M,
                tune_para_find_M_obj = min(c(tune_para_find_val_M_1$objective, tune_para_find_val_M_2$value))))
}

