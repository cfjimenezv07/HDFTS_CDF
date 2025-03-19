# convert an array to a matrix

data_pop_F_mat = data_pop_M_mat = matrix(NA, 111, 48 * 47)
for(ij in 1:47)
{
    data_pop_F_mat[,(48 * (ij-1) + 1):(48 * ij)] = t(female_prefecture_dx[[ij]])
    data_pop_M_mat[,(48 * (ij-1) + 1):(48 * ij)] = t(male_prefecture_dx[[ij]])
    rm(ij)
}

# extract functional grand mean, functional row mean, functional column mean

dum = FANOVA(data_pop1 = data_pop_F_mat, data_pop2 = data_pop_M_mat, year = 1975:2022, age = 0:110,
       n_prefectures = 47, n_populations = 2)

# main function to extract the functional residuals

dum_2 = Two_way_Residuals_means(data_pop1 = data_pop_F_mat, data_pop2 = data_pop_M_mat, year = 1975:2022, 
                              age = 0:110, n_prefectures = 47, n_populations = 2)

test_resi = matrix(NA, 48, 111)
for(iw in 1:48)
{
    test_resi[iw,] = (female_prefecture_dx[[1]])[iw,] - (dum$FGE_mean + dum$FCE_mean[1,] + dum$FRE_mean[1,])
}
  
all(round(dum_2$residuals1_mean[1:48,], 8) == round(test_resi, 8)) # TRUE

# extract functional residuals

data_pop_F_mat_resi = data_pop_M_mat_resi = array(NA, dim = c(48, 111, 47))
for(ik in 1:47)
{
    data_pop_F_mat_resi[,,ik] = dum_2$residuals1_mean[((ik-1)*48+1):(ik*48),]
    data_pop_M_mat_resi[,,ik] = dum_2$residuals2_mean[((ik-1)*48+1):(ik*48),]
}

# female data array (normalize to sum to one)
# male data array (normalize to sum to one)

female_prefecture_dx_array = male_prefecture_dx_array = array(NA, dim = c(111, 48, 47))
for(ij in 1:47)
{
    female_prefecture_dx_array[,,ij] = t(female_prefecture_dx[[ij]])/10^5
    male_prefecture_dx_array[,,ij]   = t(male_prefecture_dx[[ij]])/10^5
    rm(ij)
}

# data_F_origin: female data series
# data_M_origin: male data series
# way_ncomp: way of selecting ncomp
# horizon: forecast horizon

fore_CDF_FANOVA <- function(data_F_origin, data_M_origin, way_ncomp, horizon)
{
    # dimension
  
    n_age = dim(data_F_origin)[1]
    n_year = dim(data_F_origin)[2]
    n_state = dim(data_F_origin)[3]
      
    # cumulative sum data_pop_F_array should be [0,1]
    
    data_cumsum_dum_F = data_cumsum_dum_M = array(NA, dim = c(n_year, n_age, n_state))
    for(ik in 1:n_state)
    {
        for(ij in 1:n_year)
        {
            data_cumsum_dum_F[ij,,ik] = cumsum(data_pop_F_array[,ij,ik])
            data_cumsum_dum_M[ij,,ik] = cumsum(data_pop_M_array[,ij,ik])
            rm(ij)
        }
        rm(ik)
    }
    
    # logit can't handle zero at the beginning age
    
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
    
    # take logit transformation
    
    data_cumsum_logit_F = data_cumsum_logit_M = array(NA, dim = c(n_year, (n_age - 1), n_state))
    for(ik in 1:n_state)
    {
        for(ij in 1:n_year)
        {
            data_cumsum_logit_F[ij,,ik] = logit(data_cumsum_F[ij, 1:(n_age - 1), ik])
            data_cumsum_logit_M[ij,,ik] = logit(data_cumsum_M[ij, 1:(n_age - 1), ik])
            rm(ij)
        }
        rm(ik)
    }

    # array and matrix formats
    
    data_pop_F_logit_array = data_pop_M_logit_array = array(NA, dim = c((n_age - 1), n_year, n_state))
    data_pop_F_mat = data_pop_M_mat = matrix(NA, (n_age - 1), n_year * n_state)
    for(ij in 1:n_state)
    {
        data_pop_F_logit_array[,,ij] = data_pop_F_mat[,(n_year * (ij-1) + 1):(n_year * ij)] = t(data_cumsum_logit_F[,,ij])
        data_pop_M_logit_array[,,ij] = data_pop_M_mat[,(n_year * (ij-1) + 1):(n_year * ij)] = t(data_cumsum_logit_M[,,ij])
        rm(ij)
    }

    # FANOVA
    
    FANOVA_means = FANOVA(data_pop1 = data_pop_F_mat, data_pop2 = data_pop_M_mat, year = (1975:2022)[1:n_year], 
                          age = 0:109, n_prefectures = n_state, n_populations = 2)
        
    # means
    
    FANOVA_grand_mean = FANOVA_means$FGE_mean
    FANOVA_state_mean = FANOVA_means$FRE_mean
    FANOVA_sex_mean = FANOVA_means$FCE_mean
    
    # FANOVA residual functions
    
    FANOVA_resi = Two_way_Residuals_means(data_pop1 = data_pop_F_mat, data_pop2 = data_pop_M_mat,
                                     year = (1975:2022)[1:n_year], age = 0:109, n_prefectures = n_state, 
                                     n_populations = 2)

    # residual functions in an array format
    
    data_pop_F_mat_resi = data_pop_M_mat_resi = array(NA, dim = c(n_year, (n_age - 1), n_state))
    for(ik in 1:n_state)
    {
        data_pop_F_mat_resi[,,ik] = FANOVA_resi$residuals1_mean[((ik-1) * n_year + 1):(ik * n_year),]
        data_pop_M_mat_resi[,,ik] = FANOVA_resi$residuals2_mean[((ik-1) * n_year + 1):(ik * n_year),]
        rm(ik)
    }
    
    # curve reconstruction 
    
    recon_F_array = recon_M_array = array(NA, dim = c(n_year, (n_age - 1), n_state))
    for(ik in 1:n_state)
    {
        recon_F_array[,,ik] <- data_pop_F_mat_resi[,,ik] + matrix(rep(FANOVA_sex_mean[1,], n_year), n_year, (n_age - 1), byrow = TRUE) + 
                              matrix(rep(FANOVA_state_mean[ik,], n_year), n_year, (n_age - 1), byrow = TRUE) + 
                              matrix(rep(FANOVA_grand_mean, n_year), n_year, (n_age - 1), byrow = TRUE)
        
        recon_M_array[,,ik] <- data_pop_M_mat_resi[,,ik] + matrix(rep(FANOVA_sex_mean[2,], n_year), n_year, (n_age - 1), byrow = TRUE) + 
          matrix(rep(FANOVA_state_mean[ik,], n_year), n_year, (n_age - 1), byrow = TRUE) + 
          matrix(rep(FANOVA_grand_mean, n_year), n_year, (n_age - 1), byrow = TRUE)
        rm(ik)
    } 
    
    # exact recovery
    
    recon_F_ind = recon_M_ind = vector(, n_state)
    for(ik in 1:n_state)
    {
        recon_F_ind[ik] = all(round(data_pop_F_logit_array[,,ik], 6) == round(t(recon_F_array[,,ik]), 6))
        recon_M_ind[ik] = all(round(data_pop_M_logit_array[,,ik], 6) == round(t(recon_M_array[,,ik]), 6))
        rm(ik)
    }
    all(recon_F_ind)
    all(recon_M_ind)
    
    # forecast residual functions by multivariate FTS
    
    MFTS_F_output = MFTS_M_output = matrix(NA, (n_age - 1), n_state)
    for(iw in 1:n_state)
    {
        MFTS_data_input = array(NA, dim = c((n_age - 1), n_year, 2))
        MFTS_data_input[,,1] = t(data_pop_F_mat_resi[,,iw])
        MFTS_data_input[,,2] = t(data_pop_M_mat_resi[,,iw])
        MFTS_output = MFTS_model(data_input = MFTS_data_input, ncomp_method = way_ncomp, fh = horizon, 
                                 fore_method = "ets")
        MFTS_F_output[,iw] = MFTS_output[1:(n_age - 1)]
        MFTS_M_output[,iw] = MFTS_output[n_age:((n_age - 1) * 2)]
        rm(iw)
    }        
    
    # h-step-ahead forecast of the residual functions
    # add the fixed effect, i.e., those means
    
    FANOVA_logit_fore_F = FANOVA_logit_fore_M = matrix(NA, (n_age - 1), n_state)
    for(iw in 1:n_state)
    {
        FANOVA_logit_fore_F[,iw] = FANOVA_grand_mean + FANOVA_sex_mean[1,] + FANOVA_state_mean[iw,] + MFTS_F_output[,iw]
        FANOVA_logit_fore_M[,iw] = FANOVA_grand_mean + FANOVA_sex_mean[2,] + FANOVA_state_mean[iw,] + MFTS_M_output[,iw]
        rm(iw)
    }
    
    # inverse logit transformation
    
    data_cumsum_logit_fore_add_diff_F = data_cumsum_logit_fore_add_diff_M = matrix(NA, n_age, n_state) 
    for(iw in 1:n_state)
    {
        data_cumsum_logit_fore_add = c(invlogit(FANOVA_logit_fore_F[,iw]), 1)
        data_cumsum_logit_fore_add_diff_F[,iw] = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
        rm(data_cumsum_logit_fore_add)
        
        data_cumsum_logit_fore_add = c(invlogit(FANOVA_logit_fore_M[,iw]), 1)
        data_cumsum_logit_fore_add_diff_M[,iw] = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
        rm(data_cumsum_logit_fore_add); rm(iw)
    }
    colnames(data_cumsum_logit_fore_add_diff_F) = colnames(data_cumsum_logit_fore_add_diff_M) = state
    
    return(list(FANOVA_means_F = data_cumsum_logit_fore_add_diff_F * 10^5, 
                FANOVA_means_M = data_cumsum_logit_fore_add_diff_M * 10^5))
}

# an example

dum = fore_CDF_FANOVA(data_F_origin = female_prefecture_dx_array[,1:32,], 
                      data_M_origin = male_prefecture_dx_array[,1:32,], 
                      way_ncomp = "EVR", horizon = 1)

#############
# evaluation 
#############

# method_ncomp: EVR or provide (K = 6)
# horizon: forecast horizon from h = 1 to 16

fore_CDF_FANOVA_eval <- function(method_ncomp, horizon)
{
    n_age = dim(female_prefecture_dx_array)[1]
    n_state = dim(female_prefecture_dx_array)[3]
    
    fore_val_F = fore_val_M = array(NA, dim = c(n_age, (17 - horizon), n_state))
    for(ij in 1:(17 - horizon))
    {
        dum_output <- fore_CDF_FANOVA(data_F_origin = female_prefecture_dx_array[,1:(31 + ij),], 
                                      data_M_origin = male_prefecture_dx_array[,1:(31 + ij),], 
                                      way_ncomp = method_ncomp, horizon = horizon)
        fore_val_F[,ij,] = dum_output$FANOVA_means_F
        fore_val_M[,ij,] = dum_output$FANOVA_means_M
        rm(ij)
    }
    
    if(horizon == 16)
    {
        holdout_F = array(female_prefecture_dx_array[,48,], dim = c(n_age, 1, n_state))
        holdout_M = array(male_prefecture_dx_array[,48,], dim = c(n_age, 1, n_state))
    }
    else
    {
        holdout_F = female_prefecture_dx_array[,(32 + horizon):48,] * 10^5
        holdout_M = male_prefecture_dx_array[,(32 + horizon):48,] * 10^5
    }
    
    KL_div_val_F = JS_div_val_F = KL_div_val_M = JS_div_val_M = matrix(NA, (17 - horizon), n_state)
    for(iw in 1:n_state)
    {
        for(ij in 1:(17 - horizon))
        {
            # symmetric KL dist
          
            KL_div_val_F[ij,iw] = mean(KLdiv(cbind(fore_val_F[,ij,iw], holdout_F[,ij,iw]))[2:3])
            KL_div_val_M[ij,iw] = mean(KLdiv(cbind(fore_val_M[,ij,iw], holdout_M[,ij,iw]))[2:3])
            
            # Jensen-Shannon dist
            
            JS_div_val_F[ij,iw] = sqrt(mean(KLdiv(cbind(fore_val_F[,ij,iw], 
                          apply(cbind(fore_val_F[,ij,iw], holdout_F[,ij,iw]), 1, geometric.mean)))[2:3]))
            
            JS_div_val_M[ij,iw] = sqrt(mean(KLdiv(cbind(fore_val_M[,ij,iw], 
                          apply(cbind(fore_val_M[,ij,iw], holdout_M[,ij,iw]), 1, geometric.mean)))[2:3]))
        }
    }
    
    output_prefecture = cbind(colMeans(KL_div_val_F), colMeans(JS_div_val_F), 
                              colMeans(KL_div_val_M), colMeans(JS_div_val_M))
    return(output_prefecture)
}

# EVR

output_val = array(NA, dim = c(16, 47, 4), dimnames = list(1:16, 1:47, c("KLD (F)", "JSD (F)", "KLD (M)", "JSD (M)")))
for(iwk in 1:16)
{
    output_val[iwk,,] = fore_CDF_FANOVA_eval(method_ncomp = "EVR", horizon = iwk)
    print(iwk); rm(iwk)
}

output_val_overall = rbind(apply(output_val, c(1, 3), mean),
                            colMeans(apply(output_val, c(1, 3), mean)),
                            apply(apply(output_val, c(1, 3), mean), 2, median))
rownames(output_val_overall) = c(1:16, "Mean", "Median")

# K = 6

output_val_K6 = array(NA, dim = c(16, 47, 4), dimnames = list(1:16, 1:47, c("KLD (F)", "JSD (F)", "KLD (M)", "JSD (M)")))
for(iwk in 1:16)
{
    output_val_K6[iwk,,] = fore_CDF_FANOVA_eval(method_ncomp = "provide", horizon = iwk)
    print(iwk); rm(iwk)
}

output_val_K6_overall = rbind(apply(output_val_K6, c(1, 3), mean),
                           colMeans(apply(output_val_K6, c(1, 3), mean)),
                           apply(apply(output_val_K6, c(1, 3), mean), 2, median))
rownames(output_val_K6_overall) = c(1:16, "Mean", "Median")

