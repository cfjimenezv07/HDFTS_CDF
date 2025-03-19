# fdata_F: female data
# fdata_M: male data
# horizon: forecast horizon
# first_order: number of components
# second_order: number of components
# transformation: type of transformation
# level_sig: level of significance

int_hdfpca_fun <- function(fdata_F, fdata_M, horizon, first_order, second_order, transformation, level_sig)
{
  	n_age = ncol(fdata_F)
  	forecast_validation_F = forecast_validation_M = matrix(NA, ncol(fdata_F), (17 - horizon))
  	if(transformation == "direct")
  	{
    		for(ij in 1:(17 - horizon))
    		{
      			data_F = fdata_F[1:(15+ij),]
      			data_M = fdata_M[1:(15+ij),]
      			data_comb = list()
      			data_comb[[1]] = t(data_F)
      			data_comb[[2]] = t(data_M)
      			
      			fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
            forecast_validation_F[,ij] = (fore_val[[1]])[,horizon]
            forecast_validation_M[,ij] = (fore_val[[2]])[,horizon]
            rm(ij); rm(fore_val)
    		}
  	}
  	else if(transformation == "CLR")
  	{
    		for(ij in 1:(17 - horizon))
    		{
      			data_F = as.matrix(clr(fdata_F[1:(15+ij),]))
            data_M = as.matrix(clr(fdata_M[1:(15+ij),]))
            data_comb = list()
            data_comb[[1]] = t(data_F)
            data_comb[[2]] = t(data_M)
                
            fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
            forecast_validation_F[,ij] = as.numeric(clrInv((fore_val[[1]])[,horizon])) * 10^5
            forecast_validation_M[,ij] = as.numeric(clrInv((fore_val[[2]])[,horizon])) * 10^5
            rm(ij); rm(fore_val)
    		}	
  	}
  	else if(transformation == "CDF")
    {
          for(ijk in 1:(17 - horizon))
          {
              data_F = fdata_F[1:(15+ijk),]/10^5
              data_M = fdata_M[1:(15+ijk),]/10^5
              
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
              
              forecast_validation_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
              forecast_validation_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
              rm(ijk); rm(data_F); rm(data_M); rm(fore_val)
          }
      }
      else
      {
      	  warning("none, CLR, CDF transformation allowed only.")
      }
      
      # holdout validation data
      
      holdout_validation_F = t(matrix(fdata_F[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_F)))
      holdout_validation_M = t(matrix(fdata_M[(16 + horizon):32,], length((16 + horizon):32), ncol(fdata_M)))
      resi_mat_F = holdout_validation_F - forecast_validation_F
      resi_mat_M = holdout_validation_M - forecast_validation_M
      
   	# compute standard deviation of residuals
          
      sd_val_F = apply(resi_mat_F, 1, sd)
      sd_val_M = apply(resi_mat_M, 1, sd)
      
      # find the optimal tuning parameter
                          
      tune_para_find_val_F_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                      resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                      alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_F_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                        resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                        alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_F_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                        resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                        alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_F_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                        resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                        alpha_level = level_sig, PI_type = "pointwise")
          
      tune_para_find_val_F_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                   resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                   alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_F_6 = optim(par = 1, fn = tune_para_find_function, 
                                     resi_mat = resi_mat_F, sd_val_input = sd_val_F,
                                     alpha_level = level_sig, PI_type = "pointwise")

      obj_val = c(tune_para_find_val_F_1$objective, 
                  tune_para_find_val_F_2$objective, 
                  tune_para_find_val_F_3$objective, 
                  tune_para_find_val_F_4$objective, 
                  tune_para_find_val_F_5$value, 
                  tune_para_find_val_F_6$value)
      obj_val_min_F = min(obj_val)
      
      tune_para_find_F = c(tune_para_find_val_F_1$minimum, 
                           tune_para_find_val_F_2$minimum, 
                           tune_para_find_val_F_3$minimum, 
                           tune_para_find_val_F_4$minimum, 
                           tune_para_find_val_F_5$par,
                           tune_para_find_val_F_6$par)[which.min(obj_val)]
      rm(obj_val)
      
      # male
      
      tune_para_find_val_M_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                        resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                        alpha_level = level_sig, PI_type = "pointwise")
                                                
      tune_para_find_val_M_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                        resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                        alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_M_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                        resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                        alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_M_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                        resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                        alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_M_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                     resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                     alpha_level = level_sig, PI_type = "pointwise")
      
      tune_para_find_val_M_6 = optim(par = 1, fn = tune_para_find_function, 
                                     resi_mat = resi_mat_M, sd_val_input = sd_val_M,
                                     alpha_level = level_sig, PI_type = "pointwise")
      
      obj_val = c(tune_para_find_val_M_1$objective, 
                  tune_para_find_val_M_2$objective, 
                  tune_para_find_val_M_3$objective, 
                  tune_para_find_val_M_4$objective, 
                  tune_para_find_val_M_5$value, 
                  tune_para_find_val_M_6$value)
      obj_val_min_M = min(obj_val)
      
      tune_para_find_M = c(tune_para_find_val_M_1$minimum, 
                           tune_para_find_val_M_2$minimum, 
                           tune_para_find_val_M_3$minimum, 
                           tune_para_find_val_M_4$minimum, 
                           tune_para_find_val_M_5$par,
                           tune_para_find_val_M_6$par)[which.min(obj_val)]
      rm(obj_val)
      
      forecast_test_F = forecast_test_F_lb = forecast_test_F_ub = matrix(NA, ncol(fdata_F), (17 - horizon))
      forecast_test_M = forecast_test_M_lb = forecast_test_M_ub = matrix(NA, ncol(fdata_M), (17 - horizon))
      if(transformation == "direct")
  	  {
  		    for(ij in 1:(17 - horizon))
  		    {
        			data_F = fdata_F[1:(31+ij),]
        			data_M = fdata_M[1:(31+ij),]
        			data_comb = list()
        			data_comb[[1]] = t(data_F)
        			data_comb[[2]] = t(data_M)
        			
  			      fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
              forecast_test_F[,ij] = (fore_val[[1]])[,horizon]
              forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
              forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
              
              forecast_test_M[,ij] = (fore_val[[2]])[,horizon]
              forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
              forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
  			      rm(ij)
  		    }		
  	  }
  	else if(transformation == "CLR")
  	{
    		for(ij in 1:(17 - horizon))
    		{
  	    		data_F = as.matrix(clr(fdata_F[1:(31+ij),]))
            data_M = as.matrix(clr(fdata_M[1:(31+ij),]))
            data_comb = list()
            data_comb[[1]] = t(data_F)
            data_comb[[2]] = t(data_M)
              
            fore_val = forecast.hdfpca(hdfpca(y = data_comb, order = first_order, r = second_order), h = horizon)$forecast
            forecast_test_F[,ij] = as.numeric(clrInv((fore_val[[1]])[,horizon])) * 10^5
            forecast_test_F_lb[,ij] = forecast_test_F[,ij] - tune_para_find_F * sd_val_F
            forecast_test_F_ub[,ij] = forecast_test_F[,ij] + tune_para_find_F * sd_val_F
  
            forecast_test_M[,ij] = as.numeric(clrInv((fore_val[[2]])[,horizon])) * 10^5
            forecast_test_M_lb[,ij] = forecast_test_M[,ij] - tune_para_find_M * sd_val_M
            forecast_test_M_ub[,ij] = forecast_test_M[,ij] + tune_para_find_M * sd_val_M
            rm(ij); rm(fore_val)
  		  }	
  	}
  	else if(transformation == "CDF")
    {
        for(ijk in 1:(17 - horizon))
        {
            data_F = fdata_F[1:(31+ijk),]/10^5
            data_M = fdata_M[1:(31+ijk),]/10^5
              
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
              
              forecast_test_F[,ijk] = data_cumsum_logit_fore_add_diff_F * 10^5
              forecast_test_F_lb[,ijk] = forecast_test_F[,ijk] - tune_para_find_F * sd_val_F
              forecast_test_F_ub[,ijk] = forecast_test_F[,ijk] + tune_para_find_F * sd_val_F
  
              forecast_test_M[,ijk] = data_cumsum_logit_fore_add_diff_M * 10^5
              forecast_test_M_lb[,ijk] = forecast_test_M[,ijk] - tune_para_find_M * sd_val_M
              forecast_test_M_ub[,ijk] = forecast_test_M[,ijk] + tune_para_find_M * sd_val_M
              rm(ijk); rm(data_F); rm(data_M)
          }
      }	
  	
  	  # holdout testing data
          
      holdout_val_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
      holdout_val_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_M)))
  
      int_F_err = interval_score(holdout = holdout_val_F, lb = forecast_test_F_lb, ub = forecast_test_F_ub, alpha = (1 - level_sig))
      int_M_err = interval_score(holdout = holdout_val_M, lb = forecast_test_M_lb, ub = forecast_test_M_ub, alpha = (1 - level_sig))
      return(list(int_F_err = int_F_err, tune_para_find_F = tune_para_find_F, 
                  tune_para_find_F_obj = obj_val_min_F,
                  
                  int_M_err = int_M_err, tune_para_find_M = tune_para_find_M,
                  tune_para_find_M_obj = obj_val_min_M))    
}

###############################
## level of significance = 0.8
###############################

# CDF

hdfpca_int_fore_subnational_err_F_EVR_ETS = hdfpca_int_fore_subnational_err_M_EVR_ETS = array(NA, dim = c(47, 15, 3), dimnames = list(state, 1:15, c("ECP", "CPD", "score")))
hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para = hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj = 
hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para = hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj = matrix(NA, 47, 15)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = int_hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], 
                             horizon = iw, first_order = 6, second_order = 2, transformation = "CDF", 
                             level_sig = 0.8)
        hdfpca_int_fore_subnational_err_F_EVR_ETS[ij,iw,] = dum$int_F_err
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para[ij,iw] = dum$tune_para_find_F
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj[ij,iw] = dum$tune_para_find_F_obj
        
        hdfpca_int_fore_subnational_err_M_EVR_ETS[ij,iw,] = dum$int_M_err
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para[ij,iw] = dum$tune_para_find_M
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj[ij,iw] = dum$tune_para_find_M_obj
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS = apply(hdfpca_int_fore_subnational_err_F_EVR_ETS, c(2, 3), mean)
horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS = apply(hdfpca_int_fore_subnational_err_M_EVR_ETS, c(2, 3), mean)

# CLR

hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR = hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR = array(NA, dim = c(47, 15, 3), dimnames = list(state, 1:15, c("ECP", "CPD", "score")))
hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_CLR = hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_CLR = 
hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_CLR = hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_CLR = matrix(NA, 47, 15)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = int_hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], 
                             horizon = iw, first_order = 6, second_order = 2, transformation = "CLR", 
                             level_sig = 0.8)
        hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR[ij,iw,] = dum$int_F_err
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_CLR[ij,iw] = dum$tune_para_find_F
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_CLR[ij,iw] = dum$tune_para_find_F_obj
        
        hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR[ij,iw,] = dum$int_M_err
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_CLR[ij,iw] = dum$tune_para_find_M
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_CLR[ij,iw] = dum$tune_para_find_M_obj
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR = apply(hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR, c(2, 3), mean)
horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR = apply(hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR, c(2, 3), mean)

# direct

hdfpca_int_fore_subnational_err_F_EVR_ETS_direct = hdfpca_int_fore_subnational_err_M_EVR_ETS_direct = array(NA, dim = c(47, 15, 3), dimnames = list(state, 1:15, c("ECP", "CPD", "score")))
hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_direct = hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_direct = 
hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_direct = hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_direct = matrix(NA, 47, 15)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = int_hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], 
                             horizon = iw, first_order = 6, second_order = 2, transformation = "direct", 
                             level_sig = 0.8)
        hdfpca_int_fore_subnational_err_F_EVR_ETS_direct[ij,iw,] = dum$int_F_err
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_direct[ij,iw] = dum$tune_para_find_F
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_direct[ij,iw] = dum$tune_para_find_F_obj
        
        hdfpca_int_fore_subnational_err_M_EVR_ETS_direct[ij,iw,] = dum$int_M_err
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_direct[ij,iw] = dum$tune_para_find_M
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_direct[ij,iw] = dum$tune_para_find_M_obj
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_direct = apply(hdfpca_int_fore_subnational_err_F_EVR_ETS_direct, c(2, 3), mean)
horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_direct = apply(hdfpca_int_fore_subnational_err_M_EVR_ETS_direct, c(2, 3), mean)


################################
## level of significance = 0.95
################################

# CDF

hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95 = array(NA, dim = c(47, 15, 3), dimnames = list(state, 1:15, c("ECP", "CPD", "score")))
hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95 = hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95 = 
hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95 = matrix(NA, 47, 15)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = int_hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], 
                             horizon = iw, first_order = 6, second_order = 2, transformation = "CDF", 
                             level_sig = 0.95)
        hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95[ij,iw,] = dum$int_F_err
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95[ij,iw] = dum$tune_para_find_F
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95[ij,iw] = dum$tune_para_find_F_obj
        
        hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95[ij,iw,] = dum$int_M_err
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95[ij,iw] = dum$tune_para_find_M
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95[ij,iw] = dum$tune_para_find_M_obj
        rm(dum); print(iw); rm(iw)
    }
    print(ij); rm(ij)
}

A=c("hdfpca_int_fore_subnational_err_F_EVR_ETS",
"hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para",
"hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj",

"hdfpca_int_fore_subnational_err_M_EVR_ETS",
"hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para",
"hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj",

"hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95",
"hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_alpha_0.95",
"hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_alpha_0.95",

"hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95",
"hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_alpha_0.95",
"hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_alpha_0.95")


horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95 = apply(hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95, c(2, 3), mean)
horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95 = apply(hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95, c(2, 3), mean)

# CLR

hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR_alpha_0.95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR_alpha_0.95 = array(NA, dim = c(47, 15, 3), dimnames = list(state, 1:15, c("ECP", "CPD", "score")))
hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_CLR_alpha_0.95 = hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_CLR_alpha_0.95 = 
hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_CLR_alpha_0.95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_CLR_alpha_0.95 = matrix(NA, 47, 15)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = int_hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], 
                             horizon = iw, first_order = 6, second_order = 2, transformation = "CLR", 
                             level_sig = 0.95)
        hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR_alpha_0.95[ij,iw,] = dum$int_F_err
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_CLR_alpha_0.95[ij,iw] = dum$tune_para_find_F
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_CLR_alpha_0.95[ij,iw] = dum$tune_para_find_F_obj
        
        hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR_alpha_0.95[ij,iw,] = dum$int_M_err
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_CLR_alpha_0.95[ij,iw] = dum$tune_para_find_M
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_CLR_alpha_0.95[ij,iw] = dum$tune_para_find_M_obj
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR_alpha_0.95 = apply(hdfpca_int_fore_subnational_err_F_EVR_ETS_CLR_alpha_0.95, c(2, 3), mean)
horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR_alpha_0.95 = apply(hdfpca_int_fore_subnational_err_M_EVR_ETS_CLR_alpha_0.95, c(2, 3), mean)

# direct

hdfpca_int_fore_subnational_err_F_EVR_ETS_direct_alpha_0.95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_direct_alpha_0.95 = array(NA, dim = c(47, 15, 3), dimnames = list(state, 1:15, c("ECP", "CPD", "score")))
hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_direct_alpha_0.95 = hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_direct_alpha_0.95 = 
hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_direct_alpha_0.95 = hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_direct_alpha_0.95 = matrix(NA, 47, 15)
for(ij in 1:47)
{
    for(iw in 1:15)
    {
        dum = int_hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], 
                             horizon = iw, first_order = 6, second_order = 2, transformation = "direct", 
                             level_sig = 0.95)
        hdfpca_int_fore_subnational_err_F_EVR_ETS_direct_alpha_0.95[ij,iw,] = dum$int_F_err
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_direct_alpha_0.95[ij,iw] = dum$tune_para_find_F
        hdfpca_int_fore_subnational_err_F_EVR_ETS_tune_para_obj_direct_alpha_0.95[ij,iw] = dum$tune_para_find_F_obj
        
        hdfpca_int_fore_subnational_err_M_EVR_ETS_direct_alpha_0.95[ij,iw,] = dum$int_M_err
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_direct_alpha_0.95[ij,iw] = dum$tune_para_find_M
        hdfpca_int_fore_subnational_err_M_EVR_ETS_tune_para_obj_direct_alpha_0.95[ij,iw] = dum$tune_para_find_M_obj
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_direct_alpha_0.95 = apply(hdfpca_int_fore_subnational_err_F_EVR_ETS_direct_alpha_0.95, c(2, 3), mean)
horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_direct_alpha_0.95 = apply(hdfpca_int_fore_subnational_err_M_EVR_ETS_direct_alpha_0.95, c(2, 3), mean)

