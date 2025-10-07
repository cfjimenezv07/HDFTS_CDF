#CDF transformation

# data: scale to between 0 and 1, a data matrix of dimension (n by p) where n denotes the number of years and p denotes the number of ages

# ncomp_method: method for selecting the number of components

# fh: forecast horizon

# fmethod: forecasting method

# object_interest: point or interval forecasts

# no_boot: number of bootstrap samples

# alpha: level of significance between 0 and 1



cdf_transformation <- function(data,year)
  
{
  
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
  
  rownames(data_cumsum_logit) = year[1:nrow(data)]
  

  
  return(CDF_transformed_data= data_cumsum_logit)
}