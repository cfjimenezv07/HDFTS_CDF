# mean absolute percentage error

mape <- function(forecast, true)
{
  if (length(forecast) != length(true))
    stop("MAPE: the lengths of input vectors must be the same.")
  err = mean(100 * abs((true - forecast) / true))
  return(round(err, 6))
}

# root mean square percentage error

rmspe <- function(forecast, true)
{
  if (length(forecast) != length(true))
    stop("RMSPE: the lengths of input vectors must be the same.")
  err = sqrt(mean((100 * (true - forecast) / true)^2))
  return(round(err, 6))
}

# Kullback-Leibler divergence for densities

KLD <- function(forecast, true)
{
  if (length(forecast) != length(true))
    stop("KLD: the lengths of input vectors must be the same.")
  dat <- cbind(forecast, true)
  err = mean(as.numeric(KLdiv(dat, eps = 1e-16))[2:3])
  return(round(err, 6))
}

# Jenson-Shannon divergence for densities with simple mean

JSD_simple <- function(forecast, true)
{
  if (length(forecast) != length(true))
    stop("JSD_simple: the lengths of input vectors must be the same.")
  dat <- cbind(forecast, true)
  M <- rowMeans(dat)
  P_M <- cbind(dat[,1],M)
  E_M <- cbind(as.numeric(dat[,2]),M)
  colnames(E_M)=colnames(P_M)=c("True","M")
  err = as.numeric(0.5*KLdiv(P_M)+0.5*KLdiv(E_M))[3]
  return(round(err, 6))
}

# Jenson-Shannon divergence for densities with geometric mean

JSD_geom <- function(forecast, true)
{
  if (length(forecast) != length(true))
    stop("JSD_geom: the lengths of input vectors must be the same.")
  dat <- cbind(forecast, true)
  err = sqrt(mean(KLdiv(cbind(forecast, 
                              apply(dat, 1, geometric.mean)))[2:3]))
  return(round(err, 6))
}


