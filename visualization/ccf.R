acov_fun <- function (Y, nlags) 
{
    nt <- nrow(Y)
    nv <- ncol(Y)
    mean.Y <- colMeans(Y)
    fun.autocovariance <- list()
    for(kk in 0:nlags) 
    {
        fun.autocovariance[[paste("Lag", kk, sep = "")]] <- matrix(0, nrow = nv, ncol = nv)
        for (ind_curve in (1 + kk):nt) {
            fun.autocovariance[[paste("Lag", kk, sep = "")]] <- fun.autocovariance[[paste("Lag", kk, sep = "")]] + (as.numeric(Y[ind_curve - kk, ] - mean.Y)) %*% t(as.numeric(Y[ind_curve, ] - mean.Y))
        }
        fun.autocovariance[[paste("Lag", kk, sep = "")]] <- fun.autocovariance[[paste("Lag", kk, sep = "")]]/(nt - 1)
    }
    return(fun.autocovariance)
}

# X: n by p data matrix
# Y: n by p data matrix
# nlags: number of lags

cross_covariance_fun <- function(X, Y, nlags)
{
  	nt <- nrow(X)
  	nv <- ncol(X) 
  	mean.X <- colMeans(X)
  	mean.Y <- colMeans(Y)
  	fun.crosscovariance <- list()
  	for(kk in 0:nlags)
  	{
    		fun.crosscovariance[[paste("Lag", kk, sep = "")]] <- matrix(0, nrow = nv, ncol = nv)
    		for(ind_curve in (1 + kk):nt)
    		{
    			  fun.crosscovariance[[paste("Lag", kk, sep = "")]] <- fun.crosscovariance[[paste("Lag", kk, sep = "")]] + (as.numeric(X[ind_curve - kk, ] - mean.X)) %*% t(as.numeric(Y[ind_curve, ] - mean.Y))						
    		}
    		fun.crosscovariance[[paste("Lag", kk, sep = "")]] <- fun.crosscovariance[[paste("Lag", kk, sep = "")]]/(nt - 1)
  	}
  	
  	fun.autocovarianceX = fun.autocovarianceY = list()
  	for(kk in 0:nlags)
  	{
  	    fun.autocovarianceX[[paste("Lag", kk, sep = "")]] = 
  	    fun.autocovarianceY[[paste("Lag", kk, sep = "")]] = matrix(0, nrow = nv, ncol = nv)
  	    for(ind_curve in (1 + kk):nt)
  	    {
    		    fun.autocovarianceX[[paste("Lag", kk, sep = "")]] <- fun.autocovarianceX[[paste("Lag", kk, sep = "")]] + (as.numeric(X[ind_curve, ] - mean.X)) %*% t(as.numeric(X[ind_curve, ] - mean.X))/(nt - 1)						
    
    		    fun.autocovarianceY[[paste("Lag", kk, sep = "")]] <- fun.autocovarianceY[[paste("Lag", kk, sep = "")]] + (as.numeric(Y[ind_curve, ] - mean.Y)) %*% t(as.numeric(Y[ind_curve, ] - mean.Y))/(nt - 1)						
  	    }
  	}
  	return(list(cross_covariance = fun.crosscovariance, auto_covariance_X = fun.autocovarianceX, 
  	            auto_covariance_Y = fun.autocovarianceY))
}


obtain_suface_L2_norm <- function (v, autocovSurface) 
{
    repsize <- length(autocovSurface)
    matindex <- rep(NA, repsize)
    nt <- nrow(autocovSurface$Lag0)
    for(rr in 1:repsize) 
    {
        surf.aux <- autocovSurface[[paste("Lag", rr - 1, sep = "")]]^2
        norm.vec <- matrix(NA, nrow = nt, ncol = 1)
        for(ii in 1:nrow(surf.aux)) 
        {
            norm.vec[ii] <- pracma::trapz(v, surf.aux[ii, ])
        }
        norm.aux <- pracma::trapz(v, norm.vec)
        matindex[rr] <- norm.aux
    }
    return(matindex)
}

ccf_function <- function(v, cross_covariance)
{
    matindex <- obtain_suface_L2_norm(v, cross_covariance$cross_covariance)
    matindex <- matindex[-1]
    normalization.value_X <- pracma::trapz(v, diag(cross_covariance$auto_covariance_X$Lag0))
    normalization.value_Y <- pracma::trapz(v, diag(cross_covariance$auto_covariance_Y$Lag0))
    rho = sqrt(matindex)/(sqrt(normalization.value_X) * sqrt(normalization.value_Y))
    return(rho)
}

############
# read data
############

require(demography)
AUS_demo = read.demogdata("AUS_rate.txt", "AUS_expo.txt", type = "mortality", label = "AUS")
AUS_demo_smooth = smooth.demogdata(extract.ages(AUS_demo, 0:100))

AUS_demo_smooth_female = log(AUS_demo_smooth$rate$female, base = 10)
AUS_demo_smooth_male   = log(AUS_demo_smooth$rate$male, base = 10)

cross_cov = cross_covariance_fun(X = t(AUS_demo_smooth_female), Y = t(AUS_demo_smooth_male), nlags = 20)
ccf_value = ccf_function(v = 0:100, cross_covariance = cross_cov)

plot(1:20, ccf_value, type = "h", xlab = "Lag", ylab = "Cross correlation function", ylim = c(0, 1))

# check (okay)

cross_cov_test = cross_covariance_fun(X = t(AUS_demo_smooth_male), Y = t(AUS_demo_smooth_male), nlags = 20)
ccf_value_test = ccf_function(v = 0:100, cross_covariance = cross_cov_test)

ccf_value_female_prefecture = ccf_value_male_prefecture = ccf_value_total_prefecture = matrix(NA, 20, 47)
for(iw in 1:47)
{
    ccf_value_female_prefecture[,iw] = ccf_function(v = 0:110, cross_covariance = 
               cross_covariance_fun(X = female_prefecture_dx[[iw]], Y = Japan_female_pop, nlags = 20))
    
    ccf_value_male_prefecture[,iw] = ccf_function(v = 0:110, cross_covariance = 
               cross_covariance_fun(X = male_prefecture_dx[[iw]], Y = Japan_male_pop, nlags = 20))
    
    ccf_value_total_prefecture[,iw] = ccf_function(v = 0:110, cross_covariance = 
               cross_covariance_fun(X = total_prefecture_dx[[iw]], Y = Japan_total_pop, nlags = 20))
    print(iw); rm(iw)
}
rownames(ccf_value_female_prefecture) = rownames(ccf_value_male_prefecture) = rownames(ccf_value_total_prefecture) = 1:20  
colnames(ccf_value_female_prefecture) = colnames(ccf_value_male_prefecture) = colnames(ccf_value_total_prefecture) = state
             
savefig("CCF_female", width = 12, height = 10, toplines = 0.8, type = "png")
image(1:20, 1:47, ccf_value_female_prefecture[,rev(1:47)], ylab = "Prefecture", xlab = "Lag", 
      yaxt = "n", main = "Female")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig("CCF_male", width = 12, height = 10, toplines = 0.8, type = "png")
image(1:20, 1:47, ccf_value_male_prefecture[,rev(1:47)], ylab = "", xlab = "Lag", yaxt = "n", main = "Male")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig("Fig_3c", width = 12, height = 10, toplines = 0.8, type = "png")
#image(1:20, 1:47, ccf_value_total_prefecture[,rev(1:47)], ylab = "", xlab = "Lag", yaxt = "n", main = "Total")
#axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))

filled.contour(1:20, 1:47, ccf_value_total_prefecture[,rev(1:47)], ylab = "", xlab = "Lag", 
               yaxt = "n", main = "Total", plot.axes = {axis(1)})
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()
