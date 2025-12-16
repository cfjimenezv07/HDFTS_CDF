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
# Where heatmaps will be saved
dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/Plots/"
dir.r1 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/FANOVA/"

require(demography)
female_prefecture_dx <- readRDS(paste0(dir.r1, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(paste0(dir.r1, "male_prefecture_dx.rds"))
total_prefecture_dx  <- readRDS(paste0(dir.r1, "total_prefecture_dx.rds"))
Japan_female_pop <- readRDS(paste0(dir.r1, "Japan_female_pop.rds"))
Japan_male_pop   <- readRDS(paste0(dir.r1, "Japan_male_pop.rds"))
Japan_total_pop   <- readRDS(paste0(dir.r1, "Japan_total_pop.rds"))

ages <- 0:110
n_age <- length(ages)
year <- 1975:2023
n_year <- length(year)

state = c("Hokkaido", 
          "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
          "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
          "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
          "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

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
             
savefig(file.path(dir.p,"Fig_3a"), width = 12, height = 10, toplines = 0.8, type = "png")
image(1:20, 1:47, ccf_value_female_prefecture[,rev(1:47)], ylab = "Prefecture", xlab = "Lag", 
      yaxt = "n", main = "Female")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig(file.path(dir.p,"Fig_3b"), width = 12, height = 10, toplines = 0.8, type = "png")
image(1:20, 1:47, ccf_value_male_prefecture[,rev(1:47)], ylab = "", xlab = "Lag", yaxt = "n", main = "Male")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig(file.path(dir.p,"Fig_3c"), width = 12, height = 10, toplines = 0.8, type = "png")
#image(1:20, 1:47, ccf_value_total_prefecture[,rev(1:47)], ylab = "", xlab = "Lag", yaxt = "n", main = "Total")
#axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))

filled.contour(1:20, 1:47, ccf_value_total_prefecture[,rev(1:47)], ylab = "", xlab = "Lag", 
               yaxt = "n", main = "Total", plot.axes = {axis(1)})
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()
