###########################################################
# High-Dimensional Functional Principal Component Analysis
# transformation = "CDF"
###########################################################
source("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/auxiliary_source/auxiliary_point.R")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/HLS/"
dir.r1 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/FANOVA/"

female_prefecture_dx <- readRDS(paste0(dir.r1, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(paste0(dir.r1, "male_prefecture_dx.rds"))
total_prefecture_dx  <- readRDS(paste0(dir.r1, "total_prefecture_dx.rds"))

# K = 6, L = 2

point_fore_subnational_err_HDFPCA_F_CDF = point_fore_subnational_err_HDFPCA_M_CDF = array(NA, dim = c(47, 17, 4))
point_fore_subnational_HDFPCA_F_CDF = point_fore_subnational_HDFPCA_M_CDF = array(NA, dim = c(47, 111, 17))
for(ij in 1:47)
{
  
    for(iw in 1:17)
    {
        dum = hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], horizon = iw, 
                         first_order = 6, second_order = 2, transformation = "CDF")
        point_fore_subnational_err_HDFPCA_F_CDF[ij,iw,] = dum$err_F
        point_fore_subnational_err_HDFPCA_M_CDF[ij,iw,] = dum$err_M
        point_fore_subnational_HDFPCA_F_CDF[ij,,] = dum$forecast_pdf_F
        point_fore_subnational_HDFPCA_M_CDF[ij,,] = dum$forecast_pdf_M
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

# saveRDS(point_fore_subnational_HDFPCA_F_CDF, file = paste0(dir.r1, "HDFPCA_Forecast_F1.rds"))
# saveRDS(point_fore_subnational_HDFPCA_M_CDF,   file = paste0(dir.r1, "HDFPCA_Forecast_M1.rds"))


###########################################
# For sensitivity analaysis, K = 2, L = 2
###########################################

point_fore_subnational_err_HDFPCA_F_CDF_K2 = point_fore_subnational_err_HDFPCA_M_CDF_K2 = array(NA, dim = c(47, 17, 4))
for(ij in 1:47)
{
    for(iw in 1:17)
    {
        dum = hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], horizon = iw, 
                         first_order = 2, second_order = 2, transformation = "CDF")
        point_fore_subnational_err_HDFPCA_F_CDF_K2[ij,iw,] = dum$err_F
        point_fore_subnational_err_HDFPCA_M_CDF_K2[ij,iw,] = dum$err_M
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

point_fore_subnational_err_HDFPCA_F_CDF_K2_mean = rbind(apply(point_fore_subnational_err_HDFPCA_F_CDF_K2, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_F_CDF_K2, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_F_CDF_K2, c(1,3), mean), 2, median)) # 0.0157 0.0047
rownames(point_fore_subnational_err_HDFPCA_F_CDF_K2_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_F_CDF_K2_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

point_fore_subnational_err_HDFPCA_M_CDF_K2_mean = rbind(apply(point_fore_subnational_err_HDFPCA_M_CDF_K2, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_M_CDF_K2, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_M_CDF_K2, c(1,3), mean), 2, median)) # 0.0141 0.0041
rownames(point_fore_subnational_err_HDFPCA_M_CDF_K2_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_M_CDF_K2_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

###########
## summary
###########

# F

point_fore_subnational_err_HDFPCA_F_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(1,3), mean), 2, median)) # 0.0157 0.0047
rownames(point_fore_subnational_err_HDFPCA_F_CDF_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_F_CDF_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(2, 3), mean), 
                                              colMeans(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(2, 3), mean)),
                                              apply(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(2, 3), mean), 2, median))
rownames(horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean) = c(1:17, "Mean", "Median")

# M

point_fore_subnational_err_HDFPCA_M_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(1,3), mean), 2, median)) # 0.0141 0.0041
rownames(point_fore_subnational_err_HDFPCA_M_CDF_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_M_CDF_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")


horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(2, 3), mean),
                                            colMeans(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(2, 3), mean)),
                                            apply(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean) = c(1:17, "Mean", "Median")


saveRDS(point_fore_subnational_err_HDFPCA_F_CDF, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF.rds"))

saveRDS(point_fore_subnational_err_HDFPCA_F_CDF_K2, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF_K2.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF_K2, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF_K2.rds"))

saveRDS(point_fore_subnational_err_HDFPCA_F_CDF_mean, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF_mean.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF_mean, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF_mean.rds"))

saveRDS(point_fore_subnational_err_HDFPCA_F_CDF_K2_mean, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF_K2_mean.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF_K2_mean, file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF_K2_mean.rds"))

saveRDS(horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean, file = paste0(dir.r, "horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean, file = paste0(dir.r, "horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean.rds"))


#Read them
dir.shiny <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Shiny_App_CDF/datasets_shiny_app2/"


point_fore_subnational_err_HDFPCA_F_CDF <- readRDS(paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF.rds"))
point_fore_subnational_err_HDFPCA_M_CDF <- readRDS(paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF.rds"))

point_fore_subnational_err_HDFPCA_F_CDF_K2 <- readRDS(paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF_K2.rds"))
point_fore_subnational_err_HDFPCA_M_CDF_K2 <- readRDS(paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF_K2.rds"))




# Directory to save
dir.shiny <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Shiny_App_CDF/datasets_shiny_app2/"

# List of arrays to process
arrays_list <- list(
  F_EVR_ETS = point_fore_subnational_err_HDFPCA_F_CDF,
  M_EVR_ETS = point_fore_subnational_err_HDFPCA_M_CDF,
  F_K6_ETS  = point_fore_subnational_err_HDFPCA_F_CDF_K2,
  M_K6_ETS  = point_fore_subnational_err_HDFPCA_M_CDF_K2
)

# Loop through arrays
for(name in names(arrays_list)) {
  
  arr <- arrays_list[[name]]
  
  # Gender and component from name
  gender <- ifelse(grepl("^F", name), "female", "male")
  comp <- ifelse(grepl("EVR", name), "EVR", "K")
  
  # Loop through third dimension (1 = KLD, 2 = JSD)
  for(i in 1:dim(arr)[3]) {
    
    # Error metric
    err_metric <- ifelse(i == 1, "KLD", "JSD")
    
    # Extract matrix and transpose to 17x47
    mat <- t(arr[,,i])
    
    # Construct file name
    file_name <- paste0("HDFPCA_", err_metric, "_", gender, "_", comp, ".rds")
    
    # Save as .rds in dir.shiny
    saveRDS(mat, file = file.path(dir.shiny, file_name))
  }
}


