####################################
# Multilevel functional time series
####################################
source("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/auxiliary_source/auxiliary_point.R")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/HLS/"



point_fore_subnational_err_MLFTS_F_EVR_ETS = point_fore_subnational_err_MLFTS_M_EVR_ETS =
point_fore_subnational_err_MLFTS_F_K6_ETS  = point_fore_subnational_err_MLFTS_M_K6_ETS  = array(NA, dim = c(47, 17, 4))
for(ij in 1:47)
{
    for(iw in 1:17)
    {
        ## EVR
        
        dum = point_fore_national_cdf_MLFTS(fdata_F = female_prefecture_dx[[ij]],
                                            fdata_M = male_prefecture_dx[[ij]],
                                            fdata_common = NULL, fore_method = "CDF",
                                            horizon = iw, way_ncomp = "EVR")
        point_fore_subnational_err_MLFTS_F_EVR_ETS[ij,iw,] = dum$err_F
        point_fore_subnational_err_MLFTS_M_EVR_ETS[ij,iw,] = dum$err_M
        rm(dum)
        
        ## K = 6
        
        dum = point_fore_national_cdf_MLFTS(fdata_F = female_prefecture_dx[[ij]],
                                           fdata_M = male_prefecture_dx[[ij]],
                                           fdata_common = NULL, fore_method = "CDF",
                                           horizon = iw, way_ncomp = "provide")
        point_fore_subnational_err_MLFTS_F_K6_ETS[ij,iw,] = dum$err_F
        point_fore_subnational_err_MLFTS_M_K6_ETS[ij,iw,] = dum$err_M
        rm(dum)
    }
    print(ij); rm(ij)
}
dimnames(point_fore_subnational_err_MLFTS_F_EVR_ETS) = dimnames(point_fore_subnational_err_MLFTS_M_EVR_ETS) =
dimnames(point_fore_subnational_err_MLFTS_F_K6_ETS) = dimnames(point_fore_subnational_err_MLFTS_M_K6_ETS) = list(state, 1:17, c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2"))

#######
## ETS
#######

## female

# average across horizons

point_fore_subnational_err_MLFTS_F_EVR_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(1, 3), mean),
                                                             colMeans(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(1, 3), mean)),
                                                             apply(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(1, 3), mean), 2, median)), 4)

point_fore_subnational_err_MLFTS_F_K6_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(1, 3), mean),
                                                            colMeans(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(1, 3), mean)),
                                                            apply(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(1, 3), mean), 2, median)), 4)
rownames(point_fore_subnational_err_MLFTS_F_EVR_ETS_mean) = rownames(point_fore_subnational_err_MLFTS_F_K6_ETS_mean) = c(state, "Mean", "Median") 

# average across prefectures

horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(2, 3), mean), 2, median))

horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean) = rownames(horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean) = c(1:17, "Mean", "Median")

## male

# average across horizons

point_fore_subnational_err_MLFTS_M_EVR_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(1, 3), mean),
                                                             colMeans(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(1, 3), mean)),
                                                             apply(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(1, 3), mean), 2, median)), 4)

point_fore_subnational_err_MLFTS_M_K6_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(1, 3), mean),
                                                            colMeans(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(1, 3), mean)),
                                                            apply(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(1, 3), mean), 2, median)), 4)
rownames(point_fore_subnational_err_MLFTS_M_EVR_ETS_mean) = rownames(point_fore_subnational_err_MLFTS_M_K6_ETS_mean) = c(state, "Mean", "Median")

# average across prefectures

horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(2, 3), mean), 2, median))

horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean) = rownames(horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean) = c(1:17, "Mean", "Median")



saveRDS(point_fore_subnational_err_MLFTS_F_EVR_ETS, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_EVR_ETS.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_EVR_ETS, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_EVR_ETS.rds"))
saveRDS(point_fore_subnational_err_MLFTS_F_K6_ETS, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_K6_ETS.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_K6_ETS, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_K6_ETS.rds"))

saveRDS(point_fore_subnational_err_MLFTS_F_EVR_ETS_mean, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_EVR_ETS_mean.rds"))
saveRDS(point_fore_subnational_err_MLFTS_F_K6_ETS_mean, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_K6_ETS_mean.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_EVR_ETS_mean, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_EVR_ETS_mean.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_K6_ETS_mean, file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_K6_ETS_mean.rds"))

saveRDS(horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean, file = paste0(dir.r, "horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean, file = paste0(dir.r, "horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean, file = paste0(dir.r, "horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean, file = paste0(dir.r, "horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean.rds"))


#Read them
dir.shiny <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Shiny_App_CDF/datasets_shiny_app2/"


point_fore_subnational_err_MLFTS_F_EVR_ETS <- readRDS(paste0(dir.r, "point_fore_subnational_err_MLFTS_F_EVR_ETS.rds"))
point_fore_subnational_err_MLFTS_M_EVR_ETS <- readRDS(paste0(dir.r, "point_fore_subnational_err_MLFTS_M_EVR_ETS.rds"))
point_fore_subnational_err_MLFTS_F_K6_ETS <- readRDS(paste0(dir.r, "point_fore_subnational_err_MLFTS_F_K6_ETS.rds"))
point_fore_subnational_err_MLFTS_M_K6_ETS <- readRDS(paste0(dir.r, "point_fore_subnational_err_MLFTS_M_K6_ETS.rds"))



# Directory to save
dir.shiny <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Shiny_App_CDF/datasets_shiny_app2/"

# List of arrays to process
arrays_list <- list(
  F_EVR_ETS = point_fore_subnational_err_MLFTS_F_EVR_ETS,
  M_EVR_ETS = point_fore_subnational_err_MLFTS_M_EVR_ETS,
  F_K6_ETS  = point_fore_subnational_err_MLFTS_F_K6_ETS,
  M_K6_ETS  = point_fore_subnational_err_MLFTS_M_K6_ETS
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
    file_name <- paste0("MLFTS_", err_metric, "_", gender, "_", comp, ".rds")
    
    # Save as .rds in dir.shiny
    saveRDS(mat, file = file.path(dir.shiny, file_name))
  }
}

