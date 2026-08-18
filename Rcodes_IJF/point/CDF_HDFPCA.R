
# Load required setup and functions
source("load_packages.R")
source(here("auxiliary_source", "auxiliary_point.R"))

# Base paths using project-relative paths
dir_data    <- here("data")
dir_results <- here("results", "Results_Point")
dir_shiny   <- here("results", "Shiny_App")
dirs_to_create <- c(dir_data, dir_results, dir_shiny)
for (d in dirs_to_create) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# Load processed data
female_prefecture_dx <- readRDS(file.path(dir_data, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(file.path(dir_data, "male_prefecture_dx.rds"))

# Parameters
state <- c(
  "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
  "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa",
  "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", 
  "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", 
  "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
  "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", 
  "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa"
)

n_state   <- length(state)
horizons  <- 17
n_metrics <- 4

###########################################################
# Primary Model Estimation (K = 6, L = 2)
###########################################################
point_fore_subnational_err_HDFPCA_F_CDF = point_fore_subnational_err_HDFPCA_M_CDF = array(NA, dim = c(47, 17, 4))
point_fore_subnational_HDFPCA_F_CDF = point_fore_subnational_HDFPCA_M_CDF = array(NA, dim = c(47, 111, 17))

for(ij in 1:47)
{
  
  for(iw in 1:17)
  {
    dum = hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], horizon = iw, 
                     first_order = 6, second_order = 2, transformation = "CDF")
    print(ij)
    print(iw)
    point_fore_subnational_err_HDFPCA_F_CDF[ij,iw,] = dum$err_F
    point_fore_subnational_err_HDFPCA_M_CDF[ij,iw,] = dum$err_M
    point_fore_subnational_HDFPCA_F_CDF[ij,,1:(18-iw)] = dum$forecast_pdf_F
    point_fore_subnational_HDFPCA_M_CDF[ij,,1:(18-iw)] = dum$forecast_pdf_M
    rm(dum); rm(iw)
  }
  print(ij); rm(ij)
}

# This partial results are saved for the Figure 8 of the paper. See Figure_8.R
saveRDS(point_fore_subnational_HDFPCA_F_CDF, file = file.path(dir_results, "HDFPCA_Forecast_F1.rds"))
saveRDS(point_fore_subnational_HDFPCA_M_CDF, file = file.path(dir_results, "HDFPCA_Forecast_M1.rds"))

###########################################################
# Sensitivity Analysis (K = 2, L = 2)
###########################################################
point_fore_subnational_err_HDFPCA_F_CDF_K2 <- array(NA, dim = c(n_state, horizons, n_metrics))
point_fore_subnational_err_HDFPCA_M_CDF_K2 <- array(NA, dim = c(n_state, horizons, n_metrics))

for (ij in 1:n_state) {
  for (iw in 1:horizons) {
    dum <- hdfpca_fun(
      fdata_F        = female_prefecture_dx[[ij]],
      fdata_M        = male_prefecture_dx[[ij]],
      horizon        = iw,
      first_order    = 2,
      second_order   = 2,
      transformation = "CDF"
    )
    
    point_fore_subnational_err_HDFPCA_F_CDF_K2[ij, iw, ] <- dum$err_F
    point_fore_subnational_err_HDFPCA_M_CDF_K2[ij, iw, ] <- dum$err_M
  }
}

###########################################################
# Summaries & Metric Computation
###########################################################
metric_names <- c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

# Summary across horizons (State level summary)
compute_summary_by_horizon <- function(err_array) {
  res <- rbind(
    apply(err_array, c(1, 3), mean),
    colMeans(apply(err_array, c(1, 3), mean)),
    apply(apply(err_array, c(1, 3), mean), 2, median)
  )
  rownames(res) <- c(state, "Mean", "Median")
  colnames(res) <- metric_names
  return(res)
}

# Summary across prefectures (Horizon level summary)
compute_summary_by_prefecture <- function(err_array) {
  res <- rbind(
    apply(err_array, c(2, 3), mean),
    colMeans(apply(err_array, c(2, 3), mean)),
    apply(apply(err_array, c(2, 3), mean), 2, median)
  )
  rownames(res) <- c(1:horizons, "Mean", "Median")
  colnames(res) <- metric_names
  return(res)
}

# Generate summary matrices
point_fore_subnational_err_HDFPCA_F_CDF_mean    <- compute_summary_by_horizon(point_fore_subnational_err_HDFPCA_F_CDF)
point_fore_subnational_err_HDFPCA_M_CDF_mean    <- compute_summary_by_horizon(point_fore_subnational_err_HDFPCA_M_CDF)
point_fore_subnational_err_HDFPCA_F_CDF_K2_mean <- compute_summary_by_horizon(point_fore_subnational_err_HDFPCA_F_CDF_K2)
point_fore_subnational_err_HDFPCA_M_CDF_K2_mean <- compute_summary_by_horizon(point_fore_subnational_err_HDFPCA_M_CDF_K2)

horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean <- compute_summary_by_prefecture(point_fore_subnational_err_HDFPCA_F_CDF)
horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean <- compute_summary_by_prefecture(point_fore_subnational_err_HDFPCA_M_CDF)

###########################################################
# Save Main Execution Results
###########################################################
saveRDS(point_fore_subnational_err_HDFPCA_F_CDF, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_F_CDF.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_M_CDF.rds"))

saveRDS(point_fore_subnational_err_HDFPCA_F_CDF_K2, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_F_CDF_K2.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF_K2, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_M_CDF_K2.rds"))

saveRDS(point_fore_subnational_err_HDFPCA_F_CDF_mean, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_F_CDF_mean.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF_mean, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_M_CDF_mean.rds"))

saveRDS(point_fore_subnational_err_HDFPCA_F_CDF_K2_mean, file = file.path(dir_results, "point_fore_subnational_err_HDFPCA_F_CDF_K2_mean.rds"))
saveRDS(point_fore_subnational_err_HDFPCA_M_CDF_K2_mean, file = file.path(dir_results, "point_fore_subnational_err_M_CDF_K2_mean.rds"))

saveRDS(horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean, file = file.path(dir_results, "horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean, file = file.path(dir_results, "horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean.rds"))

###########################################################
# Export Datasets for Shiny App
###########################################################
arrays_list <- list(
  F_EVR_ETS = point_fore_subnational_err_HDFPCA_F_CDF,
  M_EVR_ETS = point_fore_subnational_err_HDFPCA_M_CDF,
  F_K6_ETS  = point_fore_subnational_err_HDFPCA_F_CDF_K2,
  M_K6_ETS  = point_fore_subnational_err_HDFPCA_M_CDF_K2
)

for (name in names(arrays_list)) {
  arr    <- arrays_list[[name]]
  gender <- ifelse(grepl("^F", name), "female", "male")
  comp   <- ifelse(grepl("EVR", name), "EVR", "K")
  
  for (i in 1:2) {
    err_metric <- ifelse(i == 1, "KLD", "JSD")
    mat        <- t(arr[, , i])
    file_name  <- paste0("HDFPCA_", err_metric, "_", gender, "_", comp, ".rds")
    
    saveRDS(mat, file = file.path(dir_shiny, file_name))
  }
}

