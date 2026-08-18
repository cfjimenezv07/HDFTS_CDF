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
ages      <- 0:110

# Initialize output containers
point_fore_subnational_err_MLFTS_F_EVR_ETS <- array(NA, dim = c(n_state, horizons, n_metrics))
point_fore_subnational_err_MLFTS_M_EVR_ETS <- array(NA, dim = c(n_state, horizons, n_metrics))
point_fore_subnational_err_MLFTS_F_K6_ETS  <- array(NA, dim = c(n_state, horizons, n_metrics))
point_fore_subnational_err_MLFTS_M_K6_ETS  <- array(NA, dim = c(n_state, horizons, n_metrics))

# Compute point forecast errors
for (ij in 1:n_state) {
  for (iw in 1:horizons) {
    
    dum_evr <- point_fore_national_cdf_MLFTS(
      fdata_F = female_prefecture_dx[[ij]],
      fdata_M = male_prefecture_dx[[ij]],
      fdata_common = NULL,
      fore_method = "CDF", 
      horizon = iw, 
      way_ncomp = "EVR"
    )
    point_fore_subnational_err_MLFTS_F_EVR_ETS[ij, iw, ] <- dum_evr$err_F
    point_fore_subnational_err_MLFTS_M_EVR_ETS[ij, iw, ] <- dum_evr$err_M
    
    dum_k6 <- point_fore_national_cdf_MLFTS(
      fdata_F = female_prefecture_dx[[ij]],
      fdata_M = male_prefecture_dx[[ij]],
      fdata_common = NULL,
      fore_method = "CDF", 
      horizon = iw, 
      way_ncomp = "provide"
    )
    point_fore_subnational_err_MLFTS_F_K6_ETS[ij, iw, ] <- dum_k6$err_F
    point_fore_subnational_err_MLFTS_M_K6_ETS[ij, iw, ] <- dum_k6$err_M
  }
}

dim_list <- list(state, 1:horizons, c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2"))
dimnames(point_fore_subnational_err_MLFTS_F_EVR_ETS) <- dim_list
dimnames(point_fore_subnational_err_MLFTS_M_EVR_ETS) <- dim_list
dimnames(point_fore_subnational_err_MLFTS_F_K6_ETS)  <- dim_list
dimnames(point_fore_subnational_err_MLFTS_M_K6_ETS)  <- dim_list

# Compute summary statistics across horizons
compute_summary_by_horizon <- function(err_array) {
  round(rbind(
    apply(err_array, c(1, 3), mean),
    colMeans(apply(err_array, c(1, 3), mean)),
    apply(apply(err_array, c(1, 3), mean), 2, median)
  ), 4)
}

point_fore_subnational_err_MLFTS_F_EVR_ETS_mean <- compute_summary_by_horizon(point_fore_subnational_err_MLFTS_F_EVR_ETS)
point_fore_subnational_err_MLFTS_F_K6_ETS_mean  <- compute_summary_by_horizon(point_fore_subnational_err_MLFTS_F_K6_ETS)
point_fore_subnational_err_MLFTS_M_EVR_ETS_mean <- compute_summary_by_horizon(point_fore_subnational_err_MLFTS_M_EVR_ETS)
point_fore_subnational_err_MLFTS_M_K6_ETS_mean  <- compute_summary_by_horizon(point_fore_subnational_err_MLFTS_M_K6_ETS)

rownames(point_fore_subnational_err_MLFTS_F_EVR_ETS_mean) <- c(state, "Mean", "Median")
rownames(point_fore_subnational_err_MLFTS_F_K6_ETS_mean)  <- c(state, "Mean", "Median")
rownames(point_fore_subnational_err_MLFTS_M_EVR_ETS_mean) <- c(state, "Mean", "Median")
rownames(point_fore_subnational_err_MLFTS_M_K6_ETS_mean)  <- c(state, "Mean", "Median")

# Compute summary statistics across prefectures
compute_summary_by_prefecture <- function(err_array) {
  rbind(
    apply(err_array, c(2, 3), mean),
    colMeans(apply(err_array, c(2, 3), mean)),
    apply(apply(err_array, c(2, 3), mean), 2, median)
  )
}

horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean <- compute_summary_by_prefecture(point_fore_subnational_err_MLFTS_F_EVR_ETS)
horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean  <- compute_summary_by_prefecture(point_fore_subnational_err_MLFTS_F_K6_ETS)
horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean <- compute_summary_by_prefecture(point_fore_subnational_err_MLFTS_M_EVR_ETS)
horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean  <- compute_summary_by_prefecture(point_fore_subnational_err_MLFTS_M_K6_ETS)

rownames(horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean) <- c(1:horizons, "Mean", "Median")
rownames(horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean)  <- c(1:horizons, "Mean", "Median")
rownames(horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean) <- c(1:horizons, "Mean", "Median")
rownames(horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean)  <- c(1:horizons, "Mean", "Median")

# Save primary result files
saveRDS(point_fore_subnational_err_MLFTS_F_EVR_ETS, file = file.path(dir_results, "point_fore_subnational_err_MLFTS_F_EVR_ETS.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_EVR_ETS, file = file.path(dir_results, "point_fore_subnational_err_MLFTS_M_EVR_ETS.rds"))
saveRDS(point_fore_subnational_err_MLFTS_F_K6_ETS,  file = file.path(dir_results, "point_fore_subnational_err_MLFTS_F_K6_ETS.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_K6_ETS,  file = file.path(dir_results, "point_fore_subnational_err_MLFTS_M_K6_ETS.rds"))

saveRDS(point_fore_subnational_err_MLFTS_F_EVR_ETS_mean, file = file.path(dir_results, "point_fore_subnational_err_MLFTS_F_EVR_ETS_mean.rds"))
saveRDS(point_fore_subnational_err_MLFTS_F_K6_ETS_mean,  file = file.path(dir_results, "point_fore_subnational_err_MLFTS_F_K6_ETS_mean.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_EVR_ETS_mean, file = file.path(dir_results, "point_fore_subnational_err_MLFTS_M_EVR_ETS_mean.rds"))
saveRDS(point_fore_subnational_err_MLFTS_M_K6_ETS_mean,  file = file.path(dir_results, "point_fore_subnational_err_MLFTS_M_K6_ETS_mean.rds"))

saveRDS(horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean, file = file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean,  file = file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean, file = file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean.rds"))
saveRDS(horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean,  file = file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean.rds"))

# Prepare interactive visualizations exports
arrays_list <- list(
  F_EVR_ETS = point_fore_subnational_err_MLFTS_F_EVR_ETS,
  M_EVR_ETS = point_fore_subnational_err_MLFTS_M_EVR_ETS,
  F_K6_ETS  = point_fore_subnational_err_MLFTS_F_K6_ETS,
  M_K6_ETS  = point_fore_subnational_err_MLFTS_M_K6_ETS
)

for (name in names(arrays_list)) {
  arr    <- arrays_list[[name]]
  gender <- ifelse(grepl("^F", name), "female", "male")
  comp   <- ifelse(grepl("EVR", name), "EVR", "K")
  
  for (i in 1:2) {
    err_metric <- ifelse(i == 1, "KLD", "JSD")
    mat        <- t(arr[, , i])
    file_name  <- paste0("MLFTS_", err_metric, "_", gender, "_", comp, ".rds")
    
    saveRDS(mat, file = file.path(dir_shiny, file_name))
  }
}

