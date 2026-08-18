###############################################################################
# Script: read_data.R
# Description: Processes Japanese subnational mortality data and saves dx matrices
###############################################################################
# Load packages 
source("load_packages.R")

# Dynamic Directories using project-relative paths
dir_data     <- here("data")
dir_data_raw <- here("data", "updated_data")
dir_results  <- here("results")


dirs_to_create <- c(dir_data , dir_data_raw, dir_results) 

for (d in dirs_to_create) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# Define study constants
state <- c(
  "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
  "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa",
  "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", 
  "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", 
  "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
  "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", 
  "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa"
)

ages   <- 0:110
n_age  <- length(ages)
years  <- 1975:2023
n_year <- length(years)

# Load raw death probability (qx) data
female_prefecture_qx <- list()
male_prefecture_qx   <- list()
total_prefecture_qx  <- list()

read_qx_matrix <- function(gender_prefix, idx) {
  file_path <- file.path(dir_data_raw, paste0(gender_prefix, "_prefecture_", idx, ".txt"))
  df <- read.table(file_path, header = TRUE, skip = 2)
  df <- df[df$Year >= 1975, ]
  mat <- t(matrix(df$qx, nrow = n_age, ncol = n_year))
  colnames(mat) <- ages
  rownames(mat) <- years
  return(mat)
}

for (ik in seq_along(state)) {
  female_prefecture_qx[[ik]] <- read_qx_matrix("female", ik)
  male_prefecture_qx[[ik]]   <- read_qx_matrix("male", ik)
  total_prefecture_qx[[ik]]  <- read_qx_matrix("total", ik)
}

# Convert death probabilities (qx) to death counts (dx)
qx_to_dx <- function(qx_matrix, radix = 10^5) {
  dx_matrix <- matrix(NA_real_, nrow = nrow(qx_matrix), ncol = ncol(qx_matrix))
  for (i in seq_len(nrow(qx_matrix))) {
    qx <- qx_matrix[i, ]
    lx <- radix * c(1, cumprod(1 - qx[-length(qx)]))
    dx_matrix[i, ] <- lx * qx
  }
  colnames(dx_matrix) <- ages
  rownames(dx_matrix) <- years
  return(dx_matrix)
}

female_prefecture_dx <- lapply(female_prefecture_qx, qx_to_dx)
male_prefecture_dx   <- lapply(male_prefecture_qx, qx_to_dx)
total_prefecture_dx  <- lapply(total_prefecture_qx, qx_to_dx)

saveRDS(female_prefecture_dx, file = file.path(dir_data, "female_prefecture_dx.rds"))
saveRDS(male_prefecture_dx,   file = file.path(dir_data, "male_prefecture_dx.rds"))
saveRDS(total_prefecture_dx,  file = file.path(dir_data, "total_prefecture_dx.rds"))

# Compute national benchmark population distributions
compute_matrix_mean <- function(mat_list) {
  Reduce("+", mat_list) / length(mat_list)
}

Japan_female_pop <- compute_matrix_mean(female_prefecture_dx)
Japan_male_pop   <- compute_matrix_mean(male_prefecture_dx)
Japan_total_pop  <- compute_matrix_mean(total_prefecture_dx)

saveRDS(Japan_female_pop, file = file.path(dir_data, "Japan_female_pop.rds"))
saveRDS(Japan_male_pop,   file = file.path(dir_data, "Japan_male_pop.rds"))
saveRDS(Japan_total_pop,  file = file.path(dir_data, "Japan_total_pop.rds"))
