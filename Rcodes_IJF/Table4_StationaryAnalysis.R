###############################################################################
# Script: Table4_StationaryAnalysis.R
# Description: Summary of the stationary analysis of the FTS
#              Produces Table 4 
###############################################################################
source("load_packages.R")

# Dynamic Directories
dir.d       <- here("data", "updated_data")
dir.results <- here("results")
dir.p       <- here("results", "Plots")
dir.aux     <- here("auxiliary_source")

# Ensure output directories exist
if (!dir.exists(dir.results)) dir.create(dir.results, recursive = TRUE)
if (!dir.exists(dir.p))       dir.create(dir.p, recursive = TRUE)

# 2. PARAMETERS & DATASET INITIALIZATION
# ------------------------------------------------------------------------------
year         <- 1975:2023
n_year       <- length(year)
age          <- 0:110
n_age        <- length(age)
n_prefectures<- 47
n_populations<- 2

# Row Partitioning (for splitting concatenated arrays back to 47 prefectures)
part_list <- lapply(1:n_prefectures, function(ik) {
  (n_year * ik - (n_year - 1)):(n_year * ik)
})

state <- c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
           "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
           "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
           "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
           "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
           "Tokushima", "Kagawa", "Ehime", "Kochi",
           "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

# 3. READ SUBNATIONAL MORTALITY DATA & CALCULATE DX
# ------------------------------------------------------------------------------
female_prefecture_qx <- male_prefecture_qx <- total_prefecture_qx <- list()

for (ik in seq_along(state)) {
  # Female
  f_data <- read.table(file.path(dir.d, paste0("female_prefecture_", ik, ".txt")), header = TRUE, skip = 2)
  female_prefecture_qx[[ik]] <- t(matrix(f_data[f_data$Year >= 1975, ]$qx, n_age, n_year))
  
  # Male
  m_data <- read.table(file.path(dir.d, paste0("male_prefecture_", ik, ".txt")), header = TRUE, skip = 2)
  male_prefecture_qx[[ik]] <- t(matrix(m_data[m_data$Year >= 1975, ]$qx, n_age, n_year))
  
  # Total
  t_data <- read.table(file.path(dir.d, paste0("total_prefecture_", ik, ".txt")), header = TRUE, skip = 2)
  total_prefecture_qx[[ik]] <- t(matrix(t_data[t_data$Year >= 1975, ]$qx, n_age, n_year))
}

# Construct Dx matrices
female_prefecture_dx <- male_prefecture_dx <- total_prefecture_dx <- list()

for (iw in seq_along(state)) {
  f_dum <- m_dum <- t_dum <- matrix(NA, n_year, n_age)
  for (ij in 1:n_year) {
    pop_f <- pop_m <- pop_t <- 10^5
    for (ik in 1:n_age) {
      f_dum[ij, ik] <- female_prefecture_qx[[iw]][ij, ik] * pop_f
      pop_f <- pop_f - f_dum[ij, ik]
      
      m_dum[ij, ik] <- male_prefecture_qx[[iw]][ij, ik] * pop_m
      pop_m <- pop_m - m_dum[ij, ik]
      
      t_dum[ij, ik] <- total_prefecture_qx[[iw]][ij, ik] * pop_t
      pop_t <- pop_t - t_dum[ij, ik]
    }
  }
  female_prefecture_dx[[iw]] <- t(f_dum)
  male_prefecture_dx[[iw]]   <- t(m_dum)
  total_prefecture_dx[[iw]]  <- t(t_dum)
}

All_Japan_female_qx <- (female_prefecture_dx)
All_Japan_male_qx   <- (male_prefecture_dx)

# 4. CDF TRANSFORMATION
# ------------------------------------------------------------------------------
source(file.path(dir.aux, "CDF_transformation.R"))

Transformed_female <- list()
Transformed_male   <- list()

for (i in 1:n_prefectures) {
  Transformed_female[[i]] <- t(cdf_transformation(t(All_Japan_female_qx[[i]]) / 10^5, year))
  Transformed_male[[i]]   <- t(cdf_transformation(t(All_Japan_male_qx[[i]]) / 10^5, year))
}

all_unconstrained_female <- t(list.cbind(Transformed_female))
all_unconstrained_male   <- t(list.cbind(Transformed_male))

# 5. FUNCTIONAL MEAN ANOVA (FM-ANOVA) DECOMPOSITION
# ------------------------------------------------------------------------------
new_age <- 0:109

# Compute decomposition and functional residuals
Residuals_means <- hdftsa::Two_way_mean_residuals(
  data_pop1 = t(all_unconstrained_male),
  data_pop2 = t(all_unconstrained_female),
  year = year,
  age = new_age,
  n_prefectures = n_prefectures,
  n_populations = n_populations
)

Res1_means <- Residuals_means$residuals1_mean
Res2_means <- Residuals_means$residuals2_mean

# Split residuals per prefecture and transpose (110 grid points x 49 time points)
male_res   <- lapply(part_list, function(k) t(Res1_means[k, ]))
female_res <- lapply(part_list, function(k) t(Res2_means[k, ]))



# --- 6. Functional Stationarity Tests (ftsa) ---
center_matrix_corrected <- function(mat) {
  col_means <- colMeans(mat)
  sweep(mat, MARGIN = 2, STATS = col_means, FUN = "-")
}

run_prefecture_stationarity <- function(data_list, alpha = 0.05) {
  bandwidth_h <- 110^0.5
  pvals <- numeric(length(data_list))
  
  for (i in seq_along(data_list)) {
    mat_centered <- center_matrix_corrected(data_list[[i]])
    output_text <- capture.output(
      T_stationary(
        sample = mat_centered, L = 49, J = 500, MC_rep = 1000,
        cumulative_var = 0.90, Ker1 = FALSE, Ker2 = TRUE,
        h = bandwidth_h, pivotal = FALSE, use_table = FALSE, significance = alpha
      )
    )
    pval_line <- grep("p-value", output_text, value = TRUE)
    pvals[i]  <- as.numeric(sub(".*p-value = ([0-9\\.eE-]+).*", "\\1", pval_line))
  }
  return(pvals)
}

male_stat_pvals   <- run_prefecture_stationarity(male_res)
female_stat_pvals <- run_prefecture_stationarity(female_res)

stat_summary_table <- data.frame(
  Group   = c("Male", "Female"),
  Min     = c(min(male_stat_pvals), min(female_stat_pvals)),
  `1st Qu.` = c(quantile(male_stat_pvals, 0.25), quantile(female_stat_pvals, 0.25)),
  Median  = c(median(male_stat_pvals), median(female_stat_pvals)),
  Mean    = c(mean(male_stat_pvals), mean(female_stat_pvals)),
  `3rd Qu.` = c(quantile(male_stat_pvals, 0.75), quantile(female_stat_pvals, 0.75)),
  Max     = c(max(male_stat_pvals), max(female_stat_pvals))
)

print(stat_summary_table)
# write.csv(stat_summary_table, file.path(dir.results, "stationarity_summary.csv"), row.names = FALSE)


# --- C. Heteroscedasticity Tests & Plots (FTSgof) ---
get_fCH_pval <- function(mat, H = 10) {
  out <- fCH_test(mat, H = H, stat_Method = "functional", pplot = FALSE)
  return(out$p_value)
}

female_fCH_pvals <- sapply(female_res, get_fCH_pval)
male_fCH_pvals   <- sapply(male_res,   get_fCH_pval)

df_fCH <- data.frame(
  p_value = c(female_fCH_pvals, male_fCH_pvals),
  sex     = rep(c("Female", "Male"), each = n_prefectures)
)

