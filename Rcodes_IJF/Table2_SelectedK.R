###############################################################################
# Script: table2_component_selection.R
# Description: EVR Component Selection for UFTS, MFTS, MLFTS, and FANOVA
#              Produces Table 2 (Selected K counts across forecasting horizons)
###############################################################################
source("load_packages.R")

# Relative Directory Paths
dir_data    <- here("data")
dir_source  <- here("auxiliary_source")
dir_results <- here("results")

# Load project auxiliary functions
source(file.path(dir_source, "auxiliary_point.R"))
source(file.path(dir_source, "CDF_transformation.R"))

# Ensure results export subdirectories exist
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

# 2. Data Loading & Constants ------------------------------------------------
female_prefecture_dx <- readRDS(file.path(dir_data, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(file.path(dir_data, "male_prefecture_dx.rds"))
total_prefecture_dx  <- readRDS(file.path(dir_data, "total_prefecture_dx.rds"))

ages   <- 0:110
year  <- 1975:2023
n_age  <- length(ages)
n_year <- length(year)
n_prefectures <- length(female_prefecture_dx) # 47
n_populations <- 2                             # Female & Male

# 3. Preprocessing Functions & Helpers ---------------------------------------

# Eigenvalue Ratio Criterion (EVR) component selector
select_K <- function(tau, eigenvalue) {
  k_max <- length(eigenvalue)
  k_all <- numeric(k_max - 1)
  for (k in 1:(k_max - 1)) {
    k_all[k] <- (eigenvalue[k + 1] / eigenvalue[k]) * ifelse(eigenvalue[k] / eigenvalue[1] > tau, 1, 0) + 
      ifelse(eigenvalue[k] / eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

# Cumulative Density Function preparation on logit scale
prep_data_safe <- function(d) {
  d_norm   <- d / 10^5
  d_cumsum <- t(apply(d_norm, 1, cumsum))
  eps      <- 1e-10
  d_cumsum[d_cumsum <= 0] <- eps
  d_cumsum[d_cumsum >= 1] <- 1 - eps
  return(t(apply(d_cumsum[, 1:(ncol(d) - 1)], 1, logit)))
}

# 4. Standard Methods Component Selection Loop (UFTS, MFTS, MLFTS) -----------

# Initialize selection matrices (47 Prefectures x 17 Horizons)
ncomp_UFTS_F     <- matrix(NA, n_prefectures, 17) # Univariate Female
ncomp_UFTS_M     <- matrix(NA, n_prefectures, 17) # Univariate Male
ncomp_MFTS       <- matrix(NA, n_prefectures, 17) # Multivariate (Joint)
ncomp_MLFTS_Agg  <- matrix(NA, n_prefectures, 17) # Multilevel: Common Trend
ncomp_MLFTS_ResF <- matrix(NA, n_prefectures, 17) # Multilevel: Female Residuals
ncomp_MLFTS_ResM <- matrix(NA, n_prefectures, 17) # Multilevel: Male Residuals

for (ij in 1:n_prefectures) {
  data_F_full <- prep_data_safe(female_prefecture_dx[[ij]])
  data_M_full <- prep_data_safe(male_prefecture_dx[[ij]])
  
  for (iw in 1:17) {
    # Expanding window starting at 33 observations
    end_row <- min(32 + iw, nrow(data_F_full))
    W_F     <- data_F_full[1:end_row, ]
    W_M     <- data_M_full[1:end_row, ]
    
    # --- UFTS (Separate) ---
    ncomp_UFTS_F[ij, iw] <- select_K(tau = 10^-3, eigenvalue = (svd(W_F)$d)^2)
    ncomp_UFTS_M[ij, iw] <- select_K(tau = 10^-3, eigenvalue = (svd(W_M)$d)^2)
    
    # --- MFTS (Jointly Stacked) ---
    comb_MFTS          <- rbind(t(W_F), t(W_M))
    ncomp_MFTS[ij, iw] <- select_K(tau = 10^-2, eigenvalue = eigen(cov(t(comb_MFTS)))$values)
    
    # --- MLFTS (Common + Specific Residuals) ---
    S_F      <- t(scale(t(W_F), center = TRUE, scale = FALSE))
    S_M      <- t(scale(t(W_M), center = TRUE, scale = FALSE))
    agg_data <- (S_F + S_M) / 2
    
    K_agg                   <- select_K(tau = 10^-3, eigenvalue = eigen(cov(t(agg_data)))$values)
    ncomp_MLFTS_Agg[ij, iw] <- K_agg
    
    ftsm_agg                 <- ftsm(fts(1:nrow(agg_data), agg_data), order = K_agg)
    ncomp_MLFTS_ResF[ij, iw] <- select_K(tau = 10^-3, eigenvalue = eigen(cov(t(S_F - ftsm_agg$fitted$y)))$values)
    ncomp_MLFTS_ResM[ij, iw] <- select_K(tau = 10^-3, eigenvalue = eigen(cov(t(S_M - ftsm_agg$fitted$y)))$values)
  }
}

# Export standard model matrices
# saveRDS(ncomp_UFTS_F,    file = file.path(dir_results, "ncomp_UFTS_F.rds"))
# saveRDS(ncomp_UFTS_M,    file = file.path(dir_results, "ncomp_UFTS_M.rds"))
# saveRDS(ncomp_MFTS,      file = file.path(dir_results, "ncomp_MFTS.rds"))
# saveRDS(ncomp_MLFTS_Agg, file = file.path(dir_results, "ncomp_MLFTS_Agg.rds"))

# 5. Functional ANOVA Component Selection (FANOVA) ---------------------------

# Compute CDF Transformations across subnational series
Transformed_female <- list()
Transformed_male   <- list()

for (i in 1:n_prefectures) {
  Transformed_female[[i]] <- t(cdf_transformation((female_prefecture_dx[[i]]) / 10^5, year))
  Transformed_male[[i]]   <- t(cdf_transformation((male_prefecture_dx[[i]]) / 10^5, year))
}

all_unconstrained_female <- t(list.cbind(Transformed_female))
all_unconstrained_male   <- t(list.cbind(Transformed_male))

# Two-way Functional ANOVA Residuals
new_age         <- 0:109
Residuals_means <- hdftsa::Two_way_mean_residuals(
  data_pop1     = t(all_unconstrained_male), 
  data_pop2     = t(all_unconstrained_female), 
  year          = year, 
  age           = new_age, 
  n_prefectures      = n_prefectures, 
  n_populations = n_populations
)

Residuals_mean <- cbind(Residuals_means$residuals1_mean, Residuals_means$residuals2_mean)

# FANOVA Selection Loop
ncomp_ANOVA           <- matrix(NA, n_prefectures, 17)
rownames(ncomp_ANOVA) <- paste0("Pref_", 1:n_prefectures)
colnames(ncomp_ANOVA) <- paste0("h_", 1:17)

for (ij in 1:n_prefectures) {
  pref_rows_anova    <- (n_year * ij - (n_year - 1)):(n_year * ij)
  Residuals_f1_anova <- Residuals_mean[pref_rows_anova, ]
  
  for (iw in 1:17) {
    end_row           <- min(32 + iw, nrow(Residuals_f1_anova))
    current_res_anova <- Residuals_f1_anova[1:end_row, ]
    evs_anova         <- eigen(cov(current_res_anova))$values
    
    ncomp_ANOVA[ij, iw] <- select_K(tau = 10^-2, eigenvalue = evs_anova)
  }
}

saveRDS(ncomp_ANOVA, file = file.path(dir_results, "ncomp_ANOVA.rds"))

# 6. Tabulation & Summary Output (Table 2 Formatting) ------------------------

all_mats <- list(
  UFTS_F = ncomp_UFTS_F,
  UFTS_M = ncomp_UFTS_M,
  MFTS   = ncomp_MFTS,
  MLFTS  = ncomp_MLFTS_Agg,
  FANOVA = ncomp_ANOVA
)

# Tabulate K values across forecasting horizons h = 1..17
horizon_results <- lapply(1:17, function(h) {
  method_counts <- sapply(all_mats, function(mat) {
    tabulate(as.numeric(mat[, h]), nbins = 2)
  })
  return(as.vector(method_counts))
})

final_tab <- do.call(rbind, horizon_results)
colnames(final_tab) <- paste0(rep(names(all_mats), each = 2), c("_K1", "_K2"))
rownames(final_tab) <- paste0("h=", 1:17)

# Output summary table
print(final_tab)
