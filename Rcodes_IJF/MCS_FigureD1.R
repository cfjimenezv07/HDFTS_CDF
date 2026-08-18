# ==============================================================================
# MODEL CONFIDENCE SET (MCS) ANALYSIS
# ==============================================================================
source("load_packages.R")

# Subdirectories using project-relative paths
dir.r <- here("results", "Results_Point")
dir.p <- here("results", "Plots")

# Ensure output directory exists
if (!dir.exists(dir.p)) dir.create(dir.p, recursive = TRUE)

# Model & Horizon Parameters
n_methods  <- 9
n_horizons <- 17
n_pref     <- 47

method_names <- c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA",
                  "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)", "FANOVA(K=6)")

alpha     <- 0.10
B         <- 2000
statistic <- "Tmax"

# 2. HELPER FUNCTIONS & ARRAY INITIALIZATION
# ------------------------------------------------------------------------------
# Helper function to safely read and transpose RDS matrix objects
load_t <- function(file_name, idx = 1) {
  full_path <- file.path(dir.r, file_name)
  if (!file.exists(full_path)) {
    stop(paste("File missing:", full_path))
  }
  t(readRDS(full_path)[,, idx])
}

load_t1 <- function(file_name) {
  full_path <- file.path(dir.r, file_name)
  if (!file.exists(full_path)) {
    stop(paste("File missing:", full_path))
  }
  t(readRDS(full_path))
}

# Pre-allocate arrays (17 horizons x 47 prefectures x 9 methods)
KLD_F <- array(NA, dim = c(n_horizons, n_pref, n_methods), 
               dimnames = list(1:n_horizons, 1:n_pref, method_names))

KLD_M <- array(NA, dim = c(n_horizons, n_pref, n_methods), 
               dimnames = list(1:n_horizons, 1:n_pref, method_names))

# --- Populate Female Array ---
KLD_F[,,1] <- load_t("point_fore_subnational_err_F_EVR_ETS.rds")
KLD_F[,,2] <- load_t("point_fore_subnational_err_MFTS_F_EVR_ETS.rds")
KLD_F[,,3] <- load_t("point_fore_subnational_err_MLFTS_F_EVR_ETS.rds")
KLD_F[,,4] <- load_t1("FANOVA_errors_KLD_female_EVR.rds") 
KLD_F[,,5] <- load_t("point_fore_subnational_err_HDFPCA_F_CDF.rds")
KLD_F[,,6] <- load_t("point_fore_subnational_err_F_K6_ETS.rds")
KLD_F[,,7] <- load_t("point_fore_subnational_err_MFTS_F_K6_ETS.rds")
KLD_F[,,8] <- load_t("point_fore_subnational_err_MLFTS_F_K6_ETS.rds")
KLD_F[,,9] <- load_t1("FANOVA_errors_KLD_female_K.rds")  

# --- Populate Male Array ---
KLD_M[,,1] <- load_t("point_fore_subnational_err_M_EVR_ETS.rds")
KLD_M[,,2] <- load_t("point_fore_subnational_err_MFTS_M_EVR_ETS.rds")
KLD_M[,,3] <- load_t("point_fore_subnational_err_MLFTS_M_EVR_ETS.rds")
KLD_M[,,4] <- load_t1("FANOVA_errors_KLD_male_EVR.rds")
KLD_M[,,5] <- load_t("point_fore_subnational_err_HDFPCA_M_CDF.rds")
KLD_M[,,6] <- load_t("point_fore_subnational_err_M_K6_ETS.rds")
KLD_M[,,7] <- load_t("point_fore_subnational_err_MFTS_M_K6_ETS.rds")
KLD_M[,,8] <- load_t("point_fore_subnational_err_MLFTS_M_K6_ETS.rds")
KLD_M[,,9] <- load_t1("FANOVA_errors_KLD_male_K.rds")

# ==============================================================================
# 3. MCS COMPUTATION 
# ==============================================================================
compute_MCS_survival <- function(loss_array) {
  surv_mat <- matrix(0, nrow = n_horizons, ncol = n_methods,
                     dimnames = list(1:n_horizons, method_names))
  
  for (h in 1:n_horizons) {
    Loss <- matrix(loss_array[h,,], nrow = n_pref, ncol = n_methods)
    colnames(Loss) <- method_names
    
    res <- MCSprocedure(Loss, alpha = alpha, B = B, statistic = statistic, verbose = FALSE)

      surv <- rownames(res@show)[rownames(res@show) %in% res@Info$included]
      surv_mat[h, surv] <- 47 
   
  }
  
  return(surv_mat)
}

# Compute MCS survival counts
counts_female <- compute_MCS_survival(KLD_F)
counts_male   <- compute_MCS_survival(KLD_M)
# ==============================================================================
# 4. PLOTTING FUNCTION (MATCHING EXACT PAPER STYLE)
# ==============================================================================
plot_MCS_final <- function(counts_mat, title_label = "") {
  
  method_labels <- c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", "HDFPCA", 
                     "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)")
  
  df_long <- as.data.frame(counts_mat) %>%
    tibble::rownames_to_column("Horizon") %>%
    pivot_longer(-Horizon, names_to = "Model", values_to = "count") %>%
    mutate(
      Horizon = factor(as.integer(Horizon), levels = 1:n_horizons),
      Model   = factor(Model, levels = method_names)
    )
  
  p <- ggplot(df_long, aes(x = Model, y = factor(Horizon, levels = rev(1:n_horizons)))) +
    geom_tile(aes(fill = count), color = "white", linewidth = 0.4) +
    geom_text(aes(label = count, color = count > 20), 
              size = 3.5, fontface = "bold") +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.1, end = 0.9, limits = c(0, 47)) +
    scale_color_manual(values = c("black", "white"), guide = "none") + 
    scale_x_discrete(labels = method_labels, position = "bottom") +
    labs(y = "Forecast Horizon", x = NULL, fill = "Count", title = title_label) +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 9, color = "black"),
      axis.title.y    = element_text(size = 10, face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = 9, face = "bold"),
      plot.title      = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  return(p)
}


# ==============================================================================
# 5. GENERATE AND SAVE FIGURES
# ==============================================================================
fig_6a <- plot_MCS_final(counts_female, "(a) Female data")
fig_6b <- plot_MCS_final(counts_male, "(b) Male data")

ggsave(file.path(dir.p, "Figure6a_Female_MCS.pdf"), plot = fig_6a, width = 7.5, height = 7)
ggsave(file.path(dir.p, "Figure6b_Male_MCS.pdf"), plot = fig_6b, width = 7.5, height = 7)