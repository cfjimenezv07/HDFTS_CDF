# The rcode below produces Figure 5. Heatmaps for point forecast results.
source("load_packages.R")

dir.results <- here("results", "Results_Point")
dir.plots   <- here("results", "Plots")
if(!dir.exists(dir.plots)) dir.create(dir.plots, recursive = TRUE)

# --- 2. SAFE DATA LOADER ---
load_res <- function(filename) {
  path <- file.path(dir.results, filename)
  if(!file.exists(path)) stop(paste("File missing:", path))
  
  obj <- readRDS(path)
  
  # Handle 3D arrays (take first slice if 3D)
  if (length(dim(obj)) == 3) {
    obj <- obj[, , 1]
  }
  
  obj <- as.matrix(obj)
  
  if (nrow(obj) == 47 && ncol(obj) == 17) {
    obj <- t(obj)
  }
  
  return(obj)
}

# --- 3. CORE COMPUTATION ENGINE ---
compute_min_counts <- function(method_list) {
  n_methods <- length(method_list)
  n_hor     <- nrow(method_list[[1]]) # 17
  n_sub     <- ncol(method_list[[1]]) # 47
  
  # Stack into [subnational, horizon, method] array
  arr <- array(NA, dim = c(n_sub, n_hor, n_methods),
               dimnames = list(NULL, paste0("H", 1:n_hor), names(method_list)))
  
  for (m in names(method_list)) {
    # Transpose [17 x 47] matrix to fill [47 subnationals x 17 horizons]
    arr[, , m] <- t(method_list[[m]])
  }
  
  # Calculate minimum winner counts per horizon
  counts <- matrix(0, nrow = n_hor, ncol = n_methods,
                   dimnames = list(1:n_hor, names(method_list)))
  
  for (h in 1:n_hor) {
    for (s in 1:n_sub) {
      m_idx <- which.min(arr[s, h, ])
      counts[h, m_idx] <- counts[h, m_idx] + 1
    }
  }
  
  # Format into ggplot long-dataframe
  df <- as.data.frame(counts) %>%
    mutate(Horizon = factor(1:n_hor, levels = 1:n_hor)) %>%
    pivot_longer(-Horizon, names_to = "Model", values_to = "Count") %>%
    mutate(Model = factor(Model, levels = names(method_list)))
  
  return(df)
}

# --- 4. GGPLOT HEATMAP RENDERER ---
plot_heatmap <- function(data, custom_labels) {
  ggplot(data, aes(x = Model, y = Horizon)) +
    geom_tile(aes(fill = Count), color = "white", linewidth = 0.3) + 
    geom_text(aes(label = Count, color = Count > (max(Count) / 1.8)), 
              size = 3.8, fontface = "bold") +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.1, end = 0.9) +
    scale_color_manual(values = c("black", "white"), guide = "none") + 
    scale_x_discrete(labels = custom_labels) +
    labs(y = "Forecast Horizon", x = NULL, fill = "Count") +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 10, color = "black"),
      axis.title.y     = element_text(size = 12, face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = 10, face = "bold")
    )
}

# --- 5. DATA LOADING & EXECUTION ---
point_female_list <- list(
  UFTS_EVR   = load_res("point_fore_subnational_err_F_EVR_ETS.rds"),
  MFTS_EVR   = load_res("point_fore_subnational_err_MFTS_F_EVR_ETS.rds"),
  MLFTS_EVR  = load_res("point_fore_subnational_err_MLFTS_F_EVR_ETS.rds"),
  FANOVA_EVR = load_res("FANOVA_errors_KLD_female_EVR.rds"),
  HDFPCA     = load_res("point_fore_subnational_err_HDFPCA_F_CDF.rds"),
  UFTS_K6    = load_res("point_fore_subnational_err_F_K6_ETS.rds"),
  MFTS_K6    = load_res("point_fore_subnational_err_MFTS_F_K6_ETS.rds"),
  MLFTS_K6   = load_res("point_fore_subnational_err_MLFTS_F_K6_ETS.rds"),
  FANOVA_K6  = load_res("FANOVA_errors_KLD_female_K.rds")
)

point_male_list <- list(
  UFTS_EVR   = load_res("point_fore_subnational_err_M_EVR_ETS.rds"),
  MFTS_EVR   = load_res("point_fore_subnational_err_MFTS_M_EVR_ETS.rds"),
  MLFTS_EVR  = load_res("point_fore_subnational_err_MLFTS_M_EVR_ETS.rds"),
  FANOVA_EVR = load_res("FANOVA_errors_KLD_male_EVR.rds"),
  HDFPCA     = load_res("point_fore_subnational_err_HDFPCA_M_CDF.rds"),
  UFTS_K6    = load_res("point_fore_subnational_err_M_K6_ETS.rds"),
  MFTS_K6    = load_res("point_fore_subnational_err_MFTS_M_K6_ETS.rds"),
  MLFTS_K6   = load_res("point_fore_subnational_err_MLFTS_M_K6_ETS.rds"),
  FANOVA_K6  = load_res("FANOVA_errors_KLD_male_K.rds")
)

labels_9m <- c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", 
               "HDFPCA", 
               "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)")

# Compute dataframes
df_p_female <- compute_min_counts(point_female_list)
df_p_male   <- compute_min_counts(point_male_list)

# Generate plots
fig_female <- plot_heatmap(df_p_female, labels_9m)
fig_male   <- plot_heatmap(df_p_male, labels_9m)

# Save PDF figures
ggsave(file.path(dir.plots, "Figure5a_Female.pdf"), plot = fig_female, width = 8, height = 9)
ggsave(file.path(dir.plots, "Figure5b_Male.pdf"),   plot = fig_male,   width = 8, height = 9)
