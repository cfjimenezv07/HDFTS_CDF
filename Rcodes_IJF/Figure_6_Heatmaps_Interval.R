# The rcode below produces Figure 5. Heatmaps for interval forecast results.

source("load_packages.R")


dir.results <- here( "results", "Results_Interval", "sd")
dir.plots   <- here( "results", "Plots")

if (!dir.exists(dir.plots)) dir.create(dir.plots, recursive = TRUE)
if(!dir.exists(dir.plots)) dir.create(dir.plots, recursive = TRUE)


# --- 2. SAFE INTERVAL DATA LOADER ---
# Loads .rds files, extracts metric_slice (2 = CPD, 3 = Interval Score), 
# and guarantees output is formatted as [16 horizons x 47 subnationals]
load_interval_res <- function(filename, metric_slice = 2) {
  path <- file.path(dir.results, filename)
  if(!file.exists(path)) stop(paste("File missing:", path))
  
  obj <- readRDS(path)
  
  # Ensure input is a 3D array [Dim1, Dim2, Dim3]
  if (length(dim(obj)) != 3) {
    stop(paste("Expected 3D array for interval results in:", filename))
  }
  
  # Extract metric slice depending on array orientation:
  # Case A: [47 subnationals, 16 horizons, 3 metrics]
  if (dim(obj)[1] == 47) {
    mat <- t(obj[, , metric_slice])
  } 
  # Case B: [16 horizons, 3 metrics, 47 subnationals]
  else if (dim(obj)[2] == 3) {
    mat <- obj[, metric_slice, ]
  } 
  # Fallback for alternative orderings
  else {
    mat <- obj[, , metric_slice]
    if (nrow(mat) == 47 && ncol(mat) == 16) {
      mat <- t(mat)
    }
  }
  
  return(as.matrix(mat))
}

# --- 3. CORE COMPUTATION ENGINE ---
compute_min_counts <- function(method_list) {
  n_methods <- length(method_list)
  n_hor     <- nrow(method_list[[1]]) # 16
  n_sub     <- ncol(method_list[[1]]) # 47
  
  # Stack into [subnational, horizon, method] array
  arr <- array(NA, dim = c(n_sub, n_hor, n_methods),
               dimnames = list(NULL, paste0("H", 1:n_hor), names(method_list)))
  
  for (m in names(method_list)) {
    # Transpose [16 x 47] matrix to fill [47 subnationals x 16 horizons]
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

# --- 5. MODEL MAPPING & EXECUTION ---
labels_9m <- c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", 
               "HDFPCA", 
               "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)")

generate_interval_lists <- function(metric_idx) {
  list(
    female = list(
      UFTS_EVR   = load_interval_res("int_fore_subnational_err_F_EVR_ETS.rds", metric_idx),
      MFTS_EVR   = load_interval_res("MFTS_int_fore_subnational_err_F_EVR_ETS.rds", metric_idx),
      MLFTS_EVR  = load_interval_res("MLFTS_int_fore_subnational_err_F_EVR_ETS.rds", metric_idx),
      FANOVA_EVR = load_interval_res("int_fore_subnational_err_F_EVR_ETS_FANOVA.rds", metric_idx),
      HDFPCA     = load_interval_res("hdfpca_int_fore_subnational_err_F_EVR_ETS_CDF.rds", metric_idx),
      UFTS_K6    = load_interval_res("int_fore_subnational_err_F_K6_ETS.rds", metric_idx),
      MFTS_K6    = load_interval_res("MFTS_int_fore_subnational_err_F_K6_ETS.rds", metric_idx),
      MLFTS_K6   = load_interval_res("MLFTS_int_fore_subnational_err_F_K6_ETS.rds", metric_idx),
      FANOVA_K6  = load_interval_res("int_fore_subnational_err_F_ETS_K6_FANOVA.rds", metric_idx)
    ),
    male = list(
      UFTS_EVR   = load_interval_res("int_fore_subnational_err_M_EVR_ETS.rds", metric_idx),
      MFTS_EVR   = load_interval_res("MFTS_int_fore_subnational_err_M_EVR_ETS.rds", metric_idx),
      MLFTS_EVR  = load_interval_res("MLFTS_int_fore_subnational_err_M_EVR_ETS.rds", metric_idx),
      FANOVA_EVR = load_interval_res("int_fore_subnational_err_M_EVR_ETS_FANOVA.rds", metric_idx),
      HDFPCA     = load_interval_res("hdfpca_int_fore_subnational_err_M_EVR_ETS_CDF.rds", metric_idx),
      UFTS_K6    = load_interval_res("int_fore_subnational_err_M_K6_ETS.rds", metric_idx),
      MFTS_K6    = load_interval_res("MFTS_int_fore_subnational_err_M_K6_ETS.rds", metric_idx),
      MLFTS_K6   = load_interval_res("MLFTS_int_fore_subnational_err_M_K6_ETS.rds", metric_idx),
      FANOVA_K6  = load_interval_res("int_fore_subnational_err_M_ETS_K6_FANOVA.rds", metric_idx)
    )
  )
}

# --- 6. RUN AND SAVE FIGURES ---

# --- CPD (Metric Index 2) ---
cpd_lists  <- generate_interval_lists(metric_idx = 2)
df_cpd_f   <- compute_min_counts(cpd_lists$female)
df_cpd_m   <- compute_min_counts(cpd_lists$male)

ggsave(file.path(dir.plots, "Figure5a_Female_CPD.pdf"), plot = plot_heatmap(df_cpd_f, labels_9m), width = 8, height = 9)
ggsave(file.path(dir.plots, "Figure5b_Male_CPD.pdf"),   plot = plot_heatmap(df_cpd_m, labels_9m), width = 8, height = 9)

# --- Interval Score (Metric Index 3) ---
score_lists <- generate_interval_lists(metric_idx = 3)
df_score_f  <- compute_min_counts(score_lists$female)
df_score_m  <- compute_min_counts(score_lists$male)

ggsave(file.path(dir.plots, "Figure5c_Female_Score.pdf"), plot = plot_heatmap(df_score_f, labels_9m), width = 8, height = 9)
ggsave(file.path(dir.plots, "Figure5d_Male_Score.pdf"),   plot = plot_heatmap(df_score_m, labels_9m), width = 8, height = 9)