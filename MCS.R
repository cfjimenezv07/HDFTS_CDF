# HLS where the HLS results are stored
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/HLS/"
# FANOVA where the FANOVA results are stored
dir.r1 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/FANOVA/"

# Where heatmaps will be saved
dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/Plots/"

# source("Heatmaps.R")

library(MCS)
library(tidyverse)
library(ggplot2)

# =========================
# PARAMETERS
# =========================
n_methods <- 9
n_horizons <- 17
n_pref <- 47

method_names <- c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA",
                  "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)", "FANOVA(K=6)")

alpha <- 0.10
B <- 2000
statistic <- "Tmax"



# =========================
# ARRAYS: KLD
# =========================
# UFTS
point_fore_subnational_err_F_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_F_EVR_ETS.rds"))[,,1])
point_fore_subnational_err_M_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_M_EVR_ETS.rds"))[,,1])
point_fore_subnational_err_F_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_F_K6_ETS.rds"))[,,1])
point_fore_subnational_err_M_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_M_K6_ETS.rds"))[,,1])

# MFTS
point_fore_subnational_err_MFTS_F_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_F_EVR_ETS.rds"))[,,1])
point_fore_subnational_err_MFTS_M_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_M_EVR_ETS.rds"))[,,1])
point_fore_subnational_err_MFTS_F_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_F_K6_ETS.rds"))[,,1])
point_fore_subnational_err_MFTS_M_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_M_K6_ETS.rds"))[,,1])

# MLFTS
point_fore_subnational_err_MLFTS_F_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_EVR_ETS.rds"))[,,1])
point_fore_subnational_err_MLFTS_M_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_EVR_ETS.rds"))[,,1])
point_fore_subnational_err_MLFTS_F_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_K6_ETS.rds"))[,,1])
point_fore_subnational_err_MLFTS_M_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_K6_ETS.rds"))[,,1])

# HDFPCA
point_fore_subnational_err_HDFPCA_F_CDF <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF.rds"))[,,1])
point_fore_subnational_err_HDFPCA_M_CDF <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF.rds"))[,,1])


# Female array
KLD_fore_subnational_err_F_array <- array(NA, dim=c(n_horizons, n_pref, n_methods),
                                          dimnames=list(1:n_horizons, 1:n_pref, method_names))

KLD_fore_subnational_err_F_array[,,1] <- point_fore_subnational_err_F_EVR_ETS
KLD_fore_subnational_err_F_array[,,2] <- point_fore_subnational_err_MFTS_F_EVR_ETS
KLD_fore_subnational_err_F_array[,,3] <- point_fore_subnational_err_MLFTS_F_EVR_ETS
KLD_fore_subnational_err_F_array[,,4] <- errors_KLD_female_EVR
KLD_fore_subnational_err_F_array[,,5] <- point_fore_subnational_err_HDFPCA_F_CDF
KLD_fore_subnational_err_F_array[,,6] <- point_fore_subnational_err_F_K6_ETS
KLD_fore_subnational_err_F_array[,,7] <- point_fore_subnational_err_MFTS_F_K6_ETS
KLD_fore_subnational_err_F_array[,,8] <- point_fore_subnational_err_MLFTS_F_K6_ETS
KLD_fore_subnational_err_F_array[,,9] <- errors_KLD_female_K

# Male array
KLD_fore_subnational_err_M_array <- array(NA, dim=c(n_horizons, n_pref, n_methods),
                                          dimnames=list(1:n_horizons, 1:n_pref, method_names))

KLD_fore_subnational_err_M_array[,,1] <- point_fore_subnational_err_M_EVR_ETS
KLD_fore_subnational_err_M_array[,,2] <- point_fore_subnational_err_MFTS_M_EVR_ETS
KLD_fore_subnational_err_M_array[,,3] <- point_fore_subnational_err_MLFTS_M_EVR_ETS
KLD_fore_subnational_err_M_array[,,4] <- errors_KLD_male_EVR
KLD_fore_subnational_err_M_array[,,5] <- point_fore_subnational_err_HDFPCA_M_CDF
KLD_fore_subnational_err_M_array[,,6] <- point_fore_subnational_err_M_K6_ETS
KLD_fore_subnational_err_M_array[,,7] <- point_fore_subnational_err_MFTS_M_K6_ETS
KLD_fore_subnational_err_M_array[,,8] <- point_fore_subnational_err_MLFTS_M_K6_ETS
KLD_fore_subnational_err_M_array[,,9] <- errors_KLD_male_K


# =========================
# ARRAYS: JSD
# =========================

# UFTS
point_fore_subnational_err_F_JSD_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_F_EVR_ETS.rds"))[,,2])
point_fore_subnational_err_M_JSD_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_M_EVR_ETS.rds"))[,,2])
point_fore_subnational_err_F_JSD_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_F_K6_ETS.rds"))[,,2])
point_fore_subnational_err_M_JSD_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_M_K6_ETS.rds"))[,,2])

# MFTS
point_fore_subnational_err_MFTS_F_JSD_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_F_EVR_ETS.rds"))[,,2])
point_fore_subnational_err_MFTS_M_JSD_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_M_EVR_ETS.rds"))[,,2])
point_fore_subnational_err_MFTS_F_JSD_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_F_K6_ETS.rds"))[,,2])
point_fore_subnational_err_MFTS_M_JSD_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MFTS_M_K6_ETS.rds"))[,,2])

# MLFTS
point_fore_subnational_err_MLFTS_F_JSD_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_EVR_ETS.rds"))[,,2])
point_fore_subnational_err_MLFTS_M_JSD_EVR_ETS <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_EVR_ETS.rds"))[,,2])
point_fore_subnational_err_MLFTS_F_JSD_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_F_K6_ETS.rds"))[,,2])
point_fore_subnational_err_MLFTS_M_JSD_K6_ETS  <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_MLFTS_M_K6_ETS.rds"))[,,2])

# HDFPCA
point_fore_subnational_err_HDFPCA_F_JSD_CDF <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_F_CDF.rds"))[,,2])
point_fore_subnational_err_HDFPCA_M_JSD_CDF <- t(readRDS(file = paste0(dir.r, "point_fore_subnational_err_HDFPCA_M_CDF.rds"))[,,2])

# =========================
# Female array (_JSD)
# =========================
JSD_fore_subnational_err_F_array <- array(NA, dim=c(n_horizons, n_pref, n_methods),
                                          dimnames=list(1:n_horizons, 1:n_pref, method_names))

JSD_fore_subnational_err_F_array[,,1] <- point_fore_subnational_err_F_JSD_EVR_ETS
JSD_fore_subnational_err_F_array[,,2] <- point_fore_subnational_err_MFTS_F_JSD_EVR_ETS
JSD_fore_subnational_err_F_array[,,3] <- point_fore_subnational_err_MLFTS_F_JSD_EVR_ETS
JSD_fore_subnational_err_F_array[,,4] <- errors_JSD_female_EVR
JSD_fore_subnational_err_F_array[,,5] <- point_fore_subnational_err_HDFPCA_F_JSD_CDF
JSD_fore_subnational_err_F_array[,,6] <- point_fore_subnational_err_F_JSD_K6_ETS
JSD_fore_subnational_err_F_array[,,7] <- point_fore_subnational_err_MFTS_F_JSD_K6_ETS
JSD_fore_subnational_err_F_array[,,8] <- point_fore_subnational_err_MLFTS_F_JSD_K6_ETS
JSD_fore_subnational_err_F_array[,,9] <- errors_JSD_female_K

# =========================
# Male array (_JSD)
# =========================
JSD_fore_subnational_err_M_array <- array(NA, dim=c(n_horizons, n_pref, n_methods),
                                          dimnames=list(1:n_horizons, 1:n_pref, method_names))

JSD_fore_subnational_err_M_array[,,1] <- point_fore_subnational_err_M_JSD_EVR_ETS
JSD_fore_subnational_err_M_array[,,2] <- point_fore_subnational_err_MFTS_M_JSD_EVR_ETS
JSD_fore_subnational_err_M_array[,,3] <- point_fore_subnational_err_MLFTS_M_JSD_EVR_ETS
JSD_fore_subnational_err_M_array[,,4] <- errors_JSD_male_EVR
JSD_fore_subnational_err_M_array[,,5] <- point_fore_subnational_err_HDFPCA_M_JSD_CDF
JSD_fore_subnational_err_M_array[,,6] <- point_fore_subnational_err_M_JSD_K6_ETS
JSD_fore_subnational_err_M_array[,,7] <- point_fore_subnational_err_MFTS_M_JSD_K6_ETS
JSD_fore_subnational_err_M_array[,,8] <- point_fore_subnational_err_MLFTS_M_JSD_K6_ETS
JSD_fore_subnational_err_M_array[,,9] <- errors_JSD_male_K


# =========================
# FUNCTION: Compute counts per horizon
# =========================
compute_MCS_counts_per_horizon <- function(loss_array, alpha=0.1, B=2000, statistic="Tmax") {
  
  n_horizons <- dim(loss_array)[1]
  n_pref <- dim(loss_array)[2]
  n_methods <- dim(loss_array)[3]
  method_names <- dimnames(loss_array)[[3]]
  
  counts_mat <- matrix(0, nrow=n_horizons, ncol=n_methods,
                       dimnames=list(1:n_horizons, method_names))
  
  for(h in 1:n_horizons) {
    # Take all prefectures at this horizon
    Loss <- matrix(loss_array[h,,], nrow=n_pref, ncol=n_methods)
    colnames(Loss) <- method_names
    
    # Run MCS
    res <- tryCatch({
      MCSprocedure(Loss, alpha=alpha, B=B, statistic=statistic, verbose=FALSE)
    }, error=function(e) NULL)
    
    if(!is.null(res)) {
      surv <- rownames(res@show)
      # Count number of prefectures contributing to each method in superior set
      counts_mat[h, surv] <- colSums(!is.na(Loss[, surv]))
    }
  }
  return(counts_mat)
}

# =========================
# RUN MCS
# =========================
counts_female <- compute_MCS_counts_per_horizon(KLD_fore_subnational_err_F_array, alpha, B, statistic)
counts_male   <- compute_MCS_counts_per_horizon(KLD_fore_subnational_err_M_array, alpha, B, statistic)

counts_female_JSD <- compute_MCS_counts_per_horizon(JSD_fore_subnational_err_F_array, alpha, B, statistic)
counts_male_JSD   <- compute_MCS_counts_per_horizon(JSD_fore_subnational_err_M_array, alpha, B, statistic)

# =========================
# FUNCTION: Plot heatmap 
# =========================
plot_MCS_heatmap <- function(counts_mat, file_name, title_text) {
  
  n_horizons <- nrow(counts_mat)
  
  # 9 methods in order
  all_methods <- c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA",
                   "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)", "FANOVA(K=6)")
  
  method_labels <- c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", "HDFPCA\n",
                     "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)")
  
  # Ensure all 9 columns exist
  missing_methods <- setdiff(all_methods, colnames(counts_mat))
  if(length(missing_methods) > 0){
    counts_mat <- cbind(counts_mat, matrix(0, nrow=n_horizons, ncol=length(missing_methods),
                                           dimnames=list(NULL, missing_methods)))
  }
  counts_mat <- counts_mat[, all_methods]  # reorder
  
  # Pivot to long format
  counts_long <- counts_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Horizon") %>%
    pivot_longer(-Horizon, names_to="Model", values_to="count") %>%
    mutate(
      Horizon = factor(as.integer(Horizon), levels=1:n_horizons),
      Model = factor(Model, levels=all_methods)  # FIX: ensure all 9 levels
    )
  
  # Plot
  p <- ggplot(counts_long, aes(x=Model, y=factor(Horizon, levels=rev(1:n_horizons)))) +
    geom_tile(aes(fill=count), color="black") +
    geom_text(aes(label=count), size=3) +
    scale_fill_gradient2(high="green", mid="white", low="red",
                         midpoint = 5, na.value="yellow") +
    scale_x_discrete(labels=method_labels, position="bottom") +
    ylab("Forecast horizon") +
    ggtitle(title_text) +
    theme_bw() +
    theme(
      legend.position="right",
      axis.text.x = element_text(angle=0, vjust=0.5)
    )
  
  # Save
  ggsave(file.path(dir.p, file_name), plot=p, width=10, height=6)
}








# =========================
# CREATE HEATMAPS
# =========================
plot_MCS_heatmap(counts_female, "Fig_6a.png", "Female Forecasts")
plot_MCS_heatmap(counts_male, "Fig_6b.png", "Male Forecasts")
