# This Rcode produces the Figure 7 in the paper about an example of the prediction intervals computed.
# --------------------------------------------------------------------------------------------------------------
source("load_packages.R")

# Subdirectories
dir.data    <- here("data", "updated_data")
dir.aux     <- here("auxiliary_source")
dir.results <- here("results", "Results_Interval")
dir.plots   <- here("results", "Plots")

# Create output directories if they do not exist
if(!dir.exists(dir.plots)) dir.create(dir.plots, recursive = TRUE)
if(!dir.exists(dir.results)) dir.create(dir.results, recursive = TRUE)

# Load external dependencies
source(file.path("load_packages.R"))
source(file.path(dir.aux, "auxiliary_interval.R"))
source(file.path(dir.aux, "CDF_transformation.R"))

# 2. DATA LOADING & PREPARATION
# ------------------------------------------------------------------------------
year          <- 1975:2023
n_year        <- length(year)
age           <- 0:110
n_age         <- length(age)
n_prefectures <- 47
n_populations <- 2
new_age       <- 0:109
nn_age        <- length(new_age)

state <- c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
           "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
           "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
           "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
           "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
           "Tokushima", "Kagawa", "Ehime", "Kochi",
           "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

part_list   <- lapply(1:n_prefectures, function(ik) (n_year*ik-(n_year-1)):(n_year*ik))
part_list_c <- lapply(1:n_populations, function(ik) (n_age*ik-(n_age-1)):(n_age*ik))

female_prefecture_qx <- male_prefecture_qx <- total_prefecture_qx <- list()
for (ik in seq_along(state)) {
  female_data <- read.table(paste0(dir.data, "female_prefecture_", ik, ".txt"), header = TRUE, skip = 2)
  male_data   <- read.table(paste0(dir.data, "male_prefecture_", ik, ".txt"), header = TRUE, skip = 2)
  total_data  <- read.table(paste0(dir.data, "total_prefecture_", ik, ".txt"), header = TRUE, skip = 2)
  
  female_prefecture_qx[[ik]] <- t(matrix(female_data[female_data$Year >= 1975, ]$qx, n_age, n_year))
  male_prefecture_qx[[ik]]   <- t(matrix(male_data[male_data$Year >= 1975, ]$qx, n_age, n_year))
  total_prefecture_qx[[ik]]  <- t(matrix(total_data[total_data$Year >= 1975, ]$qx, n_age, n_year))
}

female_prefecture_dx <- male_prefecture_dx <- total_prefecture_dx <- list()
for(iw in 1:length(state)) {
  f_dum = m_dum = t_dum = matrix(NA, n_year, n_age)
  for(ij in 1:n_year) {
    s_f = s_m = s_t = 10^5
    for(ik in 1:n_age) {
      f_dum[ij,ik] = female_prefecture_qx[[iw]][ij,ik] * s_f
      s_f = s_f - f_dum[ij,ik]
      m_dum[ij,ik] = male_prefecture_qx[[iw]][ij,ik] * s_m
      s_m = s_m - m_dum[ij,ik]
      t_dum[ij,ik] = total_prefecture_qx[[iw]][ij,ik] * s_t
      s_t = s_t - t_dum[ij,ik]
    }
  }
  female_prefecture_dx[[iw]] = t(f_dum)
  male_prefecture_dx[[iw]]   = t(m_dum)
  total_prefecture_dx[[iw]]  = t(t_dum)
}

# 3. HELPER FUNCTIONS & SHARED PLOTTING THEME
# ------------------------------------------------------------------------------
tune_para_find_function <- function(tune_para, resi_mat, sd_val_input, alpha_level, PI_type) {
  n_age_loc = nrow(resi_mat)
  if(PI_type == "pointwise") {
    ind = matrix(NA, n_age_loc, ncol(resi_mat))
    for(iw in 1:ncol(resi_mat)) { ind[,iw] = ifelse(resi_mat[,iw] >= -tune_para * sd_val_input & resi_mat[,iw] <= tune_para * sd_val_input, 1, 0) }
    ecp = sum(ind)/(n_age_loc * ncol(resi_mat))
  }
  return(abs(ecp - alpha_level))
}

plot_and_save <- function(df_h, filename) {
  p <- ggplot(df_h, aes(x = Age, group = Gender)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Gender), alpha = 0.2) +
    geom_line(aes(y = Forecast, color = Gender, linetype = Gender), linewidth = 1) +
    scale_fill_manual(values = c("Female" = "firebrick3", "Male" = "dodgerblue4")) +
    scale_color_manual(values = c("Female" = "firebrick3", "Male" = "dodgerblue4")) +
    scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +
    theme_bw(base_size = 15) +
    labs(x = "Age", y = "Forecasted dx", title = NULL, subtitle = NULL) +
    theme(
      legend.position = c(0.05, 0.95), 
      legend.justification = c("left", "top"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
      legend.key = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  
  print(p)
  ggsave(filename = paste0(dir.plots, filename), plot = p, width = 8, height = 6)
}

# 4. METHOD 1: UFTS
# ------------------------------------------------------------------------------
# Pass fdata as (Years x Ages) inside functions
get_multi_horizon_data_UFTS <- function(fdata_mat, horizons, method_ncomp, fore_method, uni_fore_method, level_sig, type_PI) {
  # Transpose input matrix if needed so rows are years and columns are ages
  if(nrow(fdata_mat) == 111) fdata <- t(fdata_mat) else fdata <- fdata_mat
  
  all_horizons_df <- data.frame()
  for(h in horizons) {
    fore_validation = matrix(NA, ncol(fdata), (18 - h))
    for(ij in 1:(18 - h)) {
      if(fore_method == "CDF") {
        fore_validation[,ij] = fore_national_cdf(data_set = fdata[1:(15 + ij),], ncomp_method = method_ncomp, fh = h, fmethod = uni_fore_method)
      } else {
        fore_validation[,ij] = as.numeric(clr_fun(fdata = fdata[1:(15 + ij),], ncomp_selection = method_ncomp, fh = h, fore_method = "ETS")$fore_count)
      }
    }
    
    holdout_validation_dum = t(matrix(fdata[(16 + h):33,], length((16 + h):33), ncol(fdata)))
    resi_mat = holdout_validation_dum - fore_validation
    sd_val_input = apply(resi_mat, 1, sd)
    tune_para = optimise(f = tune_para_find_function, interval = c(0, 20), resi_mat = resi_mat, sd_val_input = sd_val_input, alpha_level = level_sig, PI_type = type_PI)$minimum
    
    fore_val <- if(fore_method == "CDF") fore_national_cdf(data_set = fdata[1:32,], ncomp_method = method_ncomp, fh = h, fmethod = uni_fore_method) else as.numeric(clr_fun(fdata = fdata[1:32,], ncomp_selection = method_ncomp, fh = h, fore_method = "ETS")$fore_count)
    
    temp_df <- data.frame(
      Age = 0:(ncol(fdata)-1),
      Forecast = fore_val,
      Lower = fore_val - tune_para * sd_val_input,
      Upper = fore_val + tune_para * sd_val_input,
      Horizon = paste0("h = ", h)
    )
    all_horizons_df <- rbind(all_horizons_df, temp_df)
  }
  return(all_horizons_df)
}

ufts_female <- get_multi_horizon_data_UFTS(female_prefecture_dx[[1]], 1:3, "EVR", "CDF", "ets", 0.80, "pointwise") %>% mutate(Gender = "Female")
ufts_male   <- get_multi_horizon_data_UFTS(male_prefecture_dx[[1]], 1:3, "EVR", "CDF", "ets", 0.80, "pointwise") %>% mutate(Gender = "Male")
ufts_plot_data <- rbind(ufts_female, ufts_male)

for(h_val in 1:3) {
  df_h <- ufts_plot_data %>% filter(Horizon == paste0("h = ", h_val))
  plot_and_save(df_h, paste0("UFTS_h", h_val, ".pdf"))
}

# 5. METHOD 2: MFTS
# ------------------------------------------------------------------------------
get_multi_horizon_data_MFTS <- function(fdata_F_mat, fdata_M_mat, horizons, method_ncomp, fore_method, uni_fore_method, level_sig, type_PI) {
  fdata_F <- if(nrow(fdata_F_mat) == 111) t(fdata_F_mat) else fdata_F_mat
  fdata_M <- if(nrow(fdata_M_mat) == 111) t(fdata_M_mat) else fdata_M_mat
  
  all_horizons_df <- data.frame()
  for(h in horizons) {
    forecast_validation_F = forecast_validation_M = matrix(NA, ncol(fdata_F), (18 - h))
    for(ij in 1:(18 - h)) {
      dum <- fore_national_cdf_MFTS(data_set_F = fdata_F[1:(15+ij),], data_set_M = fdata_M[1:(15+ij),], fh = h, fmethod = uni_fore_method, method_ncomp = method_ncomp)
      forecast_validation_F[,ij] = dum$mfts_fore_F
      forecast_validation_M[,ij] = dum$mfts_fore_M
    }
    
    holdout_F = t(matrix(fdata_F[(16 + h):33,], length((16 + h):33), ncol(fdata_F)))
    holdout_M = t(matrix(fdata_M[(16 + h):33,], length((16 + h):33), ncol(fdata_M)))
    
    sd_val_F = apply(holdout_F - forecast_validation_F, 1, sd)
    sd_val_M = apply(holdout_M - forecast_validation_M, 1, sd)
    
    tune_F = optimise(f = tune_para_find_function, interval = c(0, 20), resi_mat = holdout_F - forecast_validation_F, sd_val_input = sd_val_F, alpha_level = level_sig, PI_type = type_PI)$minimum
    tune_M = optimise(f = tune_para_find_function, interval = c(0, 20), resi_mat = holdout_M - forecast_validation_M, sd_val_input = sd_val_M, alpha_level = level_sig, PI_type = type_PI)$minimum
    
    final_fore <- fore_national_cdf_MFTS(data_set_F = fdata_F[1:32,], data_set_M = fdata_M[1:32,], fh = h, fmethod = uni_fore_method, method_ncomp = method_ncomp)
    
    temp_F <- data.frame(Age = 0:(ncol(fdata_F)-1), Forecast = final_fore$mfts_fore_F, Lower = final_fore$mfts_fore_F - tune_F * sd_val_F, Upper = final_fore$mfts_fore_F + tune_F * sd_val_F, Horizon = paste0("h = ", h), Gender = "Female")
    temp_M <- data.frame(Age = 0:(ncol(fdata_M)-1), Forecast = final_fore$mfts_fore_M, Lower = final_fore$mfts_fore_M - tune_M * sd_val_M, Upper = final_fore$mfts_fore_M + tune_M * sd_val_M, Horizon = paste0("h = ", h), Gender = "Male")
    
    all_horizons_df <- rbind(all_horizons_df, temp_F, temp_M)
  }
  return(all_horizons_df)
}

mfts_plot_data <- get_multi_horizon_data_MFTS(female_prefecture_dx[[1]], male_prefecture_dx[[1]], 1:3, "EVR", "CDF", "ets", 0.80, "pointwise")

for(h_val in 1:3) {
  df_h <- mfts_plot_data %>% filter(Horizon == paste0("h = ", h_val))
  plot_and_save(df_h, paste0("MFTS_h", h_val, ".pdf"))
}

# 6. METHOD 3: MLFTS
# ------------------------------------------------------------------------------
get_multi_horizon_data_MLFTS <- function(fdata_F_mat, fdata_M_mat, fdata_common_mat, horizons, method_ncomp, fore_method, level_sig, type_PI) {
  fdata_F <- if(nrow(fdata_F_mat) == 111) t(fdata_F_mat) else fdata_F_mat
  fdata_M <- if(nrow(fdata_M_mat) == 111) t(fdata_M_mat) else fdata_M_mat
  fdata_common <- if(!is.null(fdata_common_mat) && nrow(fdata_common_mat) == 111) t(fdata_common_mat) else fdata_common_mat
  
  all_horizons_df <- data.frame()
  for(h in horizons) {
    forecast_validation_F = forecast_validation_M = matrix(NA, ncol(fdata_F), (18 - h))
    for(ij in 1:(18 - h)) {
      dum <- fore_national_cdf_MLFTS(
        data_set_F = fdata_F[1:(15+ij),], 
        data_set_M = fdata_M[1:(15+ij),],
        aux_variable = if(is.null(fdata_common)) NULL else fdata_common[1:(15+ij),], 
        fh = h, fmethod = "ets", method_ncomp = method_ncomp
      )
      forecast_validation_F[,ij] = dum$mlfts_fore_F
      forecast_validation_M[,ij] = dum$mlfts_fore_M
    }
    
    holdout_F = t(matrix(fdata_F[(16 + h):33,], length((16 + h):33), ncol(fdata_F)))
    holdout_M = t(matrix(fdata_M[(16 + h):33,], length((16 + h):33), ncol(fdata_M)))
    
    sd_val_F = apply(holdout_F - forecast_validation_F, 1, sd)
    sd_val_M = apply(holdout_M - forecast_validation_M, 1, sd)
    
    tune_F = optimise(f = tune_para_find_function, interval = c(0, 20), resi_mat = holdout_F - forecast_validation_F, sd_val_input = sd_val_F, alpha_level = level_sig, PI_type = type_PI)$minimum
    tune_M = optimise(f = tune_para_find_function, interval = c(0, 20), resi_mat = holdout_M - forecast_validation_M, sd_val_input = sd_val_M, alpha_level = level_sig, PI_type = type_PI)$minimum
    
    final_fore <- fore_national_cdf_MLFTS(
      data_set_F = fdata_F[1:32,], 
      data_set_M = fdata_M[1:32,], 
      aux_variable = if(is.null(fdata_common)) NULL else fdata_common[1:32,], 
      fh = h, fmethod = "ets", method_ncomp = method_ncomp
    )
    
    df_h <- rbind(
      data.frame(Age = 0:(ncol(fdata_F)-1), Forecast = final_fore$mlfts_fore_F, Lower = final_fore$mlfts_fore_F - tune_F * sd_val_F, Upper = final_fore$mlfts_fore_F + tune_F * sd_val_F, Gender = "Female", Horizon = paste0("h = ", h)),
      data.frame(Age = 0:(ncol(fdata_M)-1), Forecast = final_fore$mlfts_fore_M, Lower = final_fore$mlfts_fore_M - tune_M * sd_val_M, Upper = final_fore$mlfts_fore_M + tune_M * sd_val_M, Gender = "Male", Horizon = paste0("h = ", h))
    )
    all_horizons_df <- rbind(all_horizons_df, df_h)
  }
  return(all_horizons_df)
}

mlfts_plot_data <- get_multi_horizon_data_MLFTS(female_prefecture_dx[[1]], male_prefecture_dx[[1]], NULL, 1:3, "EVR", "CDF", 0.80, "pointwise")

for(h_val in 1:3) {
  df_h <- mlfts_plot_data %>% filter(Horizon == paste0("h = ", h_val))
  plot_and_save(df_h, paste0("MLFTS_Clean_h", h_val, ".pdf"))
}

# 7. METHOD 4: FANOVA
# ------------------------------------------------------------------------------
Transformed_female <- lapply(1:n_prefectures, function(i) t(cdf_transformation(t(female_prefecture_dx[[i]])/10^5, year)))
Transformed_male   <- lapply(1:n_prefectures, function(i) t(cdf_transformation(t(male_prefecture_dx[[i]])/10^5, year)))

all_unconstrained_female <- t(do.call(cbind, Transformed_female))
all_unconstrained_male   <- t(do.call(cbind, Transformed_male))

Residuals_means <- hdftsa::Two_way_mean_residuals(t(all_unconstrained_male), t(all_unconstrained_female), year, new_age, n_prefectures, n_populations)

Res1_means = Residuals_means$residuals1_mean
Res2_means = Residuals_means$residuals2_mean
Fixed_part_means_1 = Residuals_means$Fixed_comp_mean[, 1:nn_age]
Fixed_part_means_2 = Residuals_means$Fixed_comp_mean[, (nn_age+1):(2*nn_age)]

male_prefecture_res_means   <- lapply(1:n_prefectures, function(k) Res1_means[part_list[[k]], ])
female_prefecture_res_means <- lapply(1:n_prefectures, function(k) Res2_means[part_list[[k]], ])
male_prefecture_fixed_means <- lapply(1:n_prefectures, function(k) Fixed_part_means_1[part_list[[k]], ])
female_prefecture_fixed_means <- lapply(1:n_prefectures, function(k) Fixed_part_means_2[part_list[[k]], ])

select_k <- function(tau, eigenvalue) {
  k_all = sapply(1:(length(eigenvalue)-1), function(k) (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0))
  return(which.min(k_all))
}

Pref_forecast_curves <- function(fixed_com, Residuals_f, est_method="cov", fh=30, B=1000, prediction_method="ets", select_K="EVR", K=6) {
  med_polish_resi = t(Residuals_f)
  med_polish_resi_lrc = if(est_method == "cov") cov(t(med_polish_resi)) else long_run_covariance_estimation(med_polish_resi)
  med_polish_resi_eigen = eigen(med_polish_resi_lrc)
  retain_component = if(select_K == "EVR") select_k(10^-2, med_polish_resi_eigen$values) else K
  med_polish_resi_basis = as.matrix(med_polish_resi_eigen$vectors[,1:retain_component])
  med_polish_resi_score = crossprod(med_polish_resi, med_polish_resi_basis)
  
  med_polish_resi_score_forecast = matrix(NA, retain_component, fh)
  for(ip in 1:retain_component) {
    med_polish_resi_score_forecast[ip,] = as.vector(forecast:::forecast.ets(ets(as.vector(med_polish_resi_score[,ip])), h=fh)$mean)
  }
  med_polish_curve_forecast = (med_polish_resi_basis %*% med_polish_resi_score_forecast) + t(fixed_com)[, 1:fh]
  return(list(med_polish_curve_forecast=med_polish_curve_forecast))
}

fore_FANOVA_cdf <- function(Fixed, Residuals, ncomp_method, fh, fmethod, est_method) {
  res = Pref_forecast_curves(Fixed, Residuals, est_method, fh, 1000, fmethod, ncomp_method)$med_polish_curve_forecast
  f_add = c(invlogit(res[, fh]), 1)
  return(c(f_add[1], diff(f_add)) * 10^5)
}

pref_idx <- 1
for(h in 1:3) {
  val_years_indices <- 1:(18 - h)
  val_F = val_M = matrix(NA, n_age, length(val_years_indices))
  
  for(idx in seq_along(val_years_indices)) {
    ij <- val_years_indices[idx]
    val_F[,idx] = fore_FANOVA_cdf(female_prefecture_fixed_means[[pref_idx]][1:(15+ij),], female_prefecture_res_means[[pref_idx]][1:(15+ij),], "EVR", h, "ets", "cov")
    val_M[,idx] = fore_FANOVA_cdf(male_prefecture_fixed_means[[pref_idx]][1:(15+ij),], male_prefecture_res_means[[pref_idx]][1:(15+ij),], "EVR", h, "ets", "cov")
  }
  
  holdout_F = female_prefecture_dx[[pref_idx]][, (16 + h):33]
  holdout_M = male_prefecture_dx[[pref_idx]][, (16 + h):33]
  
  resi_F = holdout_F - val_F
  resi_M = holdout_M - val_M
  sd_F = apply(resi_F, 1, sd)
  sd_M = apply(resi_M, 1, sd)
  
  tune_F = optimise(tune_para_find_function, c(0,20), resi_mat=resi_F, sd_val_input=sd_F, alpha_level=0.8, PI_type="pointwise")$minimum
  tune_M = optimise(tune_para_find_function, c(0,20), resi_mat=resi_M, sd_val_input=sd_M, alpha_level=0.8, PI_type="pointwise")$minimum
  
  fin_F = fore_FANOVA_cdf(female_prefecture_fixed_means[[pref_idx]][1:33,], female_prefecture_res_means[[pref_idx]][1:33,], "EVR", h, "ets", "cov")
  fin_M = fore_FANOVA_cdf(male_prefecture_fixed_means[[pref_idx]][1:33,], male_prefecture_res_means[[pref_idx]][1:33,], "EVR", h, "ets", "cov")
  
  df_p <- rbind(
    data.frame(Age=0:110, Forecast=fin_F, Lower=fin_F-tune_F*sd_F, Upper=fin_F+tune_F*sd_F, Gender="Female"),
    data.frame(Age=0:110, Forecast=fin_M, Lower=fin_M-tune_M*sd_M, Upper=fin_M+tune_M*sd_M, Gender="Male")
  )
  
  plot_and_save(df_p, paste0("FANOVA_h", h, ".pdf"))
}
