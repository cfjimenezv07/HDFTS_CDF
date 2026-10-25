###############################################################################
# Script: Table4_StationaryAnalysis.R
# Description: Stationary analysis of the FTS
#              Produces Table 4 
###############################################################################
source(here::here("load_packages.R"))

# Dynamic Directories
dir.d       <- here::here("data")
dir.results <- here::here("results")
dir.p       <- here::here("results", "Plots")
dir.aux     <- here::here("auxiliary_source")

# Ensure output directories exist
if (!dir.exists(dir.results)) dir.create(dir.results, recursive = TRUE)
if (!dir.exists(dir.p))       dir.create(dir.p, recursive = TRUE)

# 2. PARAMETERS & DATASET INITIALIZATION
# ------------------------------------------------------------------------------
set.seed(12345)
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

# 3. READ SUBNATIONAL MORTALITY DX DATA FROM PRE-PROCESSED RDS
# ------------------------------------------------------------------------------
female_prefecture_dx <- readRDS(file.path(dir.d, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(file.path(dir.d, "male_prefecture_dx.rds"))

# Convert matrices from Year x Age to Age x Year for downstream CDF functions
prefecture_female_dx_t <- lapply(female_prefecture_dx, t)
prefecture_male_dx_t   <- lapply(male_prefecture_dx, t)

# 4. CDF TRANSFORMATION
# ------------------------------------------------------------------------------
source(file.path(dir.aux, "CDF_transformation.R"))

Transformed_female <- list()
Transformed_male   <- list()

for (i in 1:n_prefectures) {
  Transformed_female[[i]] <- t(cdf_transformation(t(prefecture_female_dx_t[[i]]) / 10^5, year))
  Transformed_male[[i]]   <- t(cdf_transformation(t(prefecture_male_dx_t[[i]]) / 10^5, year))
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
center_columns <- function(mat) {
  col_means <- colMeans(mat)
  sweep(mat, MARGIN = 2, STATS = col_means, FUN = "-")
}

run_prefecture_stationarity <- function(data_list, alpha = 0.05) {
  bandwidth_h <- 110^0.5
  pvals <- numeric(length(data_list))
  
  for (i in seq_along(data_list)) {
    mat_centered <- center_columns(data_list[[i]])
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
  Group     = c("Male", "Female"),
  Min       = c(min(male_stat_pvals), min(female_stat_pvals)),
  `1st Qu.` = c(quantile(male_stat_pvals, 0.25), quantile(female_stat_pvals, 0.25)),
  Median    = c(median(male_stat_pvals), median(female_stat_pvals)),
  Mean      = c(mean(male_stat_pvals), mean(female_stat_pvals)),
  `3rd Qu.` = c(quantile(male_stat_pvals, 0.75), quantile(female_stat_pvals, 0.75)),
  Max       = c(max(male_stat_pvals), max(female_stat_pvals))
)

print(stat_summary_table)

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

print(df_fCH)

# --- 7. BUILD LaTeX TABLE 4 DIRECTLY FROM THE P-VALUES -----------------------

summ_p <- function(p) {
  c(Min  = min(p),
    Q1   = unname(quantile(p, 0.25)),
    Med  = median(p),
    Mean = mean(p),
    Q3   = unname(quantile(p, 0.75)),
    Max  = max(p))
}

fmt_fixed <- function(x) sprintf("%.4f", x)
fmt_sci   <- function(x, floor_val = 2.2e-16) {
  sapply(x, function(v) {
    if (v < floor_val) sprintf("$<%.1f\\times 10^{-16}$", floor_val * 1e16)
    else if (v < 1e-3 | v >= 1e3) sprintf("%.2e", v)
    else sprintf("%.4f", v)
  })
}

make_row <- function(label, stats, fmt) {
  paste0(label, " & ", paste(fmt(stats), collapse = " & "), " \\\\")
}

st_F <- summ_p(female_stat_pvals)
st_M <- summ_p(male_stat_pvals)
ch_F <- summ_p(female_fCH_pvals)
ch_M <- summ_p(male_fCH_pvals)

n_ch_F <- sum(female_fCH_pvals < 0.05)
n_ch_M <- sum(male_fCH_pvals   < 0.05)
n_st_F <- sum(female_stat_pvals < 0.05)
n_st_M <- sum(male_stat_pvals   < 0.05)

caption_txt <- paste0(
  "\\caption{\\small Summary of $p$-values from stationarity and conditional ",
  "heteroskedasticity tests across 47 prefectures (with seed set to 12345). ",
  "The stationarity null is rejected at the 5\\% level in ", n_st_F, "/47 female and ",
  n_st_M, "/47 male prefectures. The conditional homoskedasticity null is rejected ",
  "in ", n_ch_F, "/47 female and ", n_ch_M, "/47 male prefectures.}"
)

tab_lines <- c(
  "\\begin{table}[!htb]",
  "\\centering",
  "\\tabcolsep 0.24in",
  caption_txt,
  "\\begin{tabular}{@{}lrrrrrr@{}}",
  "\\toprule",
  "Gender & Min & 1st Qu. & Median & Mean & 3rd Qu. & Max \\\\",
  "\\midrule",
  "\\multicolumn{7}{l}{\\hspace{-.23in}{\\underline{Stationarity test}}} \\\\",
  make_row("Female", st_F, fmt_fixed),
  make_row("Male  ", st_M, fmt_fixed),
  "\\\\",
  "\\multicolumn{7}{l}{\\hspace{-.23in}{\\underline{Conditional heteroskedasticity test}}} \\\\",
  make_row("Female", ch_F, fmt_sci),
  make_row("Male  ", ch_M, fmt_sci),
  "\\bottomrule",
  "\\end{tabular}",
  "\\label{tab:stationarity}",
  "\\end{table}"
)

cat(tab_lines, sep = "\n")
writeLines(tab_lines, file.path(dir.results, "Table4_StationaryAnalysis.tex"))

saveRDS(list(male_stat = male_stat_pvals, female_stat = female_stat_pvals,
             male_fCH  = male_fCH_pvals,  female_fCH  = female_fCH_pvals),
        file.path(dir.results, "Table4_pvalues.rds"))