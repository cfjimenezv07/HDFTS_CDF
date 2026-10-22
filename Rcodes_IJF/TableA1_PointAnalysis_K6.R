###############################################################################
# Script: Table1_PointAnalysis_K6.R
# Description: Reads saved point forecast error .rds files (K = 6) and updates
#              Table A.1 in exact LaTeX format.
###############################################################################
source("load_packages.R")

# Dynamic directory path for Point Forecast results
dir_results <- here("results", "Results_Point")

# 1. READ ALL RDS FILES FOR POINT FORECASTS (K = 6 CRITERION)
# ------------------------------------------------------------------------------
# UFTS
ufts_f_k6 <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_F_K6_ETS_mean.rds"))
ufts_m_k6 <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_M_K6_ETS_mean.rds"))

# MFTS
mfts_f_k6 <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MFTS_F_K6_ETS_mean.rds"))
mfts_m_k6 <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MFTS_M_K6_ETS_mean.rds"))

# MLFTS
mlfts_f_k6 <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean.rds"))
mlfts_m_k6 <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean.rds"))

# FANOVA
fanova_k6  <- readRDS(file.path(dir_results, "FANOVA_errors_K_matrix.rds"))

# 2. EXTRACT KLD & JSD VALUES
# ------------------------------------------------------------------------------
# Formats numbers to exactly 3 decimal places (e.g., 0.004)
fmt <- function(x) sprintf("%.3f", x)

# Horizons 1 to 17 + Mean + Median
rows_idx <- c(1:17, "Mean", "Median")

# KLD Data Matrix (8 columns: UFTS_F, MFTS_F, MLFTS_F, FANOVA_F, UFTS_M, MFTS_M, MLFTS_M, FANOVA_M)
kld_data <- cbind(
  ufts_f_k6[rows_idx, "KLD"],
  mfts_f_k6[rows_idx, "KLD"],
  mlfts_f_k6[rows_idx, "KLD"],
  c(fanova_k6[1:17, 1], mean(fanova_k6[1:17, 1]), median(fanova_k6[1:17, 1])),
  ufts_m_k6[rows_idx, "KLD"],
  mfts_m_k6[rows_idx, "KLD"],
  mlfts_m_k6[rows_idx, "KLD"],
  c(fanova_k6[1:17, 3], mean(fanova_k6[1:17, 3]), median(fanova_k6[1:17, 3]))
)

# JSD Data Matrix (8 columns)
jsd_data <- cbind(
  ufts_f_k6[rows_idx, "JSD (geo)"],
  mfts_f_k6[rows_idx, "JSD (geo)"],
  mlfts_f_k6[rows_idx, "JSD (geo)"],
  c(fanova_k6[1:17, 2], mean(fanova_k6[1:17, 2]), median(fanova_k6[1:17, 2])),
  ufts_m_k6[rows_idx, "JSD (geo)"],
  mfts_m_k6[rows_idx, "JSD (geo)"],
  mlfts_m_k6[rows_idx, "JSD (geo)"],
  c(fanova_k6[1:17, 4], mean(fanova_k6[1:17, 4]), median(fanova_k6[1:17, 4]))
)

# Helper function to generate row string
make_row <- function(label, h_val, vals) {
  val_str <- paste(sapply(vals, fmt), collapse = " & ")
  if (is.na(h_val)) {
    return(paste0("    & ", label, " & ", val_str, " \\\\"))
  } else {
    lbl_str <- if (label != "") paste0(label, " ") else "    "
    return(paste0(lbl_str, "& ", sprintf("%2d", h_val), " & ", val_str, " \\\\"))
  }
}

# 3. PRINT UPDATED LATEX TEMPLATE FOR TABLE A.1
# ------------------------------------------------------------------------------
cat('\\begin{center}
\\renewcommand*{\\arraystretch}{0.99}
\\tabcolsep 0.1in
\\begin{longtable}{@{}llrrrrrrrr@{}}
\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the point forecast accuracy of the functional time-series models, measured by the KLD and JSD. The number of components is determined by setting $K=6$.}\\label{tab:appendix_1} \\\\
\\toprule
&  & \\multicolumn{4}{c}{Female}   & \\multicolumn{4}{c}{Male}   \\\\
\\cmidrule{3-10}
Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & UFTS  & MFTS  & MLFTS & FANOVA \\\\
\\midrule
\\endfirsthead
\\toprule
&  & \\multicolumn{4}{c}{Female}   & \\multicolumn{4}{c}{Male}   \\\\
\\cmidrule{3-10}
Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & UFTS  & MFTS  & MLFTS & FANOVA \\\\
\\midrule
\\endhead
\\midrule
\\multicolumn{10}{r}{{Continued on next page}} \\\\
\\endfoot
\\endlastfoot\n')

# KLD Rows 1-17
for (h in 1:17) {
  lbl <- if (h == 1) "KLD" else ""
  cat(make_row(lbl, h, kld_data[h, ]), "\n")
}

# KLD Mean & Median
cat('\\cmidrule{1-10}\n')
cat(make_row("Mean", NA, kld_data[18, ]), "\n")
cat(make_row("Median", NA, kld_data[19, ]), "\n")
cat('\\midrule\n')

# JSD Rows 1-17
for (h in 1:17) {
  lbl <- if (h == 1) "JSD" else ""
  cat(make_row(lbl, h, jsd_data[h, ]), "\n")
}

# JSD Mean & Median
cat('\\cmidrule{1-10}\n')
cat(make_row("Mean", NA, jsd_data[18, ]), "\n")
cat(make_row("Median", NA, jsd_data[19, ]), "\n")

cat('\\bottomrule
\\end{longtable}
\\end{center}\n')

