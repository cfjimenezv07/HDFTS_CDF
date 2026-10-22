###############################################################################
# Script: Table1_PointAnalysis_EVR.R
# Description: Reads saved point forecast error .rds files and updates
#              Table 1 in exact LaTeX format.
###############################################################################
source("load_packages.R")

# Dynamic directory path for Point Forecast results
dir_results <- here("results", "Results_Point")

# 1. READ ALL RDS FILES FOR POINT FORECASTS
# ------------------------------------------------------------------------------
# UFTS
ufts_f_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_F_EVR_ETS_mean.rds"))
ufts_m_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_M_EVR_ETS_mean.rds"))

# MFTS
mfts_f_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MFTS_F_EVR_ETS_mean.rds"))
mfts_m_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MFTS_M_EVR_ETS_mean.rds"))

# MLFTS
mlfts_f_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean.rds"))
mlfts_m_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean.rds"))

# FANOVA
fanova_evr  <- readRDS(file.path(dir_results, "FANOVA_errors_EVR_matrix.rds"))

# HDFPCA
hdfpca_f_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean.rds"))
hdfpca_m_evr <- readRDS(file.path(dir_results, "horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean.rds"))

# 2. EXTRACT KLD & JSD VALUES
# ------------------------------------------------------------------------------
# Formats numbers to exactly 3 decimal places (e.g., 0.004)
fmt <- function(x) sprintf("%.3f", x)

# Horizons 1 to 17 + Mean + Median
rows_idx <- c(1:17, "Mean", "Median")

# KLD Data Matrix (Columns: UFTS_F, MFTS_F, MLFTS_F, FANOVA_F, HDFPCA_F, UFTS_M, MFTS_M, MLFTS_M, FANOVA_M, HDFPCA_M)
kld_data <- cbind(
  ufts_f_evr[rows_idx, "KLD"],
  mfts_f_evr[rows_idx, "KLD"],
  mlfts_f_evr[rows_idx, "KLD"],
  c(fanova_evr[1:17, "KLD_Female"], mean(fanova_evr[1:17, "KLD_Female"]), median(fanova_evr[1:17, "KLD_Female"])),
  hdfpca_f_evr[rows_idx, "KLD"],
  ufts_m_evr[rows_idx, "KLD"],
  mfts_m_evr[rows_idx, "KLD"],
  mlfts_m_evr[rows_idx, "KLD"],
  c(fanova_evr[1:17, "KLD_Male"], mean(fanova_evr[1:17, "KLD_Male"]), median(fanova_evr[1:17, "KLD_Male"])),
  hdfpca_m_evr[rows_idx, "KLD"]
)

# JSD Data Matrix
jsd_data <- cbind(
  ufts_f_evr[rows_idx, "JSD (geo)"],
  mfts_f_evr[rows_idx, "JSD (geo)"],
  mlfts_f_evr[rows_idx, "JSD (geo)"],
  c(fanova_evr[1:17, "JSD_Female"], mean(fanova_evr[1:17, "JSD_Female"]), median(fanova_evr[1:17, "JSD_Female"])),
  hdfpca_f_evr[rows_idx, "JSD (geo)"],
  ufts_m_evr[rows_idx, "JSD (geo)"],
  mfts_m_evr[rows_idx, "JSD (geo)"],
  mlfts_m_evr[rows_idx, "JSD (geo)"],
  c(fanova_evr[1:17, "JSD_Male"], mean(fanova_evr[1:17, "JSD_Male"]), median(fanova_evr[1:17, "JSD_Male"])),
  hdfpca_m_evr[rows_idx, "JSD (geo)"]
)

# Helper function to generate row string
make_row <- function(label, h_val, vals) {
  val_str <- paste(sapply(vals, fmt), collapse = " & ")
  if (is.na(h_val)) {
    return(paste0("     & ", label, " & ", val_str, " \\\\"))
  } else {
    lbl_str <- if (label != "") paste0(label, " ") else "     "
    return(paste0(lbl_str, "&  ", sprintf("%2d", h_val), " & ", val_str, " \\\\"))
  }
}

# 3. PRINT UPDATED LATEX TEMPLATE
# ------------------------------------------------------------------------------
cat('\\begin{center}
\\renewcommand*{\\arraystretch}{0.99}
\\tabcolsep 0.065in
\\begin{footnotesize}
\\begin{longtable}{@{}llrrrrrrrrrrr@{}}
\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the point forecast accuracy of the functional time-series models, measured by the KLD and JSD. The number of components is determined by the EVR criterion.}\\label{tab:1} \\\\
\\toprule
&  & \\multicolumn{5}{c}{Female}   & \\multicolumn{5}{c}{Male}   \\\\
\\cmidrule{3-13}
Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & HDFPCA & UFTS  & MFTS  & MLFTS & FANOVA & HDFPCA \\\\ 
\\midrule
\\endfirsthead
\\toprule
&  & \\multicolumn{5}{c}{Female}   & \\multicolumn{5}{c}{Male}   \\\\
\\cmidrule{3-13}
Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & HDFPCA & UFTS  & MFTS  & MLFTS & FANOVA & HDFPCA \\\\ 
\\midrule
\\endhead
\\midrule
\\multicolumn{13}{r}{{Continued on next page}} \\\\
\\endfoot
\\endlastfoot\n')

# KLD Rows 1-17
for (h in 1:17) {
  lbl <- if (h == 1) "KLD" else ""
  cat(make_row(lbl, h, kld_data[h, ]), "\n")
}

# KLD Mean & Median
cat('\\cmidrule{2-12} \n')
cat(make_row("Mean", NA, kld_data[18, ]), "\n")
cat(make_row("Median", NA, kld_data[19, ]), "\n")
cat('\\midrule\n')

# JSD Rows 1-17
for (h in 1:17) {
  lbl <- if (h == 1) "JSD" else ""
  cat(make_row(lbl, h, jsd_data[h, ]), "\n")
}

# JSD Mean & Median
cat('\\cmidrule{2-12}  \n')
cat(make_row("Mean", NA, jsd_data[18, ]), "\n")
cat(make_row("Median", NA, jsd_data[19, ]), "\n")

cat('\\bottomrule
\\end{longtable}
\\end{footnotesize}
\\end{center}\n')
