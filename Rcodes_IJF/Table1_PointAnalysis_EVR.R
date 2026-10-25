###############################################################################
# Script: Table1_PointAnalysis_EVR.R
# Description: Reads saved point forecast error .rds files and prints
#               Table 1 in the LaTeX format used in the paper, while
#               simultaneously exporting a clean .tex file.
###############################################################################
source(here("load_packages.R"))

# Dynamic directory path for Point Forecast results
dir_results <- here("results", "Results_Point")
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

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

# KLD Data Matrix
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

# Helper function to generate row string cleanly without hidden spaces
make_row <- function(label, h_val, vals) {
  val_str <- paste(sapply(vals, fmt), collapse = " & ")
  lbl_str <- sprintf("%-6s", label)
  
  if (is.na(h_val)) {
    return(paste0(lbl_str, " &        & ", val_str, " \\\\"))
  } else {
    return(paste0(lbl_str, " & ", sprintf("%4d", h_val), " & ", val_str, " \\\\"))
  }
}

# 3. GENERATE LATEX TABLE CODE
# ------------------------------------------------------------------------------
generate_table_code <- function() {
  lines <- c()
  
  lines <- c(lines, '\\begin{center}')
  lines <- c(lines, '\\renewcommand*{\\arraystretch}{0.99}')
  lines <- c(lines, '\\tabcolsep 0.065in')
  lines <- c(lines, '\\begin{footnotesize}')
  lines <- c(lines, '\\begin{longtable}{@{}llrrrrrrrrrrr@{}}')
  lines <- c(lines, '\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the point forecast accuracy of the functional time-series models, measured by the KLD and JSD. The number of components is determined by the EVR criterion.}\\label{tab:1} \\\\')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '&   & \\multicolumn{5}{c}{Female}   & \\multicolumn{5}{c}{Male}   \\\\')
  lines <- c(lines, '\\cmidrule{3-13}')
  lines <- c(lines, 'Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & HDFPCA & UFTS  & MFTS  & MLFTS & FANOVA & HDFPCA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endfirsthead')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '&   & \\multicolumn{5}{c}{Female}   & \\multicolumn{5}{c}{Male}   \\\\')
  lines <- c(lines, '\\cmidrule{3-13}')
  lines <- c(lines, 'Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & HDFPCA & UFTS  & MFTS  & MLFTS & FANOVA & HDFPCA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endhead')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\multicolumn{13}{r}{{Continued on next page}} \\\\')
  lines <- c(lines, '\\endfoot')
  lines <- c(lines, '\\endlastfoot')
  
  # KLD Rows 1-17
  for (h in 1:17) {
    lbl <- if (h == 1) "KLD" else ""
    lines <- c(lines, make_row(lbl, h, kld_data[h, ]))
  }
  
  # KLD Mean & Median
  lines <- c(lines, '\\cmidrule{2-13}')
  lines <- c(lines, make_row("Mean", NA, kld_data[18, ]))
  lines <- c(lines, make_row("Median", NA, kld_data[19, ]))
  lines <- c(lines, '\\midrule')
  
  # JSD Rows 1-17
  for (h in 1:17) {
    lbl <- if (h == 1) "JSD" else ""
    lines <- c(lines, make_row(lbl, h, jsd_data[h, ]))
  }
  
  # JSD Mean & Median
  lines <- c(lines, '\\cmidrule{2-13}')
  lines <- c(lines, make_row("Mean", NA, jsd_data[18, ]))
  lines <- c(lines, make_row("Median", NA, jsd_data[19, ]))
  
  lines <- c(lines, '\\bottomrule')
  lines <- c(lines, '\\end{longtable}')
  lines <- c(lines, '\\end{footnotesize}')
  lines <- c(lines, '\\end{center}')
  
  return(lines)
}

# 4. PRINT TO CONSOLE & EXPORT FILE
# ------------------------------------------------------------------------------
latex_table_lines <- generate_table_code()

# Print directly to R console
cat(paste(latex_table_lines, collapse = "\n"), "\n")

# Export directly to results/Table1.tex
output_tex_path <- file.path(dir_results, "Table1.tex")
writeLines(latex_table_lines, con = output_tex_path)
message("\nSuccess: Table 1 LaTeX exported directly to: ", output_tex_path)
