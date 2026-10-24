###############################################################################
# Script: TableA1_PointAnalysis_K6.R
# Description: Reads saved point forecast error .rds files and prints
#               Table A1 in the LaTeX format used in the paper, while
#               simultaneously exporting a clean .tex file.
###############################################################################
source(here("load_packages.R"))

# Dynamic directory path for Point Forecast results
dir_results <- here("results", "Results_Point")
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

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

# 3. GENERATE LATEX TABLE CODE FOR TABLE A.1
# ------------------------------------------------------------------------------
generate_table_code <- function() {
  lines <- c()
  
  lines <- c(lines, '\\begin{center}')
  lines <- c(lines, '\\renewcommand*{\\arraystretch}{0.99}')
  lines <- c(lines, '\\tabcolsep 0.1in')
  lines <- c(lines, '\\begin{longtable}{@{}llrrrrrrrr@{}}')
  lines <- c(lines, '\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the point forecast accuracy of the functional time-series models, measured by the KLD and JSD. The number of components is determined by setting $K=6$.}\\label{tab:appendix_1} \\\\')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '&   & \\multicolumn{4}{c}{Female}   & \\multicolumn{4}{c}{Male}   \\\\')
  lines <- c(lines, '\\cmidrule{3-10}')
  lines <- c(lines, 'Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & UFTS  & MFTS  & MLFTS & FANOVA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endfirsthead')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '&   & \\multicolumn{4}{c}{Female}   & \\multicolumn{4}{c}{Male}   \\\\')
  lines <- c(lines, '\\cmidrule{3-10}')
  lines <- c(lines, 'Metric & $h$ & UFTS & MFTS  & MLFTS & FANOVA & UFTS  & MFTS  & MLFTS & FANOVA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endhead')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\multicolumn{10}{r}{{Continued on next page}} \\\\')
  lines <- c(lines, '\\endfoot')
  lines <- c(lines, '\\endlastfoot')
  
  # KLD Rows 1-17
  for (h in 1:17) {
    lbl <- if (h == 1) "KLD" else ""
    lines <- c(lines, make_row(lbl, h, kld_data[h, ]))
  }
  
  # KLD Mean & Median (FIXED RANGE TO 2-10)
  lines <- c(lines, '\\cmidrule{2-10}')
  lines <- c(lines, make_row("Mean", NA, kld_data[18, ]))
  lines <- c(lines, make_row("Median", NA, kld_data[19, ]))
  lines <- c(lines, '\\midrule')
  
  # JSD Rows 1-17
  for (h in 1:17) {
    lbl <- if (h == 1) "JSD" else ""
    lines <- c(lines, make_row(lbl, h, jsd_data[h, ]))
  }
  
  # JSD Mean & Median (FIXED RANGE TO 2-10)
  lines <- c(lines, '\\cmidrule{2-10}')
  lines <- c(lines, make_row("Mean", NA, jsd_data[18, ]))
  lines <- c(lines, make_row("Median", NA, jsd_data[19, ]))
  
  lines <- c(lines, '\\bottomrule')
  lines <- c(lines, '\\end{longtable}')
  lines <- c(lines, '\\end{center}')
  
  return(lines)
}

# 4. PRINT TO CONSOLE & EXPORT FILE
# ------------------------------------------------------------------------------
latex_table_lines <- generate_table_code()

# Print directly to R console
cat(paste(latex_table_lines, collapse = "\n"), "\n")

# Export directly to results/Results_Point/TableA1.tex
output_tex_path <- file.path(dir_results, "TableA1.tex")
writeLines(latex_table_lines, con = output_tex_path)
message("\nSuccess: Table A1 LaTeX exported directly to: ", output_tex_path)