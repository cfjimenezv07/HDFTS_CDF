###############################################################################
# Script: Table3_IntervalAnalysis_EVR.R
# Description: Reads saved interval forecast error (sd) .rds files and prints
#               Table 3 in the LaTeX format used in the paper, while
#               simultaneously exporting a clean .tex file.
###############################################################################
source(here("load_packages.R"))

# Dynamic directory path for Interval Forecast (SD) results
dir_results <- here("results", "Results_Interval", "sd")
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

# 1. READ ALL RDS FILES FOR INTERVAL FORECASTS (SD, EVR CRITERION)
# ------------------------------------------------------------------------------
# UFTS
ufts_f_int_evr <- readRDS(file.path(dir_results, "int_fore_subnational_err_F_EVR_ETS_mean.rds"))
ufts_m_int_evr <- readRDS(file.path(dir_results, "int_fore_subnational_err_M_EVR_ETS_mean.rds"))

# MFTS
mfts_f_int_evr <- readRDS(file.path(dir_results, "MFTS_int_fore_subnational_err_F_EVR_ETS_mean.rds"))
mfts_m_int_evr <- readRDS(file.path(dir_results, "MFTS_int_fore_subnational_err_M_EVR_ETS_mean.rds"))

# MLFTS
mlfts_f_int_evr <- readRDS(file.path(dir_results, "MLFTS_int_fore_subnational_err_F_EVR_ETS_mean.rds"))
mlfts_m_int_evr <- readRDS(file.path(dir_results, "MLFTS_int_fore_subnational_err_M_EVR_ETS_mean.rds"))

# FANOVA
fanova_f_int_evr <- readRDS(file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_ETS_mean.rds"))
fanova_m_int_evr <- readRDS(file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_ETS_mean.rds"))

# HDFPCA
hdfpca_f_int_evr <- readRDS(file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_CDF.rds"))
hdfpca_m_int_evr <- readRDS(file.path(dir_results, "horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_CDF.rds"))

# 2. EXTRACT CPD & SCORE VALUES
# ------------------------------------------------------------------------------
# Formatting helper functions
fmt_cpd   <- function(x) sprintf("%.3f", x)
fmt_score <- function(x) sprintf("%.0f", round(x))

# 16 Horizons
h_16 <- 1:16

# Extract CPD (16 horizons + Mean + Median)
cpd_data <- cbind(
  c(ufts_f_int_evr[h_16, "CPD"],   mean(ufts_f_int_evr[h_16, "CPD"]),   median(ufts_f_int_evr[h_16, "CPD"])),
  c(mfts_f_int_evr[h_16, "CPD"],   mean(mfts_f_int_evr[h_16, "CPD"]),   median(mfts_f_int_evr[h_16, "CPD"])),
  c(mlfts_f_int_evr[h_16, "CPD"],  mean(mlfts_f_int_evr[h_16, "CPD"]),  median(mlfts_f_int_evr[h_16, "CPD"])),
  c(fanova_f_int_evr[h_16, "CPD"], mean(fanova_f_int_evr[h_16, "CPD"]), median(fanova_f_int_evr[h_16, "CPD"])),
  c(hdfpca_f_int_evr[h_16, "CPD"], mean(hdfpca_f_int_evr[h_16, "CPD"]), median(hdfpca_f_int_evr[h_16, "CPD"])),
  c(ufts_m_int_evr[h_16, "CPD"],   mean(ufts_m_int_evr[h_16, "CPD"]),   median(ufts_m_int_evr[h_16, "CPD"])),
  c(mfts_m_int_evr[h_16, "CPD"],   mean(mfts_m_int_evr[h_16, "CPD"]),   median(mfts_m_int_evr[h_16, "CPD"])),
  c(mlfts_m_int_evr[h_16, "CPD"],  mean(mlfts_m_int_evr[h_16, "CPD"]),  median(mlfts_m_int_evr[h_16, "CPD"])),
  c(fanova_m_int_evr[h_16, "CPD"], mean(fanova_m_int_evr[h_16, "CPD"]), median(fanova_m_int_evr[h_16, "CPD"])),
  c(hdfpca_m_int_evr[h_16, "CPD"], mean(hdfpca_m_int_evr[h_16, "CPD"]), median(hdfpca_m_int_evr[h_16, "CPD"]))
)

# Extract Interval Score (16 horizons + Mean + Median)
score_data <- cbind(
  c(ufts_f_int_evr[h_16, "score"],   mean(ufts_f_int_evr[h_16, "score"]),   median(ufts_f_int_evr[h_16, "score"])),
  c(mfts_f_int_evr[h_16, "score"],   mean(mfts_f_int_evr[h_16, "score"]),   median(mfts_f_int_evr[h_16, "score"])),
  c(mlfts_f_int_evr[h_16, "score"],  mean(mlfts_f_int_evr[h_16, "score"]),  median(mlfts_f_int_evr[h_16, "score"])),
  c(fanova_f_int_evr[h_16, "score"], mean(fanova_f_int_evr[h_16, "score"]), median(fanova_f_int_evr[h_16, "score"])),
  c(hdfpca_f_int_evr[h_16, "score"], mean(hdfpca_f_int_evr[h_16, "score"]), median(hdfpca_f_int_evr[h_16, "score"])),
  c(ufts_m_int_evr[h_16, "score"],   mean(ufts_m_int_evr[h_16, "score"]),   median(ufts_m_int_evr[h_16, "score"])),
  c(mfts_m_int_evr[h_16, "score"],   mean(mfts_m_int_evr[h_16, "score"]),   median(mfts_m_int_evr[h_16, "score"])),
  c(mlfts_m_int_evr[h_16, "score"],  mean(mlfts_m_int_evr[h_16, "score"]),  median(mlfts_m_int_evr[h_16, "score"])),
  c(fanova_m_int_evr[h_16, "score"], mean(fanova_m_int_evr[h_16, "score"]), median(fanova_m_int_evr[h_16, "score"])),
  c(hdfpca_m_int_evr[h_16, "score"], mean(hdfpca_m_int_evr[h_16, "score"]), median(hdfpca_m_int_evr[h_16, "score"]))
)

# Helper function to generate row string
make_row <- function(label, h_val, vals, is_score = FALSE) {
  fmt_fn  <- if (is_score) fmt_score else fmt_cpd
  val_str <- paste(sapply(vals, fmt_fn), collapse = " & ")
  lbl_str <- sprintf("%-25s", label)
  
  if (is.na(h_val)) {
    return(paste0(lbl_str, " &        & ", val_str, " \\\\"))
  } else {
    return(paste0(lbl_str, " & ", sprintf("%4d", h_val), " & ", val_str, " \\\\"))
  }
}

# 3. GENERATE LATEX TABLE CODE FOR TABLE 3
# ------------------------------------------------------------------------------
generate_table_code <- function() {
  lines <- c()
  
  lines <- c(lines, '\\begin{center}')
  lines <- c(lines, '\\renewcommand*{\\arraystretch}{1}')
  lines <- c(lines, '\\tabcolsep 0.065in')
  lines <- c(lines, '\\begin{footnotesize}')
  lines <- c(lines, '\\begin{longtable}{@{}llrrrrrrrrrrr@{}}')
  lines <- c(lines, '\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the interval forecast accuracy, as measured by the CPD and interval score, for the nine functional time series models at the significance level $\alpha=0.2$. The number of components is determined by the EVR criterion. We highlight the method with the smallest overall $\overline{\text{CPD}}$ and $\overline{S}_{0.2}$. }\\label{tab:2} \\\\')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '&   & \\multicolumn{5}{c}{Female}   & \\multicolumn{5}{c}{Male}   \\\\')
  lines <- c(lines, '\\cmidrule{3-13}')
  lines <- c(lines, 'Metric & $h$ & UFTS & MFTS  & MLFTS  & FANOVA & HDFPCA & UFTS  & MFTS  & MLFTS  & FANOVA & HDFPCA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endfirsthead')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '&   & \\multicolumn{5}{c}{Female}   & \\multicolumn{5}{c}{Male}   \\\\')
  lines <- c(lines, '\\cmidrule{3-13}')
  lines <- c(lines, 'Metric & $h$ & UFTS & MFTS  & MLFTS  & FANOVA & HDFPCA & UFTS  & MFTS  & MLFTS  & FANOVA & HDFPCA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endhead')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\multicolumn{13}{r}{{Continued on next page}} \\\\')
  lines <- c(lines, '\\endfoot')
  lines <- c(lines, '\\endlastfoot')
  
  # CPD Rows 1-16
  for (h in 1:16) {
    lbl <- if (h == 1) "$\\overline{\\text{CPD}}_{0.2}$" else ""
    lines <- c(lines, make_row(lbl, h, cpd_data[h, ], is_score = FALSE))
  }
  
  # CPD Mean & Median (FIXED RANGE TO 2-13)
  lines <- c(lines, '\\cmidrule{2-13}')
  lines <- c(lines, make_row("Mean", NA, cpd_data[17, ], is_score = FALSE))
  lines <- c(lines, make_row("Median", NA, cpd_data[18, ], is_score = FALSE))
  lines <- c(lines, '\\midrule')
  
  # Interval Score Rows 1-16
  for (h in 1:16) {
    lbl <- if (h == 1) "$\\overline{S}_{0.2}$" else ""
    lines <- c(lines, make_row(lbl, h, score_data[h, ], is_score = TRUE))
  }
  
  # Interval Score Mean & Median (FIXED RANGE TO 2-13)
  lines <- c(lines, '\\cmidrule{2-13}')
  lines <- c(lines, make_row("Mean", NA, score_data[17, ], is_score = TRUE))
  lines <- c(lines, make_row("Median", NA, score_data[18, ], is_score = TRUE))
  
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

# Export directly to results/Results_Interval/sd/Table3.tex
output_tex_path <- file.path(dir_results, "Table3.tex")
writeLines(latex_table_lines, con = output_tex_path)
message("\nSuccess: Table 3 LaTeX exported directly to: ", output_tex_path)