###############################################################################
# Script: TableB1_IntervalAnalysis_K6.R
# Description: Reads saved intervsal forecast error (K6) .rds files and prints
#               Table B1 in the LaTeX format used in the paper.
###############################################################################
source(here("load_packages.R"))

# Dynamic directory path for Interval Forecast (SD) results
dir.results <- paste0(here("results", "Results_Interval", "sd"), "/")

# 1. READ ALL RDS FILES FOR INTERVAL FORECASTS (SD, K = 6 CRITERION)
# ------------------------------------------------------------------------------
# UFTS
ufts_f_int_k6 <- readRDS(file.path(dir.results, "int_fore_subnational_err_F_K6_ETS_mean.rds"))
ufts_m_int_k6 <- readRDS(file.path(dir.results, "int_fore_subnational_err_M_K6_ETS_mean.rds"))

# MFTS
mfts_f_int_k6 <- readRDS(file.path(dir.results, "MFTS_int_fore_subnational_err_F_K6_ETS_mean.rds"))
mfts_m_int_k6 <- readRDS(file.path(dir.results, "MFTS_int_fore_subnational_err_M_K6_ETS_mean.rds"))

# MLFTS
mlfts_f_int_k6 <- readRDS(file.path(dir.results, "MLFTS_int_fore_subnational_err_F_K6_ETS_mean.rds"))
mlfts_m_int_k6 <- readRDS(file.path(dir.results, "MLFTS_int_fore_subnational_err_M_K6_ETS_mean.rds"))

# FANOVA
fanova_f_int_k6 <- readRDS(file.path(dir.results, "FANOVA_int_fore_subnational_err_F_ETS_K6_mean.rds"))
fanova_m_int_k6 <- readRDS(file.path(dir.results, "FANOVA_int_fore_subnational_err_M_ETS_K6_mean.rds"))

# 2. EXTRACT CPD & SCORE VALUES
# ------------------------------------------------------------------------------
# Formatting helper functions
fmt_cpd   <- function(x) sprintf("%.3f", x)
fmt_score <- function(x) sprintf("%.0f", round(x))

# 16 Horizons
h_16 <- 1:16

# Extract CPD (16 horizons + Mean + Median)
cpd_data <- cbind(
  c(ufts_f_int_k6[h_16, "CPD"],   mean(ufts_f_int_k6[h_16, "CPD"]),   median(ufts_f_int_k6[h_16, "CPD"])),
  c(mfts_f_int_k6[h_16, "CPD"],   mean(mfts_f_int_k6[h_16, "CPD"]),   median(mfts_f_int_k6[h_16, "CPD"])),
  c(mlfts_f_int_k6[h_16, "CPD"],  mean(mlfts_f_int_k6[h_16, "CPD"]),  median(mlfts_f_int_k6[h_16, "CPD"])),
  c(fanova_f_int_k6[h_16, "CPD"], mean(fanova_f_int_k6[h_16, "CPD"]), median(fanova_f_int_k6[h_16, "CPD"])),
  c(ufts_m_int_k6[h_16, "CPD"],   mean(ufts_m_int_k6[h_16, "CPD"]),   median(ufts_m_int_k6[h_16, "CPD"])),
  c(mfts_m_int_k6[h_16, "CPD"],   mean(mfts_m_int_k6[h_16, "CPD"]),   median(mfts_m_int_k6[h_16, "CPD"])),
  c(mlfts_m_int_k6[h_16, "CPD"],  mean(mlfts_m_int_k6[h_16, "CPD"]),  median(mlfts_m_int_k6[h_16, "CPD"])),
  c(fanova_m_int_k6[h_16, "CPD"], mean(fanova_m_int_k6[h_16, "CPD"]), median(fanova_m_int_k6[h_16, "CPD"]))
)

# Extract Interval Score (16 horizons + Mean + Median)
score_data <- cbind(
  c(ufts_f_int_k6[h_16, "score"],   mean(ufts_f_int_k6[h_16, "score"]),   median(ufts_f_int_k6[h_16, "score"])),
  c(mfts_f_int_k6[h_16, "score"],   mean(mfts_f_int_k6[h_16, "score"]),   median(mfts_f_int_k6[h_16, "score"])),
  c(mlfts_f_int_k6[h_16, "score"],  mean(mlfts_f_int_k6[h_16, "score"]),  median(mlfts_f_int_k6[h_16, "score"])),
  c(fanova_f_int_k6[h_16, "score"], mean(fanova_f_int_k6[h_16, "score"]), median(fanova_f_int_k6[h_16, "score"])),
  c(ufts_m_int_k6[h_16, "score"],   mean(ufts_m_int_k6[h_16, "score"]),   median(ufts_m_int_k6[h_16, "score"])),
  c(mfts_m_int_k6[h_16, "score"],   mean(mfts_m_int_k6[h_16, "score"]),   median(mfts_m_int_k6[h_16, "score"])),
  c(mlfts_m_int_k6[h_16, "score"],  mean(mlfts_m_int_k6[h_16, "score"]),  median(mlfts_m_int_k6[h_16, "score"])),
  c(fanova_m_int_k6[h_16, "score"], mean(fanova_m_int_k6[h_16, "score"]), median(fanova_m_int_k6[h_16, "score"]))
)

# Helper function to generate row string
make_row <- function(label, h_val, vals, is_score = FALSE) {
  fmt_fn <- if (is_score) fmt_score else fmt_cpd
  val_str <- paste(sapply(vals, fmt_fn), collapse = " & ")
  if (is.na(h_val)) {
    return(paste0("    & ", label, " & ", val_str, " \\\\"))
  } else {
    lbl_str <- if (label != "") label else ""
    return(paste0(lbl_str, " & ", h_val, " & ", val_str, " \\\\"))
  }
}

# 3. PRINT UPDATED LATEX TEMPLATE FOR TABLE B1
# ------------------------------------------------------------------------------
cat('\\begin{center}
\\renewcommand*{\\arraystretch}{0.96}
\\tabcolsep 0.1in
\\begin{longtable}{@{}llrrrrrrrr@{}}
\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the interval forecast accuracy, as measured by the CPD and interval score, among the functional time series models at the level of significance $\alpha=0.2$. The number of components is determined by setting $K=6$. For either females or males, we highlight the method with the smallest overall CPD and interval scores.}\\label{tab:appendix_2_K6}\\\\
\\toprule  
  & & \\multicolumn{4}{c}{Female}   & \\multicolumn{4}{c}{Male}   \\\\
  Metric & $h$ & UFTS & MFTS & MLFTS & FANOVA & UFTS & MFTS & MLFTS & FANOVA \\\\
\\midrule
\\endfirsthead
\\toprule  
  & & \\multicolumn{4}{c}{Female}   & \\multicolumn{4}{c}{Male}   \\\\
  Metric & $h$ & UFTS & MFTS & MLFTS & FANOVA & UFTS & MFTS & MLFTS & FANOVA \\\\
\\midrule
\\endhead
\\midrule
\\multicolumn{10}{r}{{Continued on next page}} \\\\
\\endfoot
\\endlastfoot\n')

# CPD Rows 1-16
for (h in 1:16) {
  lbl <- if (h == 1) "CPD$_{0.2}$" else ""
  cat(make_row(lbl, h, cpd_data[h, ], is_score = FALSE), "\n")
}

# CPD Mean & Median
cat('\\cmidrule{1-10}\n')
cat(make_row("Mean", NA, cpd_data[17, ], is_score = FALSE), "\n")
cat(make_row("Median", NA, cpd_data[18, ], is_score = FALSE), "\n")
cat('\\midrule\n')

# Interval Score Rows 1-16
for (h in 1:16) {
  lbl <- if (h == 1) "$\\overline{S}_{0.2}$" else ""
  cat(make_row(lbl, h, score_data[h, ], is_score = TRUE), "\n")
}

# Interval Score Mean & Median
cat('\\cmidrule{1-10}\n')
cat(make_row("Mean", NA, score_data[17, ], is_score = TRUE), "\n")
cat(make_row("Median", NA, score_data[18, ], is_score = TRUE), "\n")

cat('\\bottomrule
\\end{longtable}
\\end{center}\n')

