###############################################################################
# Script: TableC1_ConformalAnalysis.R
# Description: Reads saved conformal interval forecast error .rds files and 
#              updates Table C.1 in exact LaTeX format.
###############################################################################
source("load_packages.R")

# Dynamic directory path for Conformal Interval Forecast results
dir.results <- paste0(here("results", "Results_Interval", "conformal"), "/")

# 1. READ ALL RDS FILES FOR CONFORMAL INTERVAL FORECASTS
# ------------------------------------------------------------------------------
# --- EVR CRITERION ---
# UFTS
ufts_f_evr <- readRDS(file.path(dir.results, "int_fore_subnational_err_F_EVR_conformal_mean.rds"))
ufts_m_evr <- readRDS(file.path(dir.results, "int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# MFTS
mfts_f_evr <- readRDS(file.path(dir.results, "MFTS_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
mfts_m_evr <- readRDS(file.path(dir.results, "MFTS_int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# MLFTS
mlfts_f_evr <- readRDS(file.path(dir.results, "MLFTS_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
mlfts_m_evr <- readRDS(file.path(dir.results, "MLFTS_int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# FANOVA
fanova_f_evr <- readRDS(file.path(dir.results, "FANOVA_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
fanova_m_evr <- readRDS(file.path(dir.results, "FANOVA_int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# HDFPCA
hdfpca_f_evr <- readRDS(file.path(dir.results, "hdfpca_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
hdfpca_m_evr <- readRDS(file.path(dir.results, "hdfpca_int_fore_subnational_err_M_EVR_conformal_mean.rds"))


# --- K = 6 CRITERION ---
# UFTS
ufts_f_k6 <- readRDS(file.path(dir.results, "int_fore_subnational_err_F_K6_conformal_mean.rds"))
ufts_m_k6 <- readRDS(file.path(dir.results, "int_fore_subnational_err_M_K6_conformal_mean.rds"))

# MFTS
mfts_f_k6 <- readRDS(file.path(dir.results, "MFTS_int_fore_subnational_err_F_K6_conformal_mean.rds"))
mfts_m_k6 <- readRDS(file.path(dir.results, "MFTS_int_fore_subnational_err_M_K6_conformal_mean.rds"))

# MLFTS
mlfts_f_k6 <- readRDS(file.path(dir.results, "MLFTS_int_fore_subnational_err_F_K6_conformal_mean.rds"))
mlfts_m_k6 <- readRDS(file.path(dir.results, "MLFTS_int_fore_subnational_err_M_K6_conformal_mean.rds"))

# FANOVA
fanova_f_k6 <- readRDS(file.path(dir.results, "FANOVA_int_fore_subnational_err_F_K6_conformal_mean.rds"))
fanova_m_k6 <- readRDS(file.path(dir.results, "FANOVA_int_fore_subnational_err_M_K6_conformal_mean.rds"))

# 2. EXTRACT CPD & SCORE DATA MATRICES
# ------------------------------------------------------------------------------
fmt_cpd   <- function(x) sprintf("%.3f", x)
fmt_score <- function(x) sprintf("%.0f", round(x))

h_16 <- 1:16

# --- CPD FEMALE ---
cpd_f_matrix <- cbind(
  c(ufts_f_evr[h_16, "CPD"],   mean(ufts_f_evr[h_16, "CPD"]),   median(ufts_f_evr[h_16, "CPD"])),
  c(mfts_f_evr[h_16, "CPD"],   mean(mfts_f_evr[h_16, "CPD"]),   median(mfts_f_evr[h_16, "CPD"])),
  c(mlfts_f_evr[h_16, "CPD"],  mean(mlfts_f_evr[h_16, "CPD"]),  median(mlfts_f_evr[h_16, "CPD"])),
  c(fanova_f_evr[h_16, "CPD"], mean(fanova_f_evr[h_16, "CPD"]), median(fanova_f_evr[h_16, "CPD"])),
  c(hdfpca_f_evr[h_16, "CPD"], mean(hdfpca_f_evr[h_16, "CPD"]), median(hdfpca_f_evr[h_16, "CPD"])),
  c(ufts_f_k6[h_16, "CPD"],    mean(ufts_f_k6[h_16, "CPD"]),    median(ufts_f_k6[h_16, "CPD"])),
  c(mfts_f_k6[h_16, "CPD"],    mean(mfts_f_k6[h_16, "CPD"]),    median(mfts_f_k6[h_16, "CPD"])),
  c(mlfts_f_k6[h_16, "CPD"],   mean(mlfts_f_k6[h_16, "CPD"]),   median(mlfts_f_k6[h_16, "CPD"])),
  c(fanova_f_k6[h_16, "CPD"],  mean(fanova_f_k6[h_16, "CPD"]),  median(fanova_f_k6[h_16, "CPD"]))
)

# --- CPD MALE ---
cpd_m_matrix <- cbind(
  c(ufts_m_evr[h_16, "CPD"],   mean(ufts_m_evr[h_16, "CPD"]),   median(ufts_m_evr[h_16, "CPD"])),
  c(mfts_m_evr[h_16, "CPD"],   mean(mfts_m_evr[h_16, "CPD"]),   median(mfts_m_evr[h_16, "CPD"])),
  c(mlfts_m_evr[h_16, "CPD"],  mean(mlfts_m_evr[h_16, "CPD"]),  median(mlfts_m_evr[h_16, "CPD"])),
  c(fanova_m_evr[h_16, "CPD"], mean(fanova_m_evr[h_16, "CPD"]), median(fanova_m_evr[h_16, "CPD"])),
  c(hdfpca_m_evr[h_16, "CPD"], mean(hdfpca_m_evr[h_16, "CPD"]), median(hdfpca_m_evr[h_16, "CPD"])),
  c(ufts_m_k6[h_16, "CPD"],    mean(ufts_m_k6[h_16, "CPD"]),    median(ufts_m_k6[h_16, "CPD"])),
  c(mfts_m_k6[h_16, "CPD"],    mean(mfts_m_k6[h_16, "CPD"]),    median(mfts_m_k6[h_16, "CPD"])),
  c(mlfts_m_k6[h_16, "CPD"],   mean(mlfts_m_k6[h_16, "CPD"]),   median(mlfts_m_k6[h_16, "CPD"])),
  c(fanova_m_k6[h_16, "CPD"],  mean(fanova_m_k6[h_16, "CPD"]),  median(fanova_m_k6[h_16, "CPD"]))
)

# --- INTERVAL SCORE FEMALE ---
score_f_matrix <- cbind(
  c(ufts_f_evr[h_16, "score"],   mean(ufts_f_evr[h_16, "score"]),   median(ufts_f_evr[h_16, "score"])),
  c(mfts_f_evr[h_16, "score"],   mean(mfts_f_evr[h_16, "score"]),   median(mfts_f_evr[h_16, "score"])),
  c(mlfts_f_evr[h_16, "score"],  mean(mlfts_f_evr[h_16, "score"]),  median(mlfts_f_evr[h_16, "score"])),
  c(fanova_f_evr[h_16, "score"], mean(fanova_f_evr[h_16, "score"]), median(fanova_f_evr[h_16, "score"])),
  c(hdfpca_f_evr[h_16, "score"], mean(hdfpca_f_evr[h_16, "score"]), median(hdfpca_f_evr[h_16, "score"])),
  c(ufts_f_k6[h_16, "score"],    mean(ufts_f_k6[h_16, "score"]),    median(ufts_f_k6[h_16, "score"])),
  c(mfts_f_k6[h_16, "score"],    mean(mfts_f_k6[h_16, "score"]),    median(mfts_f_k6[h_16, "score"])),
  c(mlfts_f_k6[h_16, "score"],   mean(mlfts_f_k6[h_16, "score"]),   median(mlfts_f_k6[h_16, "score"])),
  c(fanova_f_k6[h_16, "score"],  mean(fanova_f_k6[h_16, "score"]),  median(fanova_f_k6[h_16, "score"]))
)

# --- INTERVAL SCORE MALE ---
score_m_matrix <- cbind(
  c(ufts_m_evr[h_16, "score"],   mean(ufts_m_evr[h_16, "score"]),   median(ufts_m_evr[h_16, "score"])),
  c(mfts_m_evr[h_16, "score"],   mean(mfts_m_evr[h_16, "score"]),   median(mfts_m_evr[h_16, "score"])),
  c(mlfts_m_evr[h_16, "score"],  mean(mlfts_m_evr[h_16, "score"]),  median(mlfts_m_evr[h_16, "score"])),
  c(fanova_m_evr[h_16, "score"], mean(fanova_m_evr[h_16, "score"]), median(fanova_m_evr[h_16, "score"])),
  c(hdfpca_m_evr[h_16, "score"], mean(hdfpca_m_evr[h_16, "score"]), median(hdfpca_m_evr[h_16, "score"])),
  c(ufts_m_k6[h_16, "score"],    mean(ufts_m_k6[h_16, "score"]),    median(ufts_m_k6[h_16, "score"])),
  c(mfts_m_k6[h_16, "score"],    mean(mfts_m_k6[h_16, "score"]),    median(mfts_m_k6[h_16, "score"])),
  c(mlfts_m_k6[h_16, "score"],   mean(mlfts_m_k6[h_16, "score"]),   median(mlfts_m_k6[h_16, "score"])),
  c(fanova_m_k6[h_16, "score"],  mean(fanova_m_k6[h_16, "score"]),  median(fanova_m_k6[h_16, "score"]))
)

# Helper function to generate formatted table rows
make_row <- function(metric_lbl, sex_lbl, h_val, vals, is_score = FALSE) {
  fmt_fn <- if (is_score) fmt_score else fmt_cpd
  val_str <- paste(sapply(vals, fmt_fn), collapse = " & ")
  
  m_str <- if (metric_lbl != "") metric_lbl else ""
  s_str <- if (sex_lbl != "") sex_lbl else ""
  h_str <- if (!is.na(h_val)) as.character(h_val) else ""
  
  return(paste0(m_str, " & ", s_str, " & ", h_str, " & ", val_str, " \\\\"))
}

# 3. PRINT UPDATED LATEX TEMPLATE FOR TABLE C.1
# ------------------------------------------------------------------------------
cat('\\begin{center}
\\renewcommand*{\\arraystretch}{0.934}
\\tabcolsep 0.062in
\\begin{small}
\\begin{longtable}{@{}lllrrrrrrrrrr@{}}
\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the interval forecast accuracy, measured by the CPD and interval score, for the functional time-series models at the significance level $\\alpha=0.2$. The number of components is determined by the EVR criterion or $K=6$; for each, we highlight the method with the smallest overall CPD and interval score values.}\\label{tab:appendix_3} \\\\
\\toprule
& & \\multicolumn{5}{c}{EVR} & & \\multicolumn{4}{c}{$K=6$} \\\\
Metric & Sex & $h$ & UFTS & MFTS & MLFTS & FANOVA & HDFPCA & UFTS & MFTS & MLFTS & FANOVA \\\\
\\midrule
\\endfirsthead
\\toprule
& & \\multicolumn{5}{c}{EVR} & & \\multicolumn{4}{c}{$K=6$} \\\\
Metric & Sex & $h$ & UFTS & MFTS & MLFTS & FANOVA & HDFPCA & UFTS & MFTS & MLFTS & FANOVA \\\\
\\midrule
\\endhead
\\midrule
\\multicolumn{12}{r}{{Continued on next page}} \\\\
\\endfoot
\\endlastfoot
%----------------- CPD Female -----------------\n')

# CPD Female Rows 1-16
for (h in 1:16) {
  m_lbl <- if (h == 1) "CPD$_{0.2}$" else ""
  s_lbl <- if (h == 1) "F" else ""
  cat(make_row(m_lbl, s_lbl, h, cpd_f_matrix[h, ], is_score = FALSE), "\n")
}
cat('\\cmidrule{1-12}\n')
cat(make_row("", "", "Mean", cpd_f_matrix[17, ], is_score = FALSE), "\n")
cat(make_row("", "", "Median", cpd_f_matrix[18, ], is_score = FALSE), "\n")
cat('\\cmidrule{1-12}\n%----------------- CPD Male -----------------\n')

# CPD Male Rows 1-16
for (h in 1:16) {
  s_lbl <- if (h == 1) "M" else ""
  cat(make_row("", s_lbl, h, cpd_m_matrix[h, ], is_score = FALSE), "\n")
}
cat('\\cmidrule{1-12}\n')
cat(make_row("", "", "Mean", cpd_m_matrix[17, ], is_score = FALSE), "\n")
cat(make_row("", "", "Median", cpd_m_matrix[18, ], is_score = FALSE), "\n")
cat('\\cmidrule{1-12}\n%----------------- Interval Score Female -----------------\n\\cmidrule{1-12}\n')

# Interval Score Female Rows 1-16
for (h in 1:16) {
  m_lbl <- if (h == 1) "$\\overline{S}_{0.2}$" else ""
  s_lbl <- if (h == 1) "F" else ""
  cat(make_row(m_lbl, s_lbl, h, score_f_matrix[h, ], is_score = TRUE), "\n")
}
cat('\\cmidrule{1-12}\n')
cat(make_row("", "", "Mean", score_f_matrix[17, ], is_score = TRUE), "\n")
cat(make_row("", "", "Median", score_f_matrix[18, ], is_score = TRUE), "\n")
cat('\\cmidrule{1-12}\n%----------------- Interval Score Male -----------------\n')

# Interval Score Male Rows 1-16
for (h in 1:16) {
  s_lbl <- if (h == 1) "M" else ""
  cat(make_row("", s_lbl, h, score_m_matrix[h, ], is_score = TRUE), "\n")
}
cat('\\cmidrule{1-12}\n')
cat(make_row("", "", "Mean", score_m_matrix[17, ], is_score = TRUE), "\n")
cat(make_row("", "", "Median", score_m_matrix[18, ], is_score = TRUE), "\n")

cat('\\bottomrule
\\end{longtable}
\\end{small}
\\end{center}\n')

