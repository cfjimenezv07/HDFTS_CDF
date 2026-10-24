###############################################################################
# Script: TableC1_ConformalAnalysis.R
# Description: Reads saved interval forecast error (conformal) .rds files and
#               prints Table C1 in the LaTeX format used in the paper, while
#               simultaneously exporting a clean .tex file.
###############################################################################
source(here("load_packages.R"))

# Dynamic directory path for Conformal Interval Forecast results
dir_results <- here("results", "Results_Interval", "conformal")
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

# 1. READ ALL RDS FILES FOR CONFORMAL INTERVAL FORECASTS
# ------------------------------------------------------------------------------
# --- EVR CRITERION ---
# UFTS
ufts_f_evr <- readRDS(file.path(dir_results, "int_fore_subnational_err_F_EVR_conformal_mean.rds"))
ufts_m_evr <- readRDS(file.path(dir_results, "int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# MFTS
mfts_f_evr <- readRDS(file.path(dir_results, "MFTS_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
mfts_m_evr <- readRDS(file.path(dir_results, "MFTS_int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# MLFTS
mlfts_f_evr <- readRDS(file.path(dir_results, "MLFTS_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
mlfts_m_evr <- readRDS(file.path(dir_results, "MLFTS_int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# FANOVA
fanova_f_evr <- readRDS(file.path(dir_results, "FANOVA_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
fanova_m_evr <- readRDS(file.path(dir_results, "FANOVA_int_fore_subnational_err_M_EVR_conformal_mean.rds"))

# HDFPCA
hdfpca_f_evr <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_F_EVR_conformal_mean.rds"))
hdfpca_m_evr <- readRDS(file.path(dir_results, "hdfpca_int_fore_subnational_err_M_EVR_conformal_mean.rds"))


# --- K = 6 CRITERION ---
# UFTS
ufts_f_k6 <- readRDS(file.path(dir_results, "int_fore_subnational_err_F_K6_conformal_mean.rds"))
ufts_m_k6 <- readRDS(file.path(dir_results, "int_fore_subnational_err_M_K6_conformal_mean.rds"))

# MFTS
mfts_f_k6 <- readRDS(file.path(dir_results, "MFTS_int_fore_subnational_err_F_K6_conformal_mean.rds"))
mfts_m_k6 <- readRDS(file.path(dir_results, "MFTS_int_fore_subnational_err_M_K6_conformal_mean.rds"))

# MLFTS
mlfts_f_k6 <- readRDS(file.path(dir_results, "MLFTS_int_fore_subnational_err_F_K6_conformal_mean.rds"))
mlfts_m_k6 <- readRDS(file.path(dir_results, "MLFTS_int_fore_subnational_err_M_K6_conformal_mean.rds"))

# FANOVA
fanova_f_k6 <- readRDS(file.path(dir_results, "FANOVA_int_fore_subnational_err_F_K6_conformal_mean.rds"))
fanova_m_k6 <- readRDS(file.path(dir_results, "FANOVA_int_fore_subnational_err_M_K6_conformal_mean.rds"))

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
  fmt_fn  <- if (is_score) fmt_score else fmt_cpd
  val_str <- paste(sapply(vals, fmt_fn), collapse = " & ")
  
  m_str <- sprintf("%-18s", metric_lbl)
  s_str <- sprintf("%-3s", sex_lbl)
  h_str <- if (!is.na(h_val)) sprintf("%4s", as.character(h_val)) else sprintf("%-6s", h_val)
  
  return(paste0(m_str, " & ", s_str, " & ", h_str, " & ", val_str, " \\\\"))
}

# 3. GENERATE LATEX TABLE CODE FOR TABLE C.1
# ------------------------------------------------------------------------------
generate_table_code <- function() {
  lines <- c()
  
  lines <- c(lines, '\\begin{center}')
  lines <- c(lines, '\\renewcommand*{\\arraystretch}{0.934}')
  lines <- c(lines, '\\tabcolsep 0.062in')
  lines <- c(lines, '\\begin{small}')
  lines <- c(lines, '\\begin{longtable}{@{}lllrrrrrrrrrr@{}}')
  lines <- c(lines, '\\caption{\\small Averaged across 47 prefectures, we evaluate and compare the interval forecast accuracy, measured by the CPD and interval score, for the functional time series models at the significance level $\alpha=0.2$. The number of components is determined by the EVR criterion or $K=6$; for each, we highlight the method with the smallest overall CPD and interval score values.}\\label{tab:appendix_3} \\\\')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '& & \\multicolumn{5}{c}{EVR} & & \\multicolumn{4}{c}{$K=6$} \\\\')
  lines <- c(lines, '\\cmidrule{4-8} \\cmidrule{9-12}')
  lines <- c(lines, 'Metric & Sex & $h$ & UFTS & MFTS & MLFTS & FANOVA & HDFPCA & UFTS & MFTS & MLFTS & FANOVA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endfirsthead')
  lines <- c(lines, '\\toprule')
  lines <- c(lines, '& & \\multicolumn{5}{c}{EVR} & & \\multicolumn{4}{c}{$K=6$} \\\\')
  lines <- c(lines, '\\cmidrule{4-8} \\cmidrule{9-12}')
  lines <- c(lines, 'Metric & Sex & $h$ & UFTS & MFTS & MLFTS & FANOVA & HDFPCA & UFTS & MFTS & MLFTS & FANOVA \\\\')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\endhead')
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '\\multicolumn{12}{r}{{Continued on next page}} \\\\')
  lines <- c(lines, '\\endfoot')
  lines <- c(lines, '\\endlastfoot')
  lines <- c(lines, '%----------------- CPD Female -----------------')
  
  # CPD Female Rows 1-16
  for (h in 1:16) {
    m_lbl <- if (h == 1) "CPD$_{0.2}$" else ""
    s_lbl <- if (h == 1) "F" else ""
    lines <- c(lines, make_row(m_lbl, s_lbl, h, cpd_f_matrix[h, ], is_score = FALSE))
  }
  lines <- c(lines, '\\cmidrule{3-12}')
  lines <- c(lines, make_row("", "", "Mean", cpd_f_matrix[17, ], is_score = FALSE))
  lines <- c(lines, make_row("", "", "Median", cpd_f_matrix[18, ], is_score = FALSE))
  lines <- c(lines, '\\cmidrule{1-12}')
  lines <- c(lines, '%----------------- CPD Male -----------------')
  
  # CPD Male Rows 1-16
  for (h in 1:16) {
    s_lbl <- if (h == 1) "M" else ""
    lines <- c(lines, make_row("", s_lbl, h, cpd_m_matrix[h, ], is_score = FALSE))
  }
  lines <- c(lines, '\\cmidrule{3-12}')
  lines <- c(lines, make_row("", "", "Mean", cpd_m_matrix[17, ], is_score = FALSE))
  lines <- c(lines, make_row("", "", "Median", cpd_m_matrix[18, ], is_score = FALSE))
  lines <- c(lines, '\\midrule')
  lines <- c(lines, '%----------------- Interval Score Female -----------------')
  
  # Interval Score Female Rows 1-16
  for (h in 1:16) {
    m_lbl <- if (h == 1) "$\\overline{S}_{0.2}$" else ""
    s_lbl <- if (h == 1) "F" else ""
    lines <- c(lines, make_row(m_lbl, s_lbl, h, score_f_matrix[h, ], is_score = TRUE))
  }
  lines <- c(lines, '\\cmidrule{3-12}')
  lines <- c(lines, make_row("", "", "Mean", score_f_matrix[17, ], is_score = TRUE))
  lines <- c(lines, make_row("", "", "Median", score_f_matrix[18, ], is_score = TRUE))
  lines <- c(lines, '\\cmidrule{1-12}')
  lines <- c(lines, '%----------------- Interval Score Male -----------------')
  
  # Interval Score Male Rows 1-16
  for (h in 1:16) {
    s_lbl <- if (h == 1) "M" else ""
    lines <- c(lines, make_row("", s_lbl, h, score_m_matrix[h, ], is_score = TRUE))
  }
  lines <- c(lines, '\\cmidrule{3-12}')
  lines <- c(lines, make_row("", "", "Mean", score_m_matrix[17, ], is_score = TRUE))
  lines <- c(lines, make_row("", "", "Median", score_m_matrix[18, ], is_score = TRUE))
  
  lines <- c(lines, '\\bottomrule')
  lines <- c(lines, '\\end{longtable}')
  lines <- c(lines, '\\end{small}')
  lines <- c(lines, '\\end{center}')
  
  return(lines)
}

# 4. PRINT TO CONSOLE & EXPORT FILE
# ------------------------------------------------------------------------------
latex_table_lines <- generate_table_code()

# Print directly to R console
cat(paste(latex_table_lines, collapse = "\n"), "\n")

# Export directly to results/Results_Interval/conformal/TableC1.tex
output_tex_path <- file.path(dir_results, "TableC1.tex")
writeLines(latex_table_lines, con = output_tex_path)
message("\nSuccess: Table C1 LaTeX exported directly to: ", output_tex_path)