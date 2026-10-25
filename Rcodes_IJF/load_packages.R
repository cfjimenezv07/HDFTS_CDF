################################################################################
# Script: load_packages.R
# Description: Installs and loads all required package dependencies in the session.
################################################################################

# Prevent 'rgl' / XQuartz crash on macOS when loading Compositional
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = FALSE)

# Ensure CRAN mirror is set globally for non-interactive sessions
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Ensure 'here' is installed and loaded first
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}
library(here)

# Base packages list
packages <- c(
  "ftsa",
  "LaplacesDemon",
  "flexmix",
  "psych",
  "easyCODA",
  "MortalityLaws",
  "DescTools",
  "xtable",
  "tidyverse",
  "ggplot2",
  "RColorBrewer",
  "rlist",
  "hdftsa",
  "demography",
  "vars",
  "Compositional",
  "transport",
  "wwntests",
  "FTSgof",
  "forecast",
  "boot",
  "foreach",
  "dplyr",
  "viridis",
  "ggpubr",
  "compositions",
  "MCS",
  "doParallel",
  "parallel",
  "doRNG",
  "here"
)

# Add OS-specific package for Unix/Linux/macOS
if (.Platform$OS.type == "unix") {
  packages <- c(packages, "doMC")
}

# Function to check, install with explicit repos, AND load packages
load_and_check_packages <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  } else {
    message("All required packages are installed.")
  }
  
  # Load all packages into the session
  invisible(lapply(pkgs, library, character.only = TRUE))
}

load_and_check_packages(packages)