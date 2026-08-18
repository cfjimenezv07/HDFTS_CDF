##################
# load R packages
##################
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = FALSE)

# Package names
packages <- c(
  "ftsa",
  "LaplacesDemon",
  "flexmix",
  "psych",
  "easyCODA",
  "doMC",
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
  "here"
)

# Function to check, install, AND load packages
load_and_check_packages <- function(packages) {
  missing <- packages[!packages %in% rownames(installed.packages())]
  
  if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "))
    install.packages(missing)
  } else {
    message("All packages are installed.")
  }
  
  # Load all packages into the session
  invisible(lapply(packages, library, character.only = TRUE))
}

load_and_check_packages(packages)

