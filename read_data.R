##################
# load R packages
##################

require(psych)
require(flexmix)
require(RColorBrewer)
require(LaplacesDemon)
require(ftsa)
require(easyCODA)
require(doMC)
require(MortalityLaws)
require(DescTools)
require(xtable)
require(transport)
require(Compositional)
require(compositions)
require(dplyr)

setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/")
# source("auxiliary.R")

# HLS where the HLS results are stored
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/HLS/"
# FANOVA where the FANOVA results are stored
dir.r1 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/FANOVA/"


######################
# model HDFTS of CDFs
######################

# read Japanese Subnational Human Mortality Data

state = c("Hokkaido", 
          "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
          "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
          "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
          "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

# change file directory

setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/updated_data")

# Define age and year ranges
ages <- 0:110
n_age <- length(ages)
years <- 1975:2023
n_year <- length(years)

# Create empty lists to store results
female_prefecture_qx <- list()
male_prefecture_qx   <- list()
total_prefecture_qx  <- list()

# Loop over all prefectures
for (ik in seq_along(state)) {
  
  # -------------------- FEMALE --------------------
  female_data <- read.table(
    paste0("female_prefecture_", ik, ".txt"),
    header = TRUE, skip = 2
  )
  # Keep only years >= 1975
  female_data <- female_data[female_data$Year >= 1975, ]
  
  # Reshape into n_year x n_age matrix (years x ages)
  female_prefecture_qx[[ik]] <- t(matrix(female_data$qx, n_age, n_year))
  
  
  # -------------------- MALE --------------------
  male_data <- read.table(
    paste0("male_prefecture_", ik, ".txt"),
    header = TRUE, skip = 2
  )
  male_data <- male_data[male_data$Year >= 1975, ]
  
  male_prefecture_qx[[ik]] <- t(matrix(male_data$qx, n_age, n_year))
  
  
  # -------------------- TOTAL --------------------
  total_data <- read.table(
    paste0("total_prefecture_", ik, ".txt"),
    header = TRUE, skip = 2
  )
  total_data <- total_data[total_data$Year >= 1975, ]
  
  total_prefecture_qx[[ik]] <- t(matrix(total_data$qx, n_age, n_year))
  
  # Print progress
  print(paste("Finished prefecture:", ik))
}

# Optional sanity check
dim(female_prefecture_qx[[1]])
range(female_data$Year)



# Japanese subnational data

female_prefecture_dx = male_prefecture_dx = total_prefecture_dx = list()
for(iw in 1:length(state))
{
    female_prefecture_dum = male_prefecture_dum = total_prefecture_dum = matrix(NA, n_year, n_age)
    for(ij in 1:n_year)
    {
        # set radix (normalising to 1 or keep it as 10^5)
        start_pop_female = start_pop_male = start_pop_total = 10^5
        for(ik in 1:n_age)
        {
            female_prefecture_dum[ij,ik] = (female_prefecture_qx[[iw]])[ij,ik] * start_pop_female
            start_pop_female = start_pop_female - female_prefecture_dum[ij,ik]
            
            male_prefecture_dum[ij,ik] = (male_prefecture_qx[[iw]])[ij,ik] * start_pop_male
            start_pop_male = start_pop_male - male_prefecture_dum[ij,ik]
            
            total_prefecture_dum[ij,ik] = (total_prefecture_qx[[iw]])[ij,ik] * start_pop_total
            start_pop_total = start_pop_total - total_prefecture_dum[ij,ik]
        }
    }
    female_prefecture_dx[[iw]] = female_prefecture_dum
    male_prefecture_dx[[iw]]   = male_prefecture_dum
    total_prefecture_dx[[iw]]  = total_prefecture_dum
    rm(female_prefecture_dum); rm(male_prefecture_dum); rm(total_prefecture_dum)
    print(iw); rm(iw)
}

saveRDS(female_prefecture_dx, file = paste0(dir.r1, "female_prefecture_dx.rds"))
saveRDS(male_prefecture_dx,   file = paste0(dir.r1, "male_prefecture_dx.rds"))
saveRDS(total_prefecture_dx,  file = paste0(dir.r1, "total_prefecture_dx.rds"))


########################################
# Average 47 matrices (49×111 each)
########################################

# Helper: average a list of matrices element-wise
mean_matrix_list <- function(mat_list) {
  Reduce("+", mat_list) / length(mat_list)
}

# Compute national-level matrices
Japan_female_pop <- mean_matrix_list(female_prefecture_dx)
Japan_male_pop   <- mean_matrix_list(male_prefecture_dx)
Japan_total_pop   <- mean_matrix_list(total_prefecture_dx)



saveRDS(Japan_female_pop, file = paste0(dir.r1, "Japan_female_pop.rds"))
saveRDS(Japan_male_pop,   file = paste0(dir.r1, "Japan_male_pop.rds"))
saveRDS(Japan_total_pop,   file = paste0(dir.r1, "Japan_total_pop.rds"))



# Japan national data

# Japan_female_dx = t(matrix(read.table("JPN_female_lt.txt", header = TRUE)$dx, n_age, n_year))
# Japan_male_dx   = t(matrix(read.table("JPN_male_lt.txt", header = TRUE)$dx,   n_age, n_year))
# 
# Japan_female_qx = t(matrix(read.table("JPN_female_lt.txt", header = TRUE)$qx, n_age, n_year))
# Japan_male_qx   = t(matrix(read.table("JPN_male_lt.txt", header = TRUE)$qx,   n_age, n_year))
# Japan_total_qx  = t(matrix(read.table("JPN_total_lt.txt", header = TRUE)$qx,  n_age, n_year))
# 
# Japan_female_pop = Japan_male_pop = Japan_total_pop = matrix(NA, n_year, n_age)
# for(ij in 1:n_year)
# {
#     start_pop_female = start_pop_male = start_pop_total = 10^5
#     for(ik in 1:n_age)
#     {
#         Japan_female_pop[ij,ik] = Japan_female_qx[ij,ik] * start_pop_female
#         start_pop_female = start_pop_female - Japan_female_pop[ij,ik]
#         
#         Japan_male_pop[ij,ik] = Japan_male_qx[ij,ik] * start_pop_male
#         start_pop_male = start_pop_male - Japan_male_pop[ij,ik]
#         
#         Japan_total_pop[ij,ik] = Japan_total_qx[ij,ik] * start_pop_total
#         start_pop_total = start_pop_total - Japan_total_pop[ij,ik]
#         rm(ik)
#     }
#     rm(ij)
# }
# rownames(Japan_female_pop) = rownames(Japan_male_pop) = rownames(Japan_total_pop) = years
# colnames(Japan_female_pop) = colnames(Japan_male_pop) = colnames(Japan_total_pop) = ages
# 
