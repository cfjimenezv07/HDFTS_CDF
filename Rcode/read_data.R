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

setwd("/Users/hanlinshang/Dropbox/Todos/HDFTS_CDF_subnational/code")
source("auxiliary.R")

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

setwd("/Users/hanlinshang/Dropbox/Todos/HDFTS_CDF_subnational/data/Japanese_data/Lifetables/")

ages = 0:110
n_age = length(ages)
years = 1975:2022
n_year = length(years)

female_prefecture_qx = male_prefecture_qx = total_prefecture_qx = list()
for(ik in 1:length(state))
{
    female_prefecture_qx[[ik]] = t(matrix(read.table(paste(paste("female_prefecture_", ik, sep = ""), ".txt", sep = ""), header = TRUE)$qx, n_age, n_year))
    male_prefecture_qx[[ik]]   = t(matrix(read.table(paste(paste("male_prefecture_", ik, sep = ""), ".txt", sep = ""),   header = TRUE)$qx, n_age, n_year))
    total_prefecture_qx[[ik]]  = t(matrix(read.table(paste(paste("total_prefecture_", ik, sep = ""), ".txt", sep = ""),  header = TRUE)$qx, n_age, n_year))
    print(ik); rm(ik)
}

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

# Japan national data

Japan_female_dx = t(matrix(read.table("JPN_female_lt.txt", header = TRUE)$dx, n_age, n_year))
Japan_male_dx   = t(matrix(read.table("JPN_male_lt.txt", header = TRUE)$dx,   n_age, n_year))

Japan_female_qx = t(matrix(read.table("JPN_female_lt.txt", header = TRUE)$qx, n_age, n_year))
Japan_male_qx   = t(matrix(read.table("JPN_male_lt.txt", header = TRUE)$qx,   n_age, n_year))
Japan_total_qx  = t(matrix(read.table("JPN_total_lt.txt", header = TRUE)$qx,  n_age, n_year))

Japan_female_pop = Japan_male_pop = Japan_total_pop = matrix(NA, n_year, n_age)
for(ij in 1:n_year)
{
    start_pop_female = start_pop_male = start_pop_total = 10^5
    for(ik in 1:n_age)
    {
        Japan_female_pop[ij,ik] = Japan_female_qx[ij,ik] * start_pop_female
        start_pop_female = start_pop_female - Japan_female_pop[ij,ik]
        
        Japan_male_pop[ij,ik] = Japan_male_qx[ij,ik] * start_pop_male
        start_pop_male = start_pop_male - Japan_male_pop[ij,ik]
        
        Japan_total_pop[ij,ik] = Japan_total_qx[ij,ik] * start_pop_total
        start_pop_total = start_pop_total - Japan_total_pop[ij,ik]
        rm(ik)
    }
    rm(ij)
}
rownames(Japan_female_pop) = rownames(Japan_male_pop) = rownames(Japan_total_pop) = years
colnames(Japan_female_pop) = colnames(Japan_male_pop) = colnames(Japan_total_pop) = ages

