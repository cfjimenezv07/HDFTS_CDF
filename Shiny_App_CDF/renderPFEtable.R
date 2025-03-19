
library(tidyverse)
pfe.data <- readRDS('./pfe_data.Rds')

gen_pfe_table <- function(.country, .division, .gender, .pca) {
  df <- pfe.data %>% 
    filter(country == .country, pol_division == .division, gender == tolower(.gender),
           pca == .pca) %>%
    select(c(-(gender:pol_division), -pca))
  ufts.table <- df %>% 
    filter(method == "UFTS") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(h, KLD, JSD)
  
  mfts.table <- df %>% 
    filter(method == "MFTS") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  mlfts.table <- df %>% 
    filter(method == "MLFTS") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  fanova.table <- df %>% 
    filter(method == "FANOVA") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  hdfpca.table <- df %>% 
    filter(method == "HDFPCA") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  temp <- cbind(ufts.table, mfts.table,mlfts.table,fanova.table,hdfpca.table)
  last.row1  <- format(round(colMeans(temp[1:6,]), 3), nsmall = 2)
  last.row2 <- format(round(apply(temp[1:6,],2,median), 3), nsmall = 2)

  last.row3  <- format(round(colMeans(temp[7:12,]), 3), nsmall = 2)
  last.row4 <- format(round(apply(temp[7:12,],2,median), 3), nsmall = 2)

  last.row5  <- format(round(colMeans(temp[13:16,]), 3), nsmall = 2)
  last.row6 <- format(round(apply(temp[13:16,],2,median), 3), nsmall = 2)

  temp <- temp  %>%
    apply(2, function(x) format(round(x, 3), nsmall = 2)) %>% 
    as.data.frame()
  temp$h <- as.numeric(temp$h)
  last.row1[1]  <- "Mean"
  last.row2[1] <- "Median"
  last.row3[1]  <- "Mean"
  last.row4[1] <- "Median"
  last.row5[1]  <- "Mean"
  last.row6[1] <- "Median"
  rbind(temp[1:6, ], last.row1,last.row2, 
        temp[7:12, ], last.row3,last.row4,
        temp[13:16, ], last.row5,last.row6)
  

}



