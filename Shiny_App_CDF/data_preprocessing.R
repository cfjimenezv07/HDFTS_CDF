#Making data ready for the pfe table

#gettig the names regarding the prefectures
names.Japan <- sapply(readRDS("names/names_prefectures.rds"), function(x) x)
names <- list( Japan = names.Japan)

get_pfe_dataset <- function(.country, .method, .metric, .gender, .pca) {
  #reading the data
  path <- paste0("./datasets_shiny_app/PFE/",.country, "/", 
                 paste(.method, .metric, .gender, .pca, sep = "_"), ".rds")
  df <- readRDS(path)
  colnames(df) <- names[[.country]] #labeling the departments or states
  df <- as.data.frame(df)
  #adding the identifiers such
  df$method <- .method
  df$pca <- .pca
  df$metric <- .metric
  df$h <- 1:16
  df$gender <- .gender
  df$country <- .country
  #pivoting dataset
  df %>% 
    pivot_longer(-(method:country), names_to = "pol_division")
}

#generating data
country.list <- as.list(rep(c("Japan"), each = 4))
method.list <- as.list(rep(c("UFTS", "MFTS","MLFTS","FANOVA","HDFPCA"), each = 4))
metric.list <- as.list(rep(c("JSD", "KLD"), each = 4))
gender.list <- as.list(rep(c("female", "male"), each = 4))
pca.list <- as.list(rep(c("EVR", "K"), each = 4))


lista <- expand_grid(country = c("Japan"), method = c("UFTS", "MFTS","MLFTS","FANOVA","HDFPCA"),
                     metric = c("JSD", "KLD"), gender = c("male", "female"), pca = c("EVR", "K")) %>%
  apply(2, function(x) as.list(x))

pfe.result <- mapply(get_pfe_dataset, .country = lista$country, .method = lista$method, 
                     .metric = lista$metric, .gender = lista$gender, .pca = lista$pca, SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .) %>%
  mutate(value = round(value, 4), method = method)

View(pfe.result)
# saveRDS(pfe.result, "./pfe_data.rds")





#ife dataset

aux <- function(i, df, .country, .gender, .method, .coverage, .pca) {
  temp <- df[i, , ]
  temp[, -3] <- round(temp[, -3], 3)
  temp[, 3] <- round(temp[, 3], 0)
  colnames(temp) <- c("ECP","CPD", "IS")
  temp <- as.data.frame(temp)
  temp$method <- .method 
  temp$pca <- .pca
  temp$h <- 1:15
  temp$coverage <- .coverage
  temp$gender <- .gender
  temp$country <- .country
  temp$pol_division = names[[.country]][i]
  temp
}



get_ife_dataset <- function(.country, .method, .gender, .coverage, .pca) {
  #reading the data
  path <- paste0("./datasets_shiny_app/IFE/",.country, "/", 
                 paste(.method, .gender, .pca,.coverage, sep = "_"), ".rds")
  df <- readRDS(path)
  do.call(rbind.data.frame, 
          lapply(1:47, aux, df, .country, .gender, .method, .coverage, .pca))
}

lista <- expand_grid(country = c("Japan"), method = c("UFTS", "MFTS","MLFTS","FANOVA","HDFPCA"),
                     gender = c("male", "female"), coverage = c(80, 95), pca = c("EVR", "K")) %>%
  apply(2, function(x) as.list(x))

ife.result <- mapply(get_ife_dataset, .country = lista$country, .method = lista$method, 
                     .gender = lista$gender, .coverage = lista$coverage, 
                     .pca = lista$pca, SIMPLIFY = F) %>%
  do.call(rbind.data.frame, .) %>%
  mutate(method = method) %>%
  pivot_longer(ECP:IS, names_to = "metric") %>%
  relocate(metric, everything())


saveRDS(ife.result, "./ife_data.rds")

# View(ife.result)
