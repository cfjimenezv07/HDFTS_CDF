###########

# heatmaps

###########
# Where heatmaps will be saved
dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/Plots/"

source("Heatmaps.R")

n_methods = 9

KLD_fore_subnational_err_F_array = array(NA, dim = c(17, 47, n_methods), 
                                         
                                         dimnames = list(1:17, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA",
                                                                       
                                                                       "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)", "FANOVA(K=6)")))

KLD_fore_subnational_err_F_array[,,1] = point_fore_subnational_err_F_EVR_ETS

KLD_fore_subnational_err_F_array[,,2] = point_fore_subnational_err_MFTS_F_EVR_ETS

KLD_fore_subnational_err_F_array[,,3] = point_fore_subnational_err_MLFTS_F_EVR_ETS

KLD_fore_subnational_err_F_array[,,4] = errors_KLD_female_EVR

KLD_fore_subnational_err_F_array[,,5] = point_fore_subnational_err_HDFPCA_F_CDF



KLD_fore_subnational_err_F_array[,,6] = point_fore_subnational_err_F_K6_ETS

KLD_fore_subnational_err_F_array[,,7] = point_fore_subnational_err_MFTS_F_K6_ETS

KLD_fore_subnational_err_F_array[,,8] = point_fore_subnational_err_MLFTS_F_K6_ETS

KLD_fore_subnational_err_F_array[,,9] = errors_KLD_female_K



KLD_fore_subnational_err_F_ranking_mat = matrix(0, nrow = 17, ncol = n_methods)

colnames(KLD_fore_subnational_err_F_ranking_mat) = paste("M", 1:n_methods, sep = "")

rownames(KLD_fore_subnational_err_F_ranking_mat) = 1:17



for(ih in 1:17)
  
{
  
  temp_rank = table(apply(KLD_fore_subnational_err_F_array[ih,,],1, which.min))
  
  if(length(temp_rank) == 9)
    
  {
    
    KLD_fore_subnational_err_F_ranking_mat[ih,] = temp_rank
    
  }
  
  else
    
  {
    
    temp_rank_extend = rep(0, 9)
    
    temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
    
    KLD_fore_subnational_err_F_ranking_mat[ih,] = temp_rank_extend
    
    rm(temp_rank_extend)
    
  }
  
  rm(temp_rank)
  
}



# create a heatmap



KLD_fore_subnational_err_F_ranking_mat = data.frame(Horizon = 1:17, KLD_fore_subnational_err_F_ranking_mat)



KLD_female = KLD_fore_subnational_err_F_ranking_mat %>% 
  
  as_tibble() %>%
  
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  
  mutate(
    
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:17),
    
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = ""))
    
  )



ggsave(file.path(dir.p, "Fig_4a.png"))

ggplot(KLD_female, aes(Model, ordered(Horizon, levels = 17:1))) +
  
  geom_tile(aes(fill = count)) +
  
  geom_text(aes(label = count)) +
  
  scale_x_discrete(position = "top") +
  
  theme(legend.position = "none") + 
  
  scale_fill_gradient2(high="green",mid="white",low="red", 
                       
                       na.value="yellow", midpoint = 5) + 
  
  ylab("Forecast horizon") +
  
  scale_x_discrete(labels = c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", "HDFPCA\n",
                              
                              "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)")) + 
  
  theme_bw()

dev.off()



# M



KLD_fore_subnational_err_M_array = array(NA, dim = c(17, 47, 9), 
                                         
                                         dimnames = list(1:17, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA", "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)", "FANOVA (K=6)")))

KLD_fore_subnational_err_M_array[,,1] = point_fore_subnational_err_M_EVR_ETS

KLD_fore_subnational_err_M_array[,,2] = point_fore_subnational_err_MFTS_M_EVR_ETS

KLD_fore_subnational_err_M_array[,,3] = point_fore_subnational_err_MLFTS_M_EVR_ETS

KLD_fore_subnational_err_M_array[,,4] = errors_KLD_male_EVR

KLD_fore_subnational_err_M_array[,,5] = point_fore_subnational_err_HDFPCA_M_CDF



KLD_fore_subnational_err_M_array[,,6] = point_fore_subnational_err_M_K6_ETS

KLD_fore_subnational_err_M_array[,,7] = point_fore_subnational_err_MFTS_M_K6_ETS

KLD_fore_subnational_err_M_array[,,8] = point_fore_subnational_err_MLFTS_M_K6_ETS

KLD_fore_subnational_err_M_array[,,9] = errors_KLD_male_K



KLD_fore_subnational_err_M_ranking_mat = matrix(0, nrow = 17, ncol = 9)

colnames(KLD_fore_subnational_err_M_ranking_mat) = paste("M", 1:9, sep = "")

rownames(KLD_fore_subnational_err_M_ranking_mat) = 1:17



for(ih in 1:17)
  
{
  
  temp_rank = table(apply(KLD_fore_subnational_err_M_array[ih,,],1, which.min))
  
  if(length(temp_rank) == 9)
    
  {
    
    KLD_fore_subnational_err_M_ranking_mat[ih,] = temp_rank
    
  }
  
  else
    
  {
    
    temp_rank_extend = rep(0, 9)
    
    temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
    
    KLD_fore_subnational_err_M_ranking_mat[ih,] = temp_rank_extend
    
    rm(temp_rank_extend)
    
  }
  
  rm(temp_rank)
  
}



# create a heatmap



KLD_fore_subnational_err_M_ranking_mat = data.frame(Horizon = 1:17, KLD_fore_subnational_err_M_ranking_mat)



KLD_male = KLD_fore_subnational_err_M_ranking_mat %>% 
  
  as_tibble() %>%
  
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  
  mutate(
    
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:17),
    
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = ""))
    
  )



ggsave(file.path(dir.p, "Fig_4b.png"))


ggplot(KLD_male, aes(Model, ordered(Horizon, levels = 17:1))) +
  
  geom_tile(aes(fill = count)) +
  
  geom_text(aes(label = count)) +
  
  scale_x_discrete(position = "top") +
  
  theme(legend.position = "none") + 
  
  scale_fill_gradient2(high="green",mid="white",low="red", 
                       
                       na.value="yellow", midpoint = 5) + 
  
  ylab("Forecast horizon") +
  
  scale_x_discrete(labels = c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", "HDFPCA", 
                              
                              "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)"))+
  
  theme_bw()

dev.off()

