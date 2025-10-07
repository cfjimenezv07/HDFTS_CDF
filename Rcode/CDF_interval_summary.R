######
# ECP
######

## F

ECP_F_alpha_0.8 = cbind(int_fore_subnational_err_F_EVR_ETS_mean[,1], 
                        MFTS_int_fore_subnational_err_F_EVR_ETS_mean[,1],
                        MLFTS_int_fore_subnational_err_F_EVR_ETS_mean[,1],
                        horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS[,1],
                        
                        int_fore_subnational_err_F_K6_ETS_mean[,1],
                        MFTS_int_fore_subnational_err_F_K6_ETS_mean[,1],
                        MLFTS_int_fore_subnational_err_F_K6_ETS_mean[,1])
      
ECP_F_alpha_0.95 = cbind(int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,1],
                         MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,1], 
                         MLFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,1],
                         horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95[,1],
                        
                         int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,1],
                         MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,1],
                         MLFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,1])

# M

ECP_M_alpha_0.8 = cbind(int_fore_subnational_err_M_EVR_ETS_mean[,1], 
                        MFTS_int_fore_subnational_err_M_EVR_ETS_mean[,1],
                        MLFTS_int_fore_subnational_err_M_EVR_ETS_mean[,1],
                        horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS[,1],
                        
                        int_fore_subnational_err_M_K6_ETS_mean[,1],
                        MFTS_int_fore_subnational_err_M_K6_ETS_mean[,1],
                        MLFTS_int_fore_subnational_err_M_K6_ETS_mean[,1])

ECP_M_alpha_0.95 = cbind(int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,1],
                         MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,1], 
                         MLFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,1],
                         horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95[,1],
                         
                         int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,1],
                         MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,1],
                         MLFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,1])

colnames(ECP_F_alpha_0.8) = colnames(ECP_F_alpha_0.95) = 
colnames(ECP_M_alpha_0.8) = colnames(ECP_M_alpha_0.95) = c("UFTS (EVR)", "MFTS (EVR)", "MLFTS (EVR)", "HDFPCA",
                                                           "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)")

######
# CPD
######

## F

CPD_F_alpha_0.8 = cbind(int_fore_subnational_err_F_EVR_ETS_mean[,2], 
                        MFTS_int_fore_subnational_err_F_EVR_ETS_mean[,2],
                        MLFTS_int_fore_subnational_err_F_EVR_ETS_mean[,2],
                        horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS[,2],
                        
                        int_fore_subnational_err_F_K6_ETS_mean[,2],
                        MFTS_int_fore_subnational_err_F_K6_ETS_mean[,2],
                        MLFTS_int_fore_subnational_err_F_K6_ETS_mean[,2])

CPD_F_alpha_0.95 = cbind(int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,2],
                         MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,2], 
                         MLFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,2],
                         horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95[,2],
                         
                         int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,2],
                         MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,2],
                         MLFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,2])

# M

CPD_M_alpha_0.8 = cbind(int_fore_subnational_err_M_EVR_ETS_mean[,2], 
                        MFTS_int_fore_subnational_err_M_EVR_ETS_mean[,2],
                        MLFTS_int_fore_subnational_err_M_EVR_ETS_mean[,2],
                        horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS[,2],
                        
                        int_fore_subnational_err_M_K6_ETS_mean[,2],
                        MFTS_int_fore_subnational_err_M_K6_ETS_mean[,2],
                        MLFTS_int_fore_subnational_err_M_K6_ETS_mean[,2])

CPD_M_alpha_0.95 = cbind(int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,2],
                         MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,2], 
                         MLFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,2],
                         horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95[,2],
                         
                         int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,2],
                         MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,2],
                         MLFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,2])

colnames(CPD_F_alpha_0.8) = colnames(CPD_F_alpha_0.95) = 
colnames(CPD_M_alpha_0.8) = colnames(CPD_M_alpha_0.95) = c("UFTS (EVR)", "MFTS (EVR)", "MLFTS (EVR)", "HDFPCA",
                                                           "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)")

#################
# interval score
#################

## F

score_F_alpha_0.8 = cbind(int_fore_subnational_err_F_EVR_ETS_mean[,3], 
                        MFTS_int_fore_subnational_err_F_EVR_ETS_mean[,3],
                        MLFTS_int_fore_subnational_err_F_EVR_ETS_mean[,3],
                        horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS[,3],
                        
                        int_fore_subnational_err_F_K6_ETS_mean[,3],
                        MFTS_int_fore_subnational_err_F_K6_ETS_mean[,3],
                        MLFTS_int_fore_subnational_err_F_K6_ETS_mean[,3])

score_F_alpha_0.95 = cbind(int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,3],
                         MFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,3], 
                         MLFTS_int_fore_subnational_err_F_EVR_ETS_alpha_0.95_mean[,3],
                         horizon_specific_hdfpca_int_fore_subnational_err_F_EVR_ETS_alpha_0.95[,3],
                         
                         int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,3],
                         MFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,3],
                         MLFTS_int_fore_subnational_err_F_K6_ETS_alpha_0.95_mean[,3])

# M

score_M_alpha_0.8 = cbind(int_fore_subnational_err_M_EVR_ETS_mean[,3], 
                        MFTS_int_fore_subnational_err_M_EVR_ETS_mean[,3],
                        MLFTS_int_fore_subnational_err_M_EVR_ETS_mean[,3],
                        horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS[,3],
                        
                        int_fore_subnational_err_M_K6_ETS_mean[,3],
                        MFTS_int_fore_subnational_err_M_K6_ETS_mean[,3],
                        MLFTS_int_fore_subnational_err_M_K6_ETS_mean[,3])

score_M_alpha_0.95 = cbind(int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,3],
                         MFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,3], 
                         MLFTS_int_fore_subnational_err_M_EVR_ETS_alpha_0.95_mean[,3],
                         horizon_specific_hdfpca_int_fore_subnational_err_M_EVR_ETS_alpha_0.95[,3],
                         
                         int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,3],
                         MFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,3],
                         MLFTS_int_fore_subnational_err_M_K6_ETS_alpha_0.95_mean[,3])

colnames(score_F_alpha_0.8) = colnames(score_F_alpha_0.95) = 
colnames(score_M_alpha_0.8) = colnames(score_M_alpha_0.95) = c("UFTS (EVR)", "MFTS (EVR)", "MLFTS (EVR)", "HDFPCA",
                                                               "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)")

###########
## summary
###########

# ECP

ECP_F_alpha_0.8_overall = rbind(ECP_F_alpha_0.8,  apply(ECP_F_alpha_0.8, 2, mean),  apply(ECP_F_alpha_0.8, 2, median))
ECP_M_alpha_0.8_overall = rbind(ECP_M_alpha_0.8,  apply(ECP_M_alpha_0.8, 2, mean),  apply(ECP_M_alpha_0.8, 2, median))

ECP_F_alpha_0.95_overall = rbind(ECP_F_alpha_0.95, apply(ECP_F_alpha_0.95, 2, mean), apply(ECP_F_alpha_0.95, 2, median))
ECP_M_alpha_0.95_overall = rbind(ECP_M_alpha_0.95, apply(ECP_M_alpha_0.95, 2, mean), apply(ECP_M_alpha_0.95, 2, median))

# CPD

CPD_F_alpha_0.8_overall = rbind(CPD_F_alpha_0.8,  apply(CPD_F_alpha_0.8, 2, mean),  apply(CPD_F_alpha_0.8, 2, median))
CPD_M_alpha_0.8_overall = rbind(CPD_M_alpha_0.8,  apply(CPD_M_alpha_0.8, 2, mean),  apply(CPD_M_alpha_0.8, 2, median))

CPD_F_alpha_0.95_overall = rbind(CPD_F_alpha_0.95, apply(CPD_F_alpha_0.95, 2, mean), apply(CPD_F_alpha_0.95, 2, median))
CPD_M_alpha_0.95_overall = rbind(CPD_M_alpha_0.95, apply(CPD_M_alpha_0.95, 2, mean), apply(CPD_M_alpha_0.95, 2, median))

# score

score_F_alpha_0.8_overall = rbind(score_F_alpha_0.8, apply(score_F_alpha_0.8, 2, mean), apply(score_F_alpha_0.8, 2, median))
score_M_alpha_0.8_overall = rbind(score_M_alpha_0.8, apply(score_M_alpha_0.8, 2, mean), apply(score_M_alpha_0.8, 2, median))

score_F_alpha_0.95_overall = rbind(score_F_alpha_0.95, apply(score_F_alpha_0.95, 2, mean), apply(score_F_alpha_0.95, 2, median))
score_M_alpha_0.95_overall = rbind(score_M_alpha_0.95, apply(score_M_alpha_0.95, 2, mean), apply(score_M_alpha_0.95, 2, median))

rownames(ECP_F_alpha_0.8_overall) = rownames(ECP_M_alpha_0.8_overall) = 
rownames(ECP_F_alpha_0.95_overall) = rownames(ECP_M_alpha_0.95_overall) = 
  
rownames(CPD_F_alpha_0.8_overall) = rownames(CPD_M_alpha_0.8_overall) = 
rownames(CPD_F_alpha_0.95_overall) = rownames(CPD_M_alpha_0.95_overall) = 
  
rownames(score_F_alpha_0.8_overall) = rownames(score_M_alpha_0.8_overall) = 
rownames(score_F_alpha_0.95_overall) = rownames(score_M_alpha_0.95_overall) = c(1:15, "Mean", "Median")

## alpha = 0.2

# EVR

xtable(cbind(CPD_F_alpha_0.8_overall[,1:4], CPD_M_alpha_0.8_overall[,1:4]), digits = 3)
xtable(cbind(score_F_alpha_0.8_overall[,1:4], score_M_alpha_0.8_overall[,1:4]), digits = 0)

# K = 6

xtable(cbind(CPD_F_alpha_0.8_overall[,4:6], CPD_M_alpha_0.8_overall[,4:6]), digits = 3)
xtable(cbind(score_F_alpha_0.8_overall[,4:6], score_M_alpha_0.8_overall[,4:6]), digits = 0)

#################
# ranking by CPD
# heatmaps
#################

# F

dim(int_fore_subnational_err_F_EVR_ETS)        # 15 3 47
dim(MFTS_int_fore_subnational_err_F_EVR_ETS)   # 15 3 47
dim(MLFTS_int_fore_subnational_err_F_EVR_ETS)  # 15 3 47
dim(FANOVA_int_fore_subnational_err_F_EVR_ETS) # 47 15 3
dim(hdfpca_int_fore_subnational_err_F_EVR_ETS) # 47 15 3

dim(int_fore_subnational_err_F_K6_ETS)       # 15 3 47
dim(MFTS_int_fore_subnational_err_F_K6_ETS)  # 15 3 47
dim(MLFTS_int_fore_subnational_err_F_K6_ETS) # 15 3 47
dim(FANOVA_int_fore_subnational_err_F_EVR_ETS_K6) # 47 15 3

CPD_fore_subnational_err_F_array = array(NA, dim = c(15, 47, 9), 
dimnames = list(1:15, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA", 
                              "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)", "FANOVA (K=6)")))
CPD_fore_subnational_err_F_array[,,1] = int_fore_subnational_err_F_EVR_ETS[,2,]
CPD_fore_subnational_err_F_array[,,2] = MFTS_int_fore_subnational_err_F_EVR_ETS[,2,]
CPD_fore_subnational_err_F_array[,,3] = MLFTS_int_fore_subnational_err_F_EVR_ETS[,2,]
CPD_fore_subnational_err_F_array[,,4] = t(FANOVA_int_fore_subnational_err_F_EVR_ETS[,,2])

CPD_fore_subnational_err_F_array[,,5] = t(hdfpca_int_fore_subnational_err_F_EVR_ETS[,,2])
CPD_fore_subnational_err_F_array[,,6] = int_fore_subnational_err_F_K6_ETS[,2,]
CPD_fore_subnational_err_F_array[,,7] = MFTS_int_fore_subnational_err_F_K6_ETS[,2,]
CPD_fore_subnational_err_F_array[,,8] = MLFTS_int_fore_subnational_err_F_K6_ETS[,2,]
CPD_fore_subnational_err_F_array[,,9] = t(FANOVA_int_fore_subnational_err_F_EVR_ETS_K6[,,2]) 


CPD_fore_subnational_err_F_ranking_mat = matrix(0, nrow = 15, ncol = 9)
colnames(CPD_fore_subnational_err_F_ranking_mat) = paste("M", 1:9, sep = "")
rownames(CPD_fore_subnational_err_F_ranking_mat) = 1:15
  
for(ih in 1:15)
{
    temp_rank = table(apply(CPD_fore_subnational_err_F_array[ih,,], 1, which.min))
    if(length(temp_rank) == 9)
    {
        CPD_fore_subnational_err_F_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, 9)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        CPD_fore_subnational_err_F_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

CPD_fore_subnational_err_F_ranking_mat = data.frame(Horizon = 1:15, CPD_fore_subnational_err_F_ranking_mat)

CPD_female = CPD_fore_subnational_err_F_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:15),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = ""))
)

# save figure

ggsave("Fig_5a.png")
ggplot(CPD_female, aes(Model, ordered(Horizon, levels = 15:1))) +
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

#######
# Male
#######

dim(int_fore_subnational_err_M_EVR_ETS) # 15 3 47
dim(MFTS_int_fore_subnational_err_M_EVR_ETS) # 15 3 47
dim(MLFTS_int_fore_subnational_err_M_EVR_ETS) # 15 3 47
dim(FANOVA_int_fore_subnational_err_M_EVR_ETS) # 47 15 3
dim(hdfpca_int_fore_subnational_err_M_EVR_ETS) # 47 15 3

dim(int_fore_subnational_err_M_K6_ETS) # 15 3 47
dim(MFTS_int_fore_subnational_err_M_K6_ETS) # 15 3 47
dim(MLFTS_int_fore_subnational_err_M_K6_ETS) # 15 3 47
dim(FANOVA_int_fore_subnational_err_M_EVR_ETS_K6) # 47 15 3

CPD_fore_subnational_err_M_array = array(NA, dim = c(15, 47, 9), 
        dimnames = list(1:15, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA", 
                                      "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)", "FANOVA (K=6)")))
CPD_fore_subnational_err_M_array[,,1] = int_fore_subnational_err_M_EVR_ETS[,2,]
CPD_fore_subnational_err_M_array[,,2] = MFTS_int_fore_subnational_err_M_EVR_ETS[,2,]
CPD_fore_subnational_err_M_array[,,3] = MLFTS_int_fore_subnational_err_M_EVR_ETS[,2,]
CPD_fore_subnational_err_M_array[,,4] = t(FANOVA_int_fore_subnational_err_M_EVR_ETS[,,2])
  
CPD_fore_subnational_err_M_array[,,5] = t(hdfpca_int_fore_subnational_err_M_EVR_ETS[,,2])
CPD_fore_subnational_err_M_array[,,6] = int_fore_subnational_err_M_K6_ETS[,2,]
CPD_fore_subnational_err_M_array[,,7] = MFTS_int_fore_subnational_err_M_K6_ETS[,2,]
CPD_fore_subnational_err_M_array[,,8] = MLFTS_int_fore_subnational_err_M_K6_ETS[,2,]
CPD_fore_subnational_err_M_array[,,9] = t(FANOVA_int_fore_subnational_err_M_EVR_ETS_K6[,,2])
  
CPD_fore_subnational_err_M_ranking_mat = matrix(0, nrow = 15, ncol = 9)
colnames(CPD_fore_subnational_err_M_ranking_mat) = paste("M", 1:9, sep = "")
rownames(CPD_fore_subnational_err_M_ranking_mat) = 1:15

for(ih in 1:15)
{
    temp_rank = table(apply(CPD_fore_subnational_err_M_array[ih,,],1, which.min))
    if(length(temp_rank) == 9)
    {
        CPD_fore_subnational_err_M_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, 9)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        CPD_fore_subnational_err_M_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

CPD_fore_subnational_err_M_ranking_mat = data.frame(Horizon = 1:15, CPD_fore_subnational_err_M_ranking_mat)

CPD_male = CPD_fore_subnational_err_M_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:15),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = "")))

ggsave("Fig_5b.png")
ggplot(CPD_male, aes(Model, ordered(Horizon, levels = 15:1))) +
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

############################
# ranking by interval score
# heatmaps
############################

# F

dim(int_fore_subnational_err_F_EVR_ETS) # 15 3 47
dim(MFTS_int_fore_subnational_err_F_EVR_ETS) # 15 3 47
dim(MLFTS_int_fore_subnational_err_F_EVR_ETS) # 15 3 47
dim(FANOVA_int_fore_subnational_err_F_EVR_ETS) # 47 15 3
dim(hdfpca_int_fore_subnational_err_F_EVR_ETS) # 47 15 3

dim(int_fore_subnational_err_F_K6_ETS) # 15 3 47
dim(MFTS_int_fore_subnational_err_F_K6_ETS) # 15 3 47
dim(MLFTS_int_fore_subnational_err_F_K6_ETS) # 15 3 47
dim(FANOVA_int_fore_subnational_err_F_EVR_ETS_K6) # 47 15 3

score_fore_subnational_err_F_array = array(NA, dim = c(15, 47, 9), 
    dimnames = list(1:15, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA", 
                                  "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)", "FANOVA (K=6)")))
score_fore_subnational_err_F_array[,,1] = int_fore_subnational_err_F_EVR_ETS[,3,]
score_fore_subnational_err_F_array[,,2] = MFTS_int_fore_subnational_err_F_EVR_ETS[,3,]
score_fore_subnational_err_F_array[,,3] = MLFTS_int_fore_subnational_err_F_EVR_ETS[,3,]
score_fore_subnational_err_F_array[,,4] = t(FANOVA_int_fore_subnational_err_F_EVR_ETS[,,3])

score_fore_subnational_err_F_array[,,5] = t(hdfpca_int_fore_subnational_err_F_EVR_ETS[,,3])
score_fore_subnational_err_F_array[,,6] = int_fore_subnational_err_F_K6_ETS[,3,]
score_fore_subnational_err_F_array[,,7] = MFTS_int_fore_subnational_err_F_K6_ETS[,3,]
score_fore_subnational_err_F_array[,,8] = MLFTS_int_fore_subnational_err_F_K6_ETS[,3,]
score_fore_subnational_err_F_array[,,9] = t(FANOVA_int_fore_subnational_err_F_EVR_ETS_K6[,,3])

score_fore_subnational_err_F_ranking_mat = matrix(0, nrow = 15, ncol = 9)
colnames(score_fore_subnational_err_F_ranking_mat) = paste("M", 1:9, sep = "")
rownames(score_fore_subnational_err_F_ranking_mat) = 1:15

for(ih in 1:15)
{
    temp_rank = table(apply(score_fore_subnational_err_F_array[ih,,], 1, which.min))
    if(length(temp_rank) == 9)
    {
        score_fore_subnational_err_F_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, 9)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        score_fore_subnational_err_F_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

score_fore_subnational_err_F_ranking_mat = data.frame(Horizon = 1:15, score_fore_subnational_err_F_ranking_mat)

score_female = score_fore_subnational_err_F_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:15),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = ""))
  )

# save figure

ggsave("Fig_5c.png")
ggplot(score_female, aes(Model, ordered(Horizon, levels = 15:1))) +
      geom_tile(aes(fill = count)) +
      geom_text(aes(label = count)) +
      scale_x_discrete(position = "top") +
      theme(legend.position = "none") + 
      scale_fill_gradient2(high="green",mid="white",low="red", 
                           na.value="yellow", midpoint = 5) + 
      ylab("Forecast horizon") +
      scale_x_discrete(labels = c("UFTS\n(EVR)", "MFTS\n(EVR)", "MLFTS\n(EVR)", "FANOVA\n(EVR)", "HDFPCA", 
                                  "UFTS\n(K=6)", "MFTS\n(K=6)", "MLFTS\n(K=6)", "FANOVA\n(K=6)")) + 
      theme_bw()
dev.off()

#######
# Male
#######

dim(int_fore_subnational_err_M_EVR_ETS) # 15 3 47
dim(MFTS_int_fore_subnational_err_M_EVR_ETS) # 15 3 47
dim(MLFTS_int_fore_subnational_err_M_EVR_ETS) # 15 3 47
dim(hdfpca_int_fore_subnational_err_M_EVR_ETS) # 47 15 3

dim(int_fore_subnational_err_M_K6_ETS) # 15 3 47
dim(MFTS_int_fore_subnational_err_M_K6_ETS) # 15 3 47
dim(MLFTS_int_fore_subnational_err_M_K6_ETS) # 15 3 47

score_fore_subnational_err_M_array = array(NA, dim = c(15, 47, 9), 
        dimnames = list(1:15, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA", 
                                      "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)", "FANOVA (K=6)")))
score_fore_subnational_err_M_array[,,1] = int_fore_subnational_err_M_EVR_ETS[,3,]
score_fore_subnational_err_M_array[,,2] = MFTS_int_fore_subnational_err_M_EVR_ETS[,3,]
score_fore_subnational_err_M_array[,,3] = MLFTS_int_fore_subnational_err_M_EVR_ETS[,3,]
score_fore_subnational_err_M_array[,,4] = t(FANOVA_int_fore_subnational_err_M_EVR_ETS[,,3])
  
score_fore_subnational_err_M_array[,,5] = t(hdfpca_int_fore_subnational_err_M_EVR_ETS[,,3])
score_fore_subnational_err_M_array[,,6] = int_fore_subnational_err_M_K6_ETS[,3,]
score_fore_subnational_err_M_array[,,7] = MFTS_int_fore_subnational_err_M_K6_ETS[,3,]
score_fore_subnational_err_M_array[,,8] = MLFTS_int_fore_subnational_err_M_K6_ETS[,3,]
score_fore_subnational_err_M_array[,,9] = t(FANOVA_int_fore_subnational_err_M_EVR_ETS_K6[,,3])
  
score_fore_subnational_err_M_ranking_mat = matrix(0, nrow = 15, ncol = 9)
colnames(score_fore_subnational_err_M_ranking_mat) = paste("M", 1:9, sep = "")
rownames(score_fore_subnational_err_M_ranking_mat) = 1:15

for(ih in 1:15)
{
    temp_rank = table(apply(score_fore_subnational_err_M_array[ih,,], 1, which.min))
    if(length(temp_rank) == 9)
    {
        score_fore_subnational_err_M_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, 9)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        score_fore_subnational_err_M_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

score_fore_subnational_err_M_ranking_mat = data.frame(Horizon = 1:15, score_fore_subnational_err_M_ranking_mat)

score_male = score_fore_subnational_err_M_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:15),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = "")))

# save figure

ggsave("Fig_5d.png")
ggplot(score_male, aes(Model, ordered(Horizon, levels = 15:1))) +
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

