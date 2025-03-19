#####################
# KLD, JSD
# CDF transformation
#####################

K6_res_KLD_female_male = cbind(horizon_point_fore_subnational_err_F_K6_ETS_mean[,1],
                               horizon_point_fore_subnational_err_MFTS_F_K6_ETS_mean[,1],
                               horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean[,1],
                                
                               horizon_point_fore_subnational_err_M_K6_ETS_mean[,1],
                               horizon_point_fore_subnational_err_MFTS_M_K6_ETS_mean[,1],
                               horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean[,1])

K6_res_JSD_female_male = cbind(horizon_point_fore_subnational_err_F_K6_ETS_mean[,2],
                               horizon_point_fore_subnational_err_MFTS_F_K6_ETS_mean[,2],
                               horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean[,2],
                               
                               horizon_point_fore_subnational_err_M_K6_ETS_mean[,2],
                               horizon_point_fore_subnational_err_MFTS_M_K6_ETS_mean[,2],
                               horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean[,2])

colnames(K6_res_KLD_female_male) = colnames(K6_res_JSD_female_male) = rep(c("UFTS", "MFTS", "MLFTS"), 2)

#########
# female
#########

# average across horizons

res_KLD_female = cbind(point_fore_subnational_err_F_EVR_ETS_mean[,1], 
                       point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,1],
                       point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,1],
                       point_fore_subnational_err_HDFPCA_F_CDF_mean[,1])

res_JSD_female = cbind(point_fore_subnational_err_F_EVR_ETS_mean[,2], 
                       point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,2],
                       point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,2],
                       point_fore_subnational_err_HDFPCA_F_CDF_mean[,2])

res_RMSDE1_female = cbind(point_fore_subnational_err_F_EVR_ETS_mean[,3], 
                          point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,3],
                          point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,3],
                          point_fore_subnational_err_HDFPCA_F_CDF_mean[,3])

res_RMSDE2_female = cbind(point_fore_subnational_err_F_EVR_ETS_mean[,4], 
                          point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,4],
                          point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,4],
                          point_fore_subnational_err_HDFPCA_F_CDF_mean[,4])
colnames(res_KLD_female) = colnames(res_JSD_female) = 
colnames(res_RMSDE1_female) = colnames(res_RMSDE2_female) = c("UFTS", "MFTS", "MLFTS", "HDFPCA")

# figures

savefig("res_KLD_female", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_KLD_female[1:47,], main = "KLD")
dev.off()

savefig("res_JSD_female", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_JSD_female[1:47,], main = "JSD")
dev.off()

savefig("res_RMSDE1_female", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE1_female[1:47,], main = expression(RMSDE[1]))
dev.off()

savefig("res_RMSDE2_female", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE2_female[1:47,], main = expression(RMSDE[2]))
dev.off()

# average across prefectures (UFTS, MFTS, MLFTS, HDFPCA, FANOVA)

horizon_res_KLD_female = cbind(horizon_point_fore_subnational_err_F_EVR_ETS_mean[,1], 
                               horizon_point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,1],
                               horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,1],
                               horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean[,1],
                               errors_FANOVA_EVR_F_overall[,1])

horizon_res_JSD_female = cbind(horizon_point_fore_subnational_err_F_EVR_ETS_mean[,2], 
                               horizon_point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,2],
                               horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,2],
                               horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean[,2],
                               errors_FANOVA_EVR_F_overall[,2])

horizon_res_RMSDE1_female = cbind(horizon_point_fore_subnational_err_F_EVR_ETS_mean[,3], 
                          horizon_point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,3],
                          horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,3],
                          horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean[,3])

horizon_res_RMSDE2_female = cbind(horizon_point_fore_subnational_err_F_EVR_ETS_mean[,4], 
                          horizon_point_fore_subnational_err_MFTS_F_EVR_ETS_mean[,4],
                          horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean[,4],
                          horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean[,4])
colnames(horizon_res_KLD_female) = colnames(horizon_res_JSD_female) = 
colnames(horizon_res_RMSDE1_female) = colnames(horizon_res_RMSDE2_female) = c("UFTS", "MFTS", "MLFTS", "HDFPCA")

#######
# male
#######

# average across horizons

res_KLD_male = cbind(point_fore_subnational_err_M_EVR_ETS_mean[,1],
                     point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,1],
                     point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,1],
                     point_fore_subnational_err_HDFPCA_M_CDF_mean[,1],
                     errors_FANOVA_EVR_F_overall[,1])

res_JSD_male = cbind(point_fore_subnational_err_M_EVR_ETS_mean[,2],
                     point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,2],
                     point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,2],
                     point_fore_subnational_err_HDFPCA_M_CDF_mean[,2])

res_RMSDE1_male = cbind(point_fore_subnational_err_M_EVR_ETS_mean[,3],
                        point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,3],
                        point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,3],
                        point_fore_subnational_err_HDFPCA_M_CDF_mean[,3])

res_RMSDE2_male = cbind(point_fore_subnational_err_M_EVR_ETS_mean[,4],
                        point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,4],
                        point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,4],
                        point_fore_subnational_err_HDFPCA_M_CDF_mean[,4])
colnames(res_KLD_male) = colnames(res_JSD_male) = 
colnames(res_RMSDE1_male) = colnames(res_RMSDE2_male) = c("UFTS", "MFTS", "MLFTS", "HDFPCA")

# figures

savefig("res_KLD_male", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_KLD_male[1:47,], main = "KLD")
dev.off()

savefig("res_JSD_male", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_JSD_male[1:47,], main = "JSD")
dev.off()

savefig("res_RMSDE1_male", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE1_male[1:47,], main = expression(RMSDE[1]))
dev.off()

savefig("res_RMSDE2_male", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE2_male[1:47,], main = expression(RMSDE[2]))
dev.off()

# average across prefectures (UFTS, MFTS, MLFTS, HDFPCA, FANOVA)

horizon_res_KLD_male = cbind(horizon_point_fore_subnational_err_M_EVR_ETS_mean[,1],
                             horizon_point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,1],
                             horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,1],
                             horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean[,1],
                             errors_FANOVA_EVR_M_overall[,1])

horizon_res_JSD_male = cbind(horizon_point_fore_subnational_err_M_EVR_ETS_mean[,2],
                             horizon_point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,2],
                             horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,2],
                             horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean[,2],
                             errors_FANOVA_EVR_M_overall[,2])

horizon_res_RMSDE1_male = cbind(horizon_point_fore_subnational_err_M_EVR_ETS_mean[,3],
                                horizon_point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,3],
                                horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,3],
                                horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean[,3])

horizon_res_RMSDE2_male = cbind(horizon_point_fore_subnational_err_M_EVR_ETS_mean[,4],
                                horizon_point_fore_subnational_err_MFTS_M_EVR_ETS_mean[,4],
                                horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean[,4],
                                horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean[,4])
colnames(horizon_res_KLD_male) = colnames(horizon_res_JSD_male) = 
colnames(horizon_res_RMSDE1_male) = colnames(horizon_res_RMSDE2_male) = c("UFTS", "MFTS", "MLFTS", "HDFPCA")

#################
# output results
#################

require(xtable)

xtable(cbind(horizon_res_KLD_female, horizon_res_KLD_male), digits = 3)
xtable(cbind(horizon_res_JSD_female, horizon_res_JSD_male), digits = 3)


###########
# heatmaps
###########

n_methods = 9
KLD_fore_subnational_err_F_array = array(NA, dim = c(16, 47, n_methods), 
      dimnames = list(1:16, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA",
                                    "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)", "FANOVA(K=6)")))
KLD_fore_subnational_err_F_array[,,1] = t(point_fore_subnational_err_F_EVR_ETS[,,1])
KLD_fore_subnational_err_F_array[,,2] = t(point_fore_subnational_err_MFTS_F_EVR_ETS[,,1])
KLD_fore_subnational_err_F_array[,,3] = t(point_fore_subnational_err_MLFTS_F_EVR_ETS[,,1])
KLD_fore_subnational_err_F_array[,,4] = FANOVA_KLD_F_EVR
KLD_fore_subnational_err_F_array[,,5] = t(point_fore_subnational_err_HDFPCA_F_CDF[,,1])

KLD_fore_subnational_err_F_array[,,6] = t(point_fore_subnational_err_F_K6_ETS[,,1])
KLD_fore_subnational_err_F_array[,,7] = t(point_fore_subnational_err_MFTS_F_K6_ETS[,,1])
KLD_fore_subnational_err_F_array[,,8] = t(point_fore_subnational_err_MLFTS_F_K6_ETS[,,1])
KLD_fore_subnational_err_F_array[,,9] = FANOVA_KLD_F_K6

KLD_fore_subnational_err_F_ranking_mat = matrix(0, nrow = 16, ncol = n_methods)
colnames(KLD_fore_subnational_err_F_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(KLD_fore_subnational_err_F_ranking_mat) = 1:16

for(ih in 1:16)
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

KLD_fore_subnational_err_F_ranking_mat = data.frame(Horizon = 1:16, KLD_fore_subnational_err_F_ranking_mat)

KLD_female = KLD_fore_subnational_err_F_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:16),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = ""))
  )

ggsave("Fig_4a.png")
ggplot(KLD_female, aes(Model, ordered(Horizon, levels = 16:1))) +
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

KLD_fore_subnational_err_M_array = array(NA, dim = c(16, 47, 9), 
    dimnames = list(1:16, 1:47, c("UFTS", "MFTS", "MLFTS", "FANOVA", "HDFPCA", "UFTS (K=6)", "MFTS (K=6)", "MLFTS (K=6)", "FANOVA (K=6)")))
KLD_fore_subnational_err_M_array[,,1] = t(point_fore_subnational_err_M_EVR_ETS[,,1])
KLD_fore_subnational_err_M_array[,,2] = t(point_fore_subnational_err_MFTS_M_EVR_ETS[,,1])
KLD_fore_subnational_err_M_array[,,3] = t(point_fore_subnational_err_MLFTS_M_EVR_ETS[,,1])
KLD_fore_subnational_err_M_array[,,4] = FANOVA_KLD_M_EVR
KLD_fore_subnational_err_M_array[,,5] = t(point_fore_subnational_err_HDFPCA_M_CDF[,,1])

KLD_fore_subnational_err_M_array[,,6] = t(point_fore_subnational_err_M_K6_ETS[,,1])
KLD_fore_subnational_err_M_array[,,7] = t(point_fore_subnational_err_MFTS_M_K6_ETS[,,1])
KLD_fore_subnational_err_M_array[,,8] = t(point_fore_subnational_err_MLFTS_M_K6_ETS[,,1])
KLD_fore_subnational_err_M_array[,,9] = FANOVA_KLD_M_K6

KLD_fore_subnational_err_M_ranking_mat = matrix(0, nrow = 16, ncol = 9)
colnames(KLD_fore_subnational_err_M_ranking_mat) = paste("M", 1:9, sep = "")
rownames(KLD_fore_subnational_err_M_ranking_mat) = 1:16

for(ih in 1:16)
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

KLD_fore_subnational_err_M_ranking_mat = data.frame(Horizon = 1:16, KLD_fore_subnational_err_M_ranking_mat)

KLD_male = KLD_fore_subnational_err_M_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:16),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:9, sep = ""))
  )

ggsave("Fig_4b.png")
ggplot(KLD_male, aes(Model, ordered(Horizon, levels = 16:1))) +
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

######
# clr
######

#########
# female
#########

res_KLD_female_clr = cbind(point_fore_subnational_err_F_EVR_ETS_clr_mean[,1], 
                           point_fore_subnational_err_F_EVR_MFTS_clr_mean[,1],
                           point_fore_subnational_err_F_EVR_MLFTS_clr_mean[,1],
                           point_fore_subnational_err_HDFPCA_F_clr_mean[,1])

res_JSD_female_clr = cbind(point_fore_subnational_err_F_EVR_ETS_clr_mean[,2], 
                           point_fore_subnational_err_F_EVR_MFTS_clr_mean[,2],
                           point_fore_subnational_err_F_EVR_MLFTS_clr_mean[,2],
                           point_fore_subnational_err_HDFPCA_F_clr_mean[,2])

res_RMSDE1_female_clr = cbind(point_fore_subnational_err_F_EVR_ETS_clr_mean[,3], 
                              point_fore_subnational_err_F_EVR_MFTS_clr_mean[,3],
                              point_fore_subnational_err_F_EVR_MLFTS_clr_mean[,3],
                              point_fore_subnational_err_HDFPCA_F_clr_mean[,3])

res_RMSDE2_female_clr = cbind(point_fore_subnational_err_F_EVR_ETS_clr_mean[,4], 
                              point_fore_subnational_err_F_EVR_MFTS_clr_mean[,4],
                              point_fore_subnational_err_F_EVR_MLFTS_clr_mean[,4],
                              point_fore_subnational_err_HDFPCA_F_clr_mean[,4])
colnames(res_KLD_female_clr) = colnames(res_JSD_female_clr) = 
colnames(res_RMSDE1_female_clr) = colnames(res_RMSDE2_female_clr) = c("UFTS", "MFTS", "MLFTS", "HDFPCA")

# figures

savefig("res_KLD_female_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_KLD_female_clr[1:47,], main = "KLD")
dev.off()

savefig("res_JSD_female_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_JSD_female_clr[1:47,], main = "JSD")
dev.off()

savefig("res_RMSDE1_female_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE1_female_clr[1:47,], main = expression(RMSDE[1]))
dev.off()

savefig("res_RMSDE2_female_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE2_female_clr[1:47,], main = expression(RMSDE[2]))
dev.off()

#######
# male
#######

res_KLD_male_clr = cbind(point_fore_subnational_err_M_EVR_ETS_clr_mean[,1], 
                         point_fore_subnational_err_M_EVR_MFTS_clr_mean[,1],
                         point_fore_subnational_err_M_EVR_MLFTS_clr_mean[,1],
                         point_fore_subnational_err_HDFPCA_M_clr_mean[,1])

res_JSD_male_clr = cbind(point_fore_subnational_err_M_EVR_ETS_clr_mean[,2], 
                         point_fore_subnational_err_M_EVR_MFTS_clr_mean[,2],
                         point_fore_subnational_err_M_EVR_MLFTS_clr_mean[,2],
                         point_fore_subnational_err_HDFPCA_M_clr_mean[,2])

res_RMSDE1_male_clr = cbind(point_fore_subnational_err_M_EVR_ETS_clr_mean[,3], 
                            point_fore_subnational_err_M_EVR_MFTS_clr_mean[,3],
                            point_fore_subnational_err_M_EVR_MLFTS_clr_mean[,3],
                            point_fore_subnational_err_HDFPCA_M_clr_mean[,3])

res_RMSDE2_male_clr = cbind(point_fore_subnational_err_M_EVR_ETS_clr_mean[,4], 
                            point_fore_subnational_err_M_EVR_MFTS_clr_mean[,4],
                            point_fore_subnational_err_M_EVR_MLFTS_clr_mean[,4],
                            point_fore_subnational_err_HDFPCA_M_clr_mean[,4])
colnames(res_KLD_male_clr) = colnames(res_JSD_male_clr) = 
colnames(res_RMSDE1_male_clr) = colnames(res_RMSDE2_male_clr) = c("UFTS", "MFTS", "MLFTS", "HDFPCA")

# figures

savefig("res_KLD_male_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_KLD_male_clr[1:47,], main = "KLD")
dev.off()

savefig("res_JSD_male_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_JSD_male_clr[1:47,], main = "JSD")
dev.off()

savefig("res_RMSDE1_male_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE1_male_clr[1:47,], main = expression(RMSDE[1]))
dev.off()

savefig("res_RMSDE2_male_clr", width = 12, height = 10, toplines = 0.8, type = "png")
vioplot(res_RMSDE2_male_clr[1:47,], main = expression(RMSDE[2]))
dev.off()

res_female_summary = cbind(cbind(res_KLD_female[48,], res_KLD_female_clr[48,]),
                           cbind(res_JSD_female[48,], res_JSD_female_clr[48,]),
                           cbind(res_RMSDE1_female[48,], res_RMSDE1_female_clr[48,]),
                           cbind(res_RMSDE2_female[48,], res_RMSDE2_female_clr[48,]))

res_male_summary = cbind(cbind(res_KLD_male[48,], res_KLD_male_clr[48,]), 
                         cbind(res_JSD_male[48,], res_JSD_male_clr[48,]),
                         cbind(res_RMSDE1_male[48,], res_RMSDE1_male_clr[48,]),
                         cbind(res_RMSDE2_male[48,], res_RMSDE2_male_clr[48,]))
colnames(res_female_summary) = colnames(res_male_summary) = rep(c("CDF", "clr"), 4)

require(xtable)
xtable(rbind(res_female_summary, res_male_summary), digits = 4)

