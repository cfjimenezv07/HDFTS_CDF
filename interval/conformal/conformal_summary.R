######
# CPD
######

# F

conformal_F_mean = cbind(int_fore_subnational_err_F_EVR_conformal_mean[,2],
                         MFTS_int_fore_subnational_err_F_EVR_conformal_mean[,2],
                         MLFTS_int_fore_subnational_err_F_EVR_conformal_mean[,2],
                         hdfpca_int_fore_subnational_err_F_EVR_conformal_mean[,2],
                         
                         int_fore_subnational_err_F_K6_conformal_mean[,2],
                         MFTS_int_fore_subnational_err_F_K6_conformal_mean[,2],
                         MLFTS_int_fore_subnational_err_F_K6_conformal_mean[,2])

conformal_F_mean_overall = rbind(conformal_F_mean, colMeans(conformal_F_mean), apply(conformal_F_mean, 2, median))
rownames(conformal_F_mean_overall) = c(1:15, "Mean", "Median")
colnames(conformal_F_mean_overall) = c("UFTS", "MFTS", "MLFTS", "HDFPCA", "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)")

# M

conformal_M_mean = cbind(int_fore_subnational_err_M_EVR_conformal_mean[,2],
                         MFTS_int_fore_subnational_err_M_EVR_conformal_mean[,2],
                         MLFTS_int_fore_subnational_err_M_EVR_conformal_mean[,2],
                         hdfpca_int_fore_subnational_err_M_EVR_conformal_mean[,2],
                         
                         int_fore_subnational_err_M_K6_conformal_mean[,2],
                         MFTS_int_fore_subnational_err_M_K6_conformal_mean[,2],
                         MLFTS_int_fore_subnational_err_M_K6_conformal_mean[,2])

conformal_M_mean_overall = rbind(conformal_M_mean, colMeans(conformal_M_mean), apply(conformal_M_mean, 2, median))
rownames(conformal_M_mean_overall) = c(1:15, "Mean", "Median")
colnames(conformal_M_mean_overall) = c("UFTS", "MFTS", "MLFTS", "HDFPCA", "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)")

#################
# Interval score
#################

# F

score_conformal_F_mean = cbind(int_fore_subnational_err_F_EVR_conformal_mean[,3],
                               MFTS_int_fore_subnational_err_F_EVR_conformal_mean[,3],
                               MLFTS_int_fore_subnational_err_F_EVR_conformal_mean[,3],
                               hdfpca_int_fore_subnational_err_F_EVR_conformal_mean[,3],
                               
                               int_fore_subnational_err_F_K6_conformal_mean[,3],
                               MFTS_int_fore_subnational_err_F_K6_conformal_mean[,3],
                               MLFTS_int_fore_subnational_err_F_K6_conformal_mean[,3])

score_conformal_F_mean_overall = rbind(score_conformal_F_mean, colMeans(score_conformal_F_mean), apply(score_conformal_F_mean, 2, median))
rownames(score_conformal_F_mean_overall) = c(1:15, "Mean", "Median")
colnames(score_conformal_F_mean_overall) = c("UFTS", "MFTS", "MLFTS", "HDFPCA", "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)")

# M

score_conformal_M_mean = cbind(int_fore_subnational_err_M_EVR_conformal_mean[,3],
                               MFTS_int_fore_subnational_err_M_EVR_conformal_mean[,3],
                               MLFTS_int_fore_subnational_err_M_EVR_conformal_mean[,3],
                               hdfpca_int_fore_subnational_err_M_EVR_conformal_mean[,3],
                               
                               int_fore_subnational_err_M_K6_conformal_mean[,3],
                               MFTS_int_fore_subnational_err_M_K6_conformal_mean[,3],
                               MLFTS_int_fore_subnational_err_M_K6_conformal_mean[,3])

score_conformal_M_mean_overall = rbind(score_conformal_M_mean, colMeans(score_conformal_M_mean), apply(score_conformal_M_mean, 2, median))
rownames(score_conformal_M_mean_overall) = c(1:15, "Mean", "Median")
colnames(score_conformal_M_mean_overall) = c("UFTS", "MFTS", "MLFTS", "HDFPCA", "UFTS(K=6)", "MFTS(K=6)", "MLFTS(K=6)")

# summary output

xtable(conformal_F_mean_overall, digits = 3)
xtable(conformal_M_mean_overall, digits = 3)
xtable(score_conformal_F_mean_overall, digits = 0)
xtable(score_conformal_M_mean_overall, digits = 0)

