###########################################################
# High-Dimensional Functional Principal Component Analysis
# transformation = "CDF"
###########################################################

# K = 6, L = 2

point_fore_subnational_err_HDFPCA_F_CDF = point_fore_subnational_err_HDFPCA_M_CDF = array(NA, dim = c(47, 16, 4))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        dum = hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], horizon = iw, 
                         first_order = 6, second_order = 2, transformation = "CDF")
        point_fore_subnational_err_HDFPCA_F_CDF[ij,iw,] = dum$err_F
        point_fore_subnational_err_HDFPCA_M_CDF[ij,iw,] = dum$err_M
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

###########################################
# For sensitivity analaysis, K = 2, L = 2
###########################################

point_fore_subnational_err_HDFPCA_F_CDF_K2 = point_fore_subnational_err_HDFPCA_M_CDF_K2 = array(NA, dim = c(47, 16, 4))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        dum = hdfpca_fun(fdata_F = female_prefecture_dx[[ij]], fdata_M = male_prefecture_dx[[ij]], horizon = iw, 
                         first_order = 2, second_order = 2, transformation = "CDF")
        point_fore_subnational_err_HDFPCA_F_CDF_K2[ij,iw,] = dum$err_F
        point_fore_subnational_err_HDFPCA_M_CDF_K2[ij,iw,] = dum$err_M
        rm(dum); rm(iw)
    }
    print(ij); rm(ij)
}

point_fore_subnational_err_HDFPCA_F_CDF_K2_mean = rbind(apply(point_fore_subnational_err_HDFPCA_F_CDF_K2, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_F_CDF_K2, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_F_CDF_K2, c(1,3), mean), 2, median)) # 0.0157 0.0047
rownames(point_fore_subnational_err_HDFPCA_F_CDF_K2_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_F_CDF_K2_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

point_fore_subnational_err_HDFPCA_M_CDF_K2_mean = rbind(apply(point_fore_subnational_err_HDFPCA_M_CDF_K2, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_M_CDF_K2, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_M_CDF_K2, c(1,3), mean), 2, median)) # 0.0141 0.0041
rownames(point_fore_subnational_err_HDFPCA_M_CDF_K2_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_M_CDF_K2_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

###########
## summary
###########

# F

point_fore_subnational_err_HDFPCA_F_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(1,3), mean), 2, median)) # 0.0157 0.0047
rownames(point_fore_subnational_err_HDFPCA_F_CDF_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_F_CDF_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")

horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(2, 3), mean), 
                                              colMeans(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(2, 3), mean)),
                                              apply(apply(point_fore_subnational_err_HDFPCA_F_CDF, c(2, 3), mean), 2, median))
rownames(horizon_point_fore_subnational_err_HDFPCA_F_CDF_mean) = c(1:16, "Mean", "Median")

# M

point_fore_subnational_err_HDFPCA_M_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(1,3), mean),
                                                     colMeans(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(1,3), mean)),
                                                     apply(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(1,3), mean), 2, median)) # 0.0141 0.0041
rownames(point_fore_subnational_err_HDFPCA_M_CDF_mean) = c(state, "Mean", "Median")
colnames(point_fore_subnational_err_HDFPCA_M_CDF_mean) = c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2")


horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean = rbind(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(2, 3), mean),
                                            colMeans(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(2, 3), mean)),
                                            apply(apply(point_fore_subnational_err_HDFPCA_M_CDF, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_HDFPCA_M_CDF_mean) = c(1:16, "Mean", "Median")

