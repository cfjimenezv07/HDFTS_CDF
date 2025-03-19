####################################
# Univariate functional time series
####################################

# point forecast error for h = 1,...,16

point_fore_subnational_err_F_EVR_ETS   = point_fore_subnational_err_M_EVR_ETS =
point_fore_subnational_err_F_K6_ETS    = point_fore_subnational_err_M_K6_ETS = array(NA, dim = c(47, 16, 4))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## EVR
        
        point_fore_subnational_err_F_EVR_ETS[ij,iw,] = point_fore_national_cdf(fdata = female_prefecture_dx[[ij]], method_ncomp = "EVR",
                                                                               horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                                               variable_interest = "point")$err
        point_fore_subnational_err_M_EVR_ETS[ij,iw,] = point_fore_national_cdf(fdata = male_prefecture_dx[[ij]], method_ncomp = "EVR",
                                                                               horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                                               variable_interest = "point")$err
        ## K = 6
        
        point_fore_subnational_err_F_K6_ETS[ij,iw,] = point_fore_national_cdf(fdata = female_prefecture_dx[[ij]], method_ncomp = "provide",
                                                                              horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                                              variable_interest = "point")$err
        point_fore_subnational_err_M_K6_ETS[ij,iw,] = point_fore_national_cdf(fdata = male_prefecture_dx[[ij]], method_ncomp = "provide",
                                                                              horizon = iw, fore_method = "CDF", uni_fore_method = "ets",
                                                                              variable_interest = "point")$err
        rm(iw)
    }
    print(ij); rm(ij)
}
dimnames(point_fore_subnational_err_F_EVR_ETS)   = dimnames(point_fore_subnational_err_M_EVR_ETS) =
dimnames(point_fore_subnational_err_F_K6_ETS)    = dimnames(point_fore_subnational_err_M_K6_ETS) = list(state, 1:16, c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2"))

#######
## ETS
#######

## female

# average across horizons

point_fore_subnational_err_F_EVR_ETS_mean = round(rbind(apply(point_fore_subnational_err_F_EVR_ETS, c(1, 3), mean),
                                                        colMeans(apply(point_fore_subnational_err_F_EVR_ETS, c(1, 3), mean)),
                                                        apply(apply(point_fore_subnational_err_F_EVR_ETS, c(1, 3), mean),2,median)), 4)

point_fore_subnational_err_F_K6_ETS_mean = round(rbind(apply(point_fore_subnational_err_F_K6_ETS, c(1, 3), mean),
                                                       colMeans(apply(point_fore_subnational_err_F_K6_ETS, c(1, 3), mean)),
                                                       apply(apply(point_fore_subnational_err_F_K6_ETS, c(1, 3), mean), 2, median)), 4)

rownames(point_fore_subnational_err_F_EVR_ETS_mean) = rownames(point_fore_subnational_err_F_K6_ETS_mean) = c(state, "Mean", "Median")

# average across prefectures

horizon_point_fore_subnational_err_F_EVR_ETS_mean = rbind(apply(point_fore_subnational_err_F_EVR_ETS, c(2, 3), mean),
                                        colMeans(apply(point_fore_subnational_err_F_EVR_ETS, c(2, 3), mean)),
                                        apply(apply(point_fore_subnational_err_F_EVR_ETS, c(2, 3), mean), 2, median))

horizon_point_fore_subnational_err_F_K6_ETS_mean = rbind(apply(point_fore_subnational_err_F_K6_ETS, c(2, 3), mean),
                                        colMeans(apply(point_fore_subnational_err_F_K6_ETS, c(2, 3), mean)),
                                        apply(apply(point_fore_subnational_err_F_K6_ETS, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_F_EVR_ETS_mean) = rownames(horizon_point_fore_subnational_err_F_K6_ETS_mean) = c(1:16, "Mean", "Median")

## male

# average across horizons

point_fore_subnational_err_M_EVR_ETS_mean = round(rbind(apply(point_fore_subnational_err_M_EVR_ETS, c(1, 3), mean),
                                                        colMeans(apply(point_fore_subnational_err_M_EVR_ETS, c(1, 3), mean)),
                                                        apply(apply(point_fore_subnational_err_M_EVR_ETS, c(1, 3), mean), 2, median)), 4)

point_fore_subnational_err_M_K6_ETS_mean = round(rbind(apply(point_fore_subnational_err_M_K6_ETS, c(1, 3), mean),
                                                       colMeans(apply(point_fore_subnational_err_M_K6_ETS, c(1, 3), mean)),
                                                       apply(apply(point_fore_subnational_err_M_K6_ETS, c(1, 3), mean), 2, median)), 4)

rownames(point_fore_subnational_err_M_EVR_ETS_mean) = rownames(point_fore_subnational_err_M_K6_ETS_mean) = c(state, "Mean", "Median")

# average across prefectures

horizon_point_fore_subnational_err_M_EVR_ETS_mean = rbind(apply(point_fore_subnational_err_M_EVR_ETS, c(2, 3), mean), 
                                                colMeans(apply(point_fore_subnational_err_M_EVR_ETS, c(2, 3), mean)),
                                                apply(apply(point_fore_subnational_err_M_EVR_ETS, c(2, 3), mean), 2, median))

horizon_point_fore_subnational_err_M_K6_ETS_mean = rbind(apply(point_fore_subnational_err_M_K6_ETS, c(2, 3), mean),
                                                colMeans(apply(point_fore_subnational_err_M_K6_ETS, c(2, 3), mean)),
                                                apply(apply(point_fore_subnational_err_M_K6_ETS, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_M_EVR_ETS_mean) = rownames(horizon_point_fore_subnational_err_M_K6_ETS_mean) = c(1:16, "Mean", "Median")

