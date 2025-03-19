####################################
# Multilevel functional time series
####################################

point_fore_subnational_err_MLFTS_F_EVR_ETS = point_fore_subnational_err_MLFTS_M_EVR_ETS =
point_fore_subnational_err_MLFTS_F_K6_ETS  = point_fore_subnational_err_MLFTS_M_K6_ETS  = array(NA, dim = c(47, 16, 4))
for(ij in 1:47)
{
    for(iw in 1:16)
    {
        ## EVR
        
        dum = point_fore_national_cdf_MLFTS(fdata_F = female_prefecture_dx[[ij]],
                                            fdata_M = male_prefecture_dx[[ij]],
                                            fdata_common = NULL, fore_method = "CDF",
                                            horizon = iw, way_ncomp = "EVR")
        point_fore_subnational_err_MLFTS_F_EVR_ETS[ij,iw,] = dum$err_F
        point_fore_subnational_err_MLFTS_M_EVR_ETS[ij,iw,] = dum$err_M
        rm(dum)
        
        ## K = 6
        
        dum = point_fore_national_cdf_MLFTS(fdata_F = female_prefecture_dx[[ij]],
                                           fdata_M = male_prefecture_dx[[ij]],
                                           fdata_common = NULL, fore_method = "CDF",
                                           horizon = iw, way_ncomp = "provide")
        point_fore_subnational_err_MLFTS_F_K6_ETS[ij,iw,] = dum$err_F
        point_fore_subnational_err_MLFTS_M_K6_ETS[ij,iw,] = dum$err_M
        rm(dum)
    }
    print(ij); rm(ij)
}
dimnames(point_fore_subnational_err_MLFTS_F_EVR_ETS) = dimnames(point_fore_subnational_err_MLFTS_M_EVR_ETS) =
dimnames(point_fore_subnational_err_MLFTS_F_K6_ETS) = dimnames(point_fore_subnational_err_MLFTS_M_K6_ETS) = list(state, 1:16, c("KLD", "JSD (geo)", "Wasserstein L1", "Wasserstein L2"))

#######
## ETS
#######

## female

# average across horizons

point_fore_subnational_err_MLFTS_F_EVR_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(1, 3), mean),
                                                             colMeans(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(1, 3), mean)),
                                                             apply(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(1, 3), mean), 2, median)), 4)

point_fore_subnational_err_MLFTS_F_K6_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(1, 3), mean),
                                                            colMeans(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(1, 3), mean)),
                                                            apply(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(1, 3), mean), 2, median)), 4)
rownames(point_fore_subnational_err_MLFTS_F_EVR_ETS_mean) = rownames(point_fore_subnational_err_MLFTS_F_K6_ETS_mean) = c(state, "Mean", "Median") 

# average across prefectures

horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_F_EVR_ETS, c(2, 3), mean), 2, median))

horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_F_K6_ETS, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_MLFTS_F_EVR_ETS_mean) = rownames(horizon_point_fore_subnational_err_MLFTS_F_K6_ETS_mean) = c(1:16, "Mean", "Median")

## male

# average across horizons

point_fore_subnational_err_MLFTS_M_EVR_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(1, 3), mean),
                                                             colMeans(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(1, 3), mean)),
                                                             apply(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(1, 3), mean), 2, median)), 4)

point_fore_subnational_err_MLFTS_M_K6_ETS_mean = round(rbind(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(1, 3), mean),
                                                            colMeans(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(1, 3), mean)),
                                                            apply(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(1, 3), mean), 2, median)), 4)
rownames(point_fore_subnational_err_MLFTS_M_EVR_ETS_mean) = rownames(point_fore_subnational_err_MLFTS_M_K6_ETS_mean) = c(state, "Mean", "Median")

# average across prefectures

horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_M_EVR_ETS, c(2, 3), mean), 2, median))

horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean = rbind(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(2, 3), mean),
                                    colMeans(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(2, 3), mean)),
                                    apply(apply(point_fore_subnational_err_MLFTS_M_K6_ETS, c(2, 3), mean), 2, median))

rownames(horizon_point_fore_subnational_err_MLFTS_M_EVR_ETS_mean) = rownames(horizon_point_fore_subnational_err_MLFTS_M_K6_ETS_mean) = c(1:16, "Mean", "Median")


