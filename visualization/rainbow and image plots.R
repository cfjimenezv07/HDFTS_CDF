################
# rainbow plots
################

savefig("Fig_1a", width = 12, height = 10, type = "png", toplines = 0.8)
plot(fts(0:110, t(Japan_female_pop)), xlab = "Age", ylab = "Life-table death count", ylim = c(0, 5150))
dev.off()

savefig("Fig_1b", width = 12, height = 10, type = "png", toplines = 0.8)
plot(fts(0:110, t(Japan_male_pop)), xlab = "Age", ylab = "", ylim = c(0, 5150))
dev.off()

############################################
# image plots (Kullback-Leibler divergence)
############################################

# by year

female_KLdiv_year = male_KLdiv_year = total_KLdiv_year = matrix(NA, length(state), n_year)
for(iw in 1:length(state))
{
    for(ij in 1:n_year)
    {
        female_KLdiv_year[iw,ij] = mean(KLdiv(cbind((female_prefecture_dx[[iw]])[ij,], Japan_female_pop[ij,]))[2:3])
        male_KLdiv_year[iw,ij]   = mean(KLdiv(cbind((male_prefecture_dx[[iw]])[ij,],   Japan_male_pop[ij,]))[2:3])
        total_KLdiv_year[iw,ij]  = mean(KLdiv(cbind((total_prefecture_dx[[iw]])[ij,],  Japan_total_pop[ij,]))[2:3])
    }
}
rownames(female_KLdiv_year) = rownames(male_KLdiv_year) = rownames(total_KLdiv_year) = state
colnames(female_KLdiv_year) = colnames(male_KLdiv_year) = colnames(total_KLdiv_year) = year

savefig("Fig_2a", width = 12, height = 10, type = "png", toplines = 0.5)
image(years, 1:47, t(female_KLdiv_year[rev(1:47),]), xlab = "Year", ylab = "Prefecture", 
      zlim = c(0,0.16), main = "Female", yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig("Fig_2b", width = 12, height = 10, type = "png", toplines = 0.5)
image(years, 1:47, t(male_KLdiv_year[rev(1:47),]),   xlab = "Year", ylab = "", 
      zlim = c(0,0.16), main = "Male",yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig("Fig_2c", width = 12, height = 10, type = "png", toplines = 0.5)
image(years, 1:47, t(total_KLdiv_year[rev(1:47),]),  xlab = "Year", ylab = "", 
      zlim = c(0,0.16), main = "Total", yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

# by age

female_KLdiv_age = male_KLdiv_age = total_KLdiv_age = matrix(NA, length(state), n_age)
for(iw in 1:length(state))
{
    for(ij in 1:n_age)
    {
        female_KLdiv_age[iw,ij] = mean(KLdiv(cbind((female_prefecture_dx[[iw]])[,ij], Japan_female_pop[,ij]))[2:3])
        male_KLdiv_age[iw,ij]   = mean(KLdiv(cbind((male_prefecture_dx[[iw]])[,ij],   Japan_male_pop[,ij]))[2:3])
        total_KLdiv_age[iw,ij]  = mean(KLdiv(cbind((total_prefecture_dx[[iw]])[,ij],  Japan_total_pop[,ij]))[2:3])
    }
}
rownames(female_KLdiv_age) = rownames(male_KLdiv_age) = rownames(total_KLdiv_age) = state
colnames(female_KLdiv_age) = colnames(male_KLdiv_age) = colnames(total_KLdiv_age) = 0:110

savefig("Fig_2d", width = 12, height = 10, type = "png", toplines = 0.5)
image(0:110, 1:47, t(female_KLdiv_age[rev(1:47),]), xlab = "Age", ylab = "Prefecture", zlim = c(0, 1.26),
      yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig("Fig_2e", width = 12, height = 10, type = "png", toplines = 0.5)
image(0:110, 1:47, t(male_KLdiv_age[rev(1:47),]),   xlab = "Age", ylab = "", zlim = c(0, 1.26), yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig("Fig_2f", width = 12, height = 10, type = "png", toplines = 0.5)
image(0:110, 1:47, t(total_KLdiv_age[rev(1:47),]),  xlab = "Age", ylab = "", zlim = c(0, 1.26), yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

