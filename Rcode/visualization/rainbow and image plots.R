
# Where heatmaps will be saved
dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/Plots/"
dir.r1 <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/HDFTS_CDF/HDFTS_CDF/code_submission/Results_new/FANOVA/"


Japan_female_pop <- readRDS(paste0(dir.r1, "Japan_female_pop.rds"))
Japan_male_pop   <- readRDS(paste0(dir.r1, "Japan_male_pop.rds"))


savefig <- function (filename, height=10, width = (1 + sqrt(5))/2*height, type=c("eps","pdf","jpg","png"), pointsize = 10, family = "Helvetica", sublines = 0, toplines = 0, leftlines = 0, res=300) 
{
  type <- match.arg(type)
  filename <- paste(filename, ".", type, sep = "")
  if(type=="eps")
  {
    postscript(file = filename, horizontal = FALSE, 
               width = width/2.54, height = height/2.54, pointsize = pointsize, 
               family = family, onefile = FALSE, print.it = FALSE)
  }
  else if(type=="pdf")
  {
    pdf(file = filename, width=width/2.54, height=height/2.54, pointsize=pointsize,
        family=family, onefile=TRUE)
  }
  else if(type=="jpg")
  {
    jpeg(filename=filename, width=width, height=height, res=res,quality=100, units="cm")#, pointsize=pointsize*50)
  }
  else if(type=="png")
  {
    png(filename=filename, width=width, height=height, res=res, units="cm")#, pointsize=pointsize*50)
  }
  else
    stop("Unknown file type")
  par(mgp = c(2.2, 0.45, 0), tcl = -0.4, mar = c(3.2 + sublines + 0.25 * (sublines > 0), 
                                                 3.5 + leftlines, 1 + toplines, 1) + 0.1)
  par(pch = 1)
  invisible()
}
################
# rainbow plots
################

savefig(file.path(dir.p,"Fig_1a"), width = 12, height = 10, type = "png", toplines = 0.8)
plot(fts(0:110, t(Japan_female_pop)), xlab = "Age", ylab = "Life-table death count", ylim = c(0, 5150))
dev.off()

savefig(file.path(dir.p,"Fig_1b"), width = 12, height = 10, type = "png", toplines = 0.8)
plot(fts(0:110, t(Japan_male_pop)), xlab = "Age", ylab = "", ylim = c(0, 5150))
dev.off()

############################################
# image plots (Kullback-Leibler divergence)
############################################
female_prefecture_dx <- readRDS(paste0(dir.r1, "female_prefecture_dx.rds"))
male_prefecture_dx   <- readRDS(paste0(dir.r1, "male_prefecture_dx.rds"))
total_prefecture_dx  <- readRDS(paste0(dir.r1, "total_prefecture_dx.rds"))
Japan_female_pop <- readRDS(paste0(dir.r1, "Japan_female_pop.rds"))
Japan_male_pop   <- readRDS(paste0(dir.r1, "Japan_male_pop.rds"))
Japan_total_pop   <- readRDS(paste0(dir.r1, "Japan_total_pop.rds"))

ages <- 0:110
n_age <- length(ages)
year <- 1975:2023
n_year <- length(year)

state = c("Hokkaido", 
          "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
          "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
          "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
          "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

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

savefig(file.path(dir.p,"Fig_2a"), width = 12, height = 10, type = "png", toplines = 0.5)
image(year, 1:47, t(female_KLdiv_year[rev(1:47),]), xlab = "Year", ylab = "Prefecture", 
      zlim = c(0,0.16), main = "Female", yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig(file.path(dir.p,"Fig_2b"), width = 12, height = 10, type = "png", toplines = 0.5)
image(year, 1:47, t(male_KLdiv_year[rev(1:47),]),   xlab = "Year", ylab = "", 
      zlim = c(0,0.16), main = "Male",yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig(file.path(dir.p,"Fig_2c"), width = 12, height = 10, type = "png", toplines = 0.5)
image(year, 1:47, t(total_KLdiv_year[rev(1:47),]),  xlab = "Year", ylab = "", 
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

savefig(file.path(dir.p,"Fig_2d"), width = 12, height = 10, type = "png", toplines = 0.5)
image(0:110, 1:47, t(female_KLdiv_age[rev(1:47),]), xlab = "Age", ylab = "Prefecture", zlim = c(0, 1.26),
      yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig(file.path(dir.p,"Fig_2e"), width = 12, height = 10, type = "png", toplines = 0.5)
image(0:110, 1:47, t(male_KLdiv_age[rev(1:47),]),   xlab = "Age", ylab = "", zlim = c(0, 1.26), yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

savefig(file.path(dir.p,"Fig_2f"), width = 12, height = 10, type = "png", toplines = 0.5)
image(0:110, 1:47, t(total_KLdiv_age[rev(1:47),]),  xlab = "Age", ylab = "", zlim = c(0, 1.26), yaxt = "n")
axis(2, at = seq(10, 40, by = 10), labels = c(40, 30, 20, 10))
dev.off()

