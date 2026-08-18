# This Rcode produces Figure 8 in the paper. It shows the equivalences between different mortality instruments.

source("load_packages.R")

# Dynamic Directories
dir.d       <- here("data", "updated_data")
dir.results <- here("results", "Results_Point")
dir.p       <- here("results", "Plots")
dir.aux     <- here("auxiliary_source")

# Ensure output directories exist
if (!dir.exists(dir.results)) dir.create(dir.results, recursive = TRUE)
if (!dir.exists(dir.p))       dir.create(dir.p, recursive = TRUE)



state = c("Hokkaido", 
          "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
          "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", 
          "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
          "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

# Define age and year ranges
ages <- 0:110
n_age <- length(ages)
years <- 1975:2023
n_year <- length(years)


# Convert from death table forecast to other mortality forecast for HPFCA for fh=1 for the 47 prefectures.

point_fore_subnational_HDFPCA_F_CDF <- readRDS(file.path(dir.results, "HDFPCA_Forecast_F1.rds"))
point_fore_subnational_HDFPCA_M_CDF <- readRDS(file.path(dir.results, "HDFPCA_Forecast_M1.rds"))

i=1 #Hokkaido

Forecast_F <- point_fore_subnational_HDFPCA_F_CDF[1,,]
Forecast_M <- point_fore_subnational_HDFPCA_M_CDF[1,,]



mx_F<-convertFx(ages,Forecast_F,from = "dx",to="mx")
mx_M<-convertFx(ages,Forecast_M,from = "dx",to="mx")

qx_F<-convertFx(ages,Forecast_F,from = "dx",to="qx")
qx_M<-convertFx(ages,Forecast_M,from = "dx",to="qx")

ex_F<-convertFx(ages,Forecast_F,from = "dx",to="ex")
ex_M<-convertFx(ages,Forecast_M,from = "dx",to="ex")


## ============================
## Common y-axis limits
## ============================

# dx: common scale for Female & Male
ylim_dx <- range(c(Forecast_F, Forecast_M), na.rm = TRUE)

# ex: common scale for Female & Male
ylim_ex <- range(c(ex_F, ex_M), na.rm = TRUE)

# (Optional) qx common scale if ever needed
# ylim_qx <- range(c(qx_F, qx_M), na.rm = TRUE)


## ============================
## Plot-saving function
## ============================

savefig <- function (filename, height = 10,
                     width = (1 + sqrt(5))/2 * height,
                     type = c("eps", "pdf", "jpg", "png"),
                     pointsize = 10, family = "Helvetica",
                     sublines = 0, toplines = 0, leftlines = 0,
                     res = 300) 
{
  type <- match.arg(type)
  filename <- paste(filename, ".", type, sep = "")
  
  if (type == "eps") {
    postscript(file = filename, horizontal = FALSE,
               width = width/2.54, height = height/2.54,
               pointsize = pointsize, family = family,
               onefile = FALSE, print.it = FALSE)
  } else if (type == "pdf") {
    pdf(file = filename, width = width/2.54,
        height = height/2.54, pointsize = pointsize,
        family = family, onefile = TRUE)
  } else if (type == "jpg") {
    jpeg(filename = filename, width = width, height = height,
         res = res, quality = 100, units = "cm")
  } else if (type == "png") {
    png(filename = filename, width = width, height = height,
        res = res, units = "cm")
  } else {
    stop("Unknown file type")
  }
  
  par(
    mgp = c(2.2, 0.45, 0),
    tcl = -0.4,
    mar = c(3.2 + sublines + 0.25 * (sublines > 0),
            3.5 + leftlines,
            1 + toplines, 1) + 0.1,
    pch = 1
  )
  
  invisible()
}


## ============================
## Forecasted dx
## ============================

# Female
savefig(paste0(dir.p, "/Forecast_dx_F"), type = "png")
matplot(
  ages, Forecast_F, type = "l",
  ylim = ylim_dx,
  main = "Forecasted dx (F) 2007–2023",
  xlab = "Age", ylab = "dx"
)
dev.off()

# Male
savefig(paste0(dir.p, "/Forecast_dx_M"), type = "png")
matplot(
  ages, Forecast_M, type = "l",
  ylim = ylim_dx,
  main = "Forecasted dx (M) 2007–2023",
  xlab = "Age", ylab = "dx"
)
dev.off()


## ============================
## Forecasted qx
## ============================

# Female
savefig(paste0(dir.p, "/Forecast_qx_F"), type = "png")
matplot(
  ages, qx_F, type = "l",
  main = "Forecasted qx (F) 2007–2023",
  xlab = "Age", ylab = "qx"
)
dev.off()

# Male
savefig(paste0(dir.p, "/Forecast_qx_M"), type = "png")
matplot(
  ages, qx_M, type = "l",
  main = "Forecasted qx (M) 2007–2023",
  xlab = "Age", ylab = "qx"
)
dev.off()


## ============================
## Forecasted ex
## ============================

# Female
savefig(paste0(dir.p, "/Forecast_ex_F"), type = "png")
matplot(
  ages, ex_F, type = "l",
  ylim = ylim_ex,
  main = "Forecasted ex (F) 2007–2023",
  xlab = "Age", ylab = "ex"
)
dev.off()

# Male
savefig(paste0(dir.p, "/Forecast_ex_M"), type = "png")
matplot(
  ages, ex_M, type = "l",
  ylim = ylim_ex,
  main = "Forecasted ex (M) 2007–2023",
  xlab = "Age", ylab = "ex"
)
dev.off()
