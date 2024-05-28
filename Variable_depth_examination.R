library(terra)
require(tidyverse)
require(sf)
library(factoextra)
setwd("~/Documents/SANBI_Pelagic_Ecosystems_PostDoc")

files <- dir(paste("R Project/Variable_Maps/"), pattern="*.nc$")
files <- paste("R Project/Variable_Maps/", files, sep="")
comb <- rast(files)

#dir("data", "EEZ_South_Africa_buffered_beyond_allEEZversions1.shp", full.names = T, recursive = T_

EEZ_vect <- vect(EEZ)

num_centers <- 15

comb_sel <- comb[[c("SST_mean", "SST_max", "SST_var", "SST_sobel",
                    #"Chl_mean", "Chl_var", #"Chl_sobel", ## Chl_Sobel removed because of large blocks of missing data.
                    "NPP_mean", "NPP_var", "NPP_sobel",
                    "Curr_mean", "Curr_var", "Curr_sobel",
                    "Tur_mean", "Tur_var", "Tur_sobel"
                    #"Slope",
                    #"OkuboWeiss"
)]]


Bath_lyr <- comb[["Bathymetry"]]
Bath_cont <- as.contour(Bath_lyr, nlevels=50)
bath_line_10 <- Bath_cont[values(Bath_cont)==-1000]
bath_line_2 <- Bath_cont[values(Bath_cont)==-200]

Bath_lyr_10 <- Bath_lyr
Bath_lyr_10[Bath_lyr_10 < -1000] <- -1000
Bath_lyr_10_na <- Bath_lyr
Bath_lyr_10_na[Bath_lyr_10_na < -1000] <- NA

Bath_lyr_2 <- Bath_lyr
Bath_lyr_2[Bath_lyr_2 < -200] <- -200
Bath_lyr_2_na <- Bath_lyr
Bath_lyr_2_na[Bath_lyr_2_na < -200] <- NA

Bath_lyr <- math((abs(Bath_lyr)+1), "log")
Bath_lyr_10 <- math((abs(Bath_lyr_10)+1), "log")
Bath_lyr_10_na <- math((abs(Bath_lyr_10_na)+1), "log")
Bath_lyr_2 <- math((abs(Bath_lyr_2)+1), "log")
Bath_lyr_2_na <- math((abs(Bath_lyr_2_na)+1), "log")

comb_sel_st <- scale(comb_sel, center = T, scale=T)
comb_sel_2 <- c(comb_sel_st, Bath_lyr)
comb_sel_2_10 <- c(comb_sel_st, Bath_lyr_10)
comb_sel_2_10_na <- c(comb_sel_st, Bath_lyr_10_na)
comb_sel_2_2 <- c(comb_sel_st, Bath_lyr_2)
comb_sel_2_2_na <- c(comb_sel_st, Bath_lyr_2_na)

pca_rasts_all <- comb_sel_2
pca_result_all <- prcomp(pca_rasts_all) 
pred_rast_all <- predict(pca_rasts_all, pca_result_all)
allD <- k_means(pred_rast_all, centers = num_centers, iter.max=50)

pca_rasts_10 <- comb_sel_2_10
pca_result_10 <- prcomp(pca_rasts_10) 
pred_rast_10 <- predict(pca_rasts_10, pca_result_10)
all_10 <- k_means(pred_rast_10, centers = num_centers, iter.max=50)

# pca_rasts_10_na <- comb_sel_2_10_na
# pca_result_10_na <- terra::prcomp(na.exclude(pca_rasts_10_na))
# pred_rast_10_na <- terra::predict(pca_rasts_10_na, pca_result_10_na, na.rm=T)
# all_10_na <- k_means(pred_rast_10_na, centers = num_centers, iter.max=50)

pca_rasts_2 <- comb_sel_2_2
pca_result_2 <- prcomp(pca_rasts_2) 
pred_rast_2 <- predict(pca_rasts_2, pca_result_2)
all_2 <- k_means(pred_rast_2, centers = num_centers, iter.max=50)

# pca_rasts_4_na <- comb_sel_2_4_na
# pca_result_4_na <- prcomp(pca_rasts_4_na) 
# pred_rast_4_na <- predict(pca_rasts_4_na, pca_result_4_na)
# all_4_na <- k_means(pred_rast_4_na, centers = num_centers, iter.max=50)

pca_rasts_no <- comb_sel_st
pca_result_no <- prcomp(pca_rasts_no) 
pred_rast_no <- predict(pca_rasts_no, pca_result_no)
noD <- k_means(pred_rast_no, centers = num_centers, iter.max=50)



dev.off()
par(mfrow=c(2,2))
plot(allD, main="Analysis including all depths.", col = pal(15))
lines(bath_line_10)
lines(bath_line_2)
lines(EEZ, col="blue")

plot(noD, main="Analysis with no bathymetry.", col = pal(15))
lines(bath_line_10)
lines(bath_line_2)
lines(EEZ, col="blue")

plot(all_10, main="Analysis including bathymetry down to -1000m (flattened).", col = pal(15))
lines(bath_line_10)
lines(EEZ, col="blue")
# plot(all_10_na, main="Analysis including bathymetry down to -1000m (cut-off).", col = pal(15))
# lines(bath_line_10)
# lines(EEZ, col="blue")

plot(all_2, main="Analyis including bathymetry down to -200m (flattened).", col = pal(15))
lines(bath_line_2)
lines(EEZ, col="blue")
# plot(all_4_na, main="Analyis including bathymetry down to -400m (cut-off).", col = pal(15))
# lines(bath_line_4)
# lines(EEZ, col="blue")

