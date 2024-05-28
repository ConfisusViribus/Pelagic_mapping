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
crs(EEZ_vect) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
#Env_EEZ_2 <- mask(Env_EEZ, EEZ_vect)




comb_sel <- comb[[c("SST_mean", "SST_max", "SST_var", "SST_sobel",
                    #"Chl_mean", "Chl_var", #"Chl_sobel", ## Chl_Sobel removed because of large blocks of missing data.
                    "NPP_mean", "NPP_var", "NPP_sobel",
                    "Curr_mean", "Curr_var", "Curr_sobel",
                    "Tur_mean", "Tur_var", "Tur_sobel"
                    #"Slope",
                    #"OkuboWeiss"
                  )]]
Bath_lyr <- comb[["Bathymetry"]]
Bath_lyr[Bath_lyr < -1000] <- -1000
Bath_cont <- as.contour(Bath_lyr)
Bath_line <- Bath_cont[values(Bath_cont)==-1000]
#Bath_lyr[Bath_lyr < -1000] <- -1000
Bath_lyr <- math((abs(Bath_lyr)+1), "log")

comb_sel_st <- scale(comb_sel, center = T, scale=T)
comb_sel_2 <- c(comb_sel_st, Bath_lyr)

msk <- ext(15, 20, -40, -35)
comb_crop <- crop(comb_sel_2, msk)
comb_crop <- mask(comb_sel_2, EEZ_vect)


num_centers <- 15



pca_rasts <- comb_sel_2
#pca_rasts <- comb_sel_st
pca_rasts <- comb_crop

pca_result <- prcomp(pca_rasts) 
pca_result_2 <- princomp(pca_rasts)

p3 <- predict(pca_rasts, pca_result)
k3 <- k_means(p3, centers = num_centers, iter.max=50)


# Calculate the eigenvalues of individual PCs. 
ev <- pca_result$sdev^2
fviz_eig(pca_result)
ev_1 <- ev[ev>1]
format(ev_1, scientific = FALSE)


fviz_nbclust(as.data.frame(pca_rasts, na.rm=T), kmeans, method = "silhouette", k.max = 40, iter.max = 50)


pca_data <- as.data.frame(pca_result$x)

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Results", x = "Principal Component 1", y = "Principal Component 2", color = "Gradient1")


kmeans_results_pca <- kmeans(as.data.frame(pca_result$x), centers = num_centers, iter.max=100)
pca_result[["Kmeans"]] <- kmeans_results_pca$cluster
comb_sel_2[["Kmeans"]] <- kmeans_results_pca$cluster

kmeans_result <- k_means(terra::na.omit(pca_rasts), centers = num_centers, iter.max=100)  # You can adjust the number of centers as needed
pca_data["Kmeans"] <- as.data.frame(kmeans_result$lyr1)

ggplot() +
   #geom_point(data=pca_data[,-"Kmeans"],aes(PC1,PC2),color="grey", alpha=0.3)+
   geom_point(data=pca_data,aes(PC1,PC2,color=as.factor(Kmeans)))+
   #facet_wrap(~ Kmeans) + 
   labs(title = "PCA Results", x = "Principal Component 1", y = "Principal Component 2", color = "Gradient1")

#data.frame(kmeans_result$centers)
par(mfrow=c(1,2))
plot(k3, main="SpatRaster > PCA > Kmeans", col = pal(25))
lines(EEZ)
plot(kmeans_result, main="SpatRaster > Kmeans", col = pal(16))
lines(EEZ)

pca_rasts_post <- pca_rasts
pca_rasts_post[["Kmeans"]] <- kmeans_result$lyr1


f2 <- focal(kmeans_result, w=3, modal, na.rm=TRUE)
f2 <- focal(f2, w=3, modal, na.rm=TRUE)
plot(f2)


library(RColorBrewer)

cols <- brewer.pal(9, "Set1")
cols2 <- brewer.pal(12, "Set3")
cols <- c(cols, cols2)
pal <- colorRampPalette(cols)

plot(pca_rasts_post$Kmeans, col = pal(15))
lines(EEZ)
plot(f2, col = pal(15))
f2_EEZ <- mask(pca_rasts_post$Kmeans, EEZ_vect)
plot(f2_EEZ,col = pal(15))


clusplot(comb_sel_2, kmeans_result, color=TRUE, shade=TRUE, 
         labels=2, lines=0)




