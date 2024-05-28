#####################################################################################
##
## Script name: Resampling raster objects based on preset size and resolution.
##
## Purpose of script: Provide two options to resample raster objects based on a preset 
## size and resolution of the required new object.
##
## Author: RGA Watson
##
## Date Created: 2024-01-29
##
##
## Notes: First determine the outline of the new raster (in this case 10, 40, -40, -25). 
## This is placed in a data.frame to automate the calculation of the number of columns 
## and rows in relation to the resolution. The new (sample) raster is then created and the  
## coordinate system fixed to this. An example raster is subsequently resampled using the 
## raster-function from the terra-package, using bilinearity. 
## The example raster is a Chl-a raster set to the same coordinates as the sample raster. 
## The land section (southern Africa) is switched to NA to ensure there are no blurry lines 
## along the coastline as a result of the bilinearity function. 
##   
## The function created at the end (ResSampl_func) requires control implementations. 
##
##
#####################################################################################
### packages & functions

# require(tidyverse)
# require(sf)
require(terra)

# source("functions/packages.R")       # loads up all the packages we need

#####################################################################################
### settings

### filenames & directories

# options(scipen = 6, digits = 4) # prefer non-scientific notation

#####################################################################################
###

# Outline the extent of the study area
Outline <- data.frame(
  Par = c("xmin", "xmax", "ymin", "ymax"),
  Coor = c(10, 40, -40, -25)
)

# Coordinates in a Spat format
e <- terra::ext(Outline[,2])

# Determine the resolution of the new raster
res = (1/60)*3 # resolution of 0.05°, or 3'

# Determine number of columns/rows for the new raster
x_res <- diff(c(Outline[1,2], Outline[2,2]))/res
y_res <- diff(c(Outline[3,2], Outline[4,2]))/res

# Create the resampling raster
i=x_res; j=y_res
r1 <- rast(matrix(runif(i*j), nrow=j, ncol=i))

# Fix the coordinates to the new raster
ext(r1) <- e

# Example raster to be resampled. Ensure land is NA, otherwise border become blurry 
#rt <- rast_clip[[10]]
#rt[rt<0] <- NA

# resample example raster using bilinearity
#rs <- terra::resample(rt, r1, method="bilinear")

# Display differences between rasters
#par(mfrow=c(2,1))
#plot(rt[[1]])
#plot(rs[[1]])



# Information from above is re-iterated here:
ext = e
#ras = rt
res = 0.05

# the resampling function based on the work above. Some rows of codes are reworked to smooth out the operations.
ResSampl_func <- function(ext,ras,res) { 
  
  # "ext" is a SpatExtent object of the extent of the sample raster
  # "ras" is the raster object that requires resampling
  # "res" is the resolution of the new raster object.
  
  x_res <- diff(c(ext[1], ext[2]))*(1/res)
  x_res <- x_res[[1]]
  y_res <- diff(c(ext[3], ext[4]))*(1/res)
  y_res <- y_res[[1]]
  
  r1 <- rast(matrix(runif(x_res*y_res), nrow=y_res, ncol=x_res))
  ext(r1) <- ext
  
  rs <- terra::resample(ras, r1, method="bilinear")
  rs
}

# Implementation of the function
# rt_new <- ResSampl_func(ext=e, ras=rt, res=0.05)




#####################################################################################
### unload packages

# detach("package:xxx", unload=TRUE)
