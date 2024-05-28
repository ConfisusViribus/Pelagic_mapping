#####################################################################################
##
## Script name: Data Import 
##
## Purpose of script: Environmental data import data
## 
## Author: RGA Watson
##
## Date Created: 2024-03-05
##
##
## Notes:
##  Complications surrounding the handling of the landmass. NA generates problems in the analyses. 
##
#####################################################################################
### packages & functions

require(tidyverse)
# require(sf)
require(terra)
library(svMisc) # For the progress bar to track the Sobel analyses for CHL

# require(raster)    
# library(sf)
# library(grec)
# source("functions/packages.R")       # loads up all the packages we need

#####################################################################################
### settings

### filenames & directories

setwd("~/Documents/SANBI_Pelagic_Ecosystems_PostDoc")
source("R Project/pelagic_map/resampling_function.R")

# options(scipen = 6, digits = 4) # prefer non-scientific notation

#####################################################################################
###

# The resolution set for resampling
new_res <- 0.05

#### Import Bathymetry ####
# Source: GEBCO 2023

print("Start Bathymetry process.")

b = terra::rast("RAW_Data/GEBCO_2023/gebco_2023_n-25.0_s-40.0_w10.0_e40.0.nc")
b <- ifel(b>0, NA, b)
ext(b) <- c(10, 40, -40, -25)
crs(b) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
b_new <- ResSampl_func(ext=e, ras=b, res=new_res)
writeCDF(b_new, "R Project/Variable_Maps/Bathymetry.nc", overwrite=TRUE, varname="Bathymetry", unit = "m")

#b_new[b_new < -200] <- -200
#log_b_new <- math((abs(b_new)+1), "log")


#### Slope ####
# Source: Derived from GEBCO 2023

print("Start Slope process.")
sl <- terrain(b, v="slope", neighbors=8, unit="degrees")

sl_new <- ResSampl_func(ext=e, ras=sl, res=new_res)
writeCDF(sl_new, "R Project/Variable_Maps/Slope.nc", overwrite=TRUE, varname="Slope", unit = "Deg")


#### Import SST ####
# Source: Copernicus 
# Time from Jan 2011 - Dec 2021

print("Start SST process.")

s = terra::rast("RAW_Data/aphiwe_layers/SST/Data/METOFFICE-GLO-SST-L4-REP-OBS-SST_multi-vars_10.02E-39.97E_39.97S-25.02S_2011-01-01-2021-12-31.nc", subds = "analysed_sst")
s <- s - 273.15 
ext(s) <- c(10, 40, -40, -25)
crs(s) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

print("SST mean")
#plot(s[[100]])
s_mean = terra::mean(s, na.rm=T)
s_mean_new <- ResSampl_func(ext=e, ras=s_mean, res=new_res)
writeCDF(s_mean_new, "R Project/Variable_Maps/SST_mean.nc", overwrite=TRUE, varname="SST_mean", unit = "°C")

dts = time(s) |> lubridate::ymd()

yearly_i = paste(
  lubridate::year(dts),
  sep = '-') |> factor()

print("SST min")

yearly_s_min = terra::tapp(s, index = yearly_i, fun = min)
s_min <- mean(yearly_s_min, na.rm=T)
s_min_new <- ResSampl_func(ext=e, ras=s_min, res=new_res)
writeCDF(s_min_new, "R Project/Variable_Maps/SST_min.nc", overwrite=TRUE, varname="SST_min", unit = "°C")
#s_min_2 = min(s, na.rm=T)

# par(mfrow = c(1,2))
# plot(s_min, main="Minimum SST, averaged per year.")
# plot(s_min_2, main="Absolute minimum SST.")

print("SST max")
yearly_s_max = terra::tapp(s, index = yearly_i, fun = max)
s_max = mean(yearly_s_max, na.rm=T)
s_max_new <- ResSampl_func(ext=e, ras=s_max, res=new_res)
writeCDF(s_max_new, "R Project/Variable_Maps/SST_max.nc", overwrite=TRUE, varname="SST_max", unit = "°C")
print("SST var")
yearly_s_var = terra::tapp(s, index = yearly_i, fun = sd)^2
s_var = mean(yearly_s_var, na.rm=T)
s_var_new <- ResSampl_func(ext=e, ras=s_var, res=new_res)
writeCDF(s_var_new, "R Project/Variable_Maps/SST_var.nc", overwrite=TRUE, varname="SST_var", unit = "°C")
print("SST range")
s_range = s_max - s_min
s_range_new <- ResSampl_func(ext=e, ras=s_range, res=new_res)
writeCDF(s_range_new, "R Project/Variable_Maps/SST_range.nc", overwrite=TRUE, varname="SST_range", unit = "°C")

print("Start SST Sobel process.")
# monthly stack (sum of daily values)
# make a monthly index first
weekly_i = paste(
  lubridate::isoweek(dts),
  lubridate::year(dts),
  sep = '-') |> factor()



# aggregate the stack
weekly_s_mean = terra::tapp(s, index = weekly_i, fun = mean)
# monthly_s_max = terra::tapp(s, index = monthly_i, fun = max)
# monthly_s_var = terra::tapp(s, index = monthly_i, fun = var)
# monthly_s_range = terra::tapp(s, index = monthly_i, fun = range)

stack_l <- nlyr(weekly_s_mean)
SST_sobel_stack <- NULL

for(i in c(1:stack_l)){
  w = matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  sobel_st <- terra::focal(weekly_s_mean[[i]], w = w, na.policy = "omit")
  vsobel_st <- terra::focal(weekly_s_mean[[i]], w = t(w), na.policy = "omit")
  comb_sobel_st <- sqrt(sobel_st^2 + vsobel_st^2)
  SST_sobel_stack <- c(SST_sobel_stack,comb_sobel_st)
  
  progress(i, stack_l)
  
  rm(sobel_st)
  rm(vsobel_st)
  rm(comb_sobel_st)
  gc()
  
}

SST_sobel_stack <- rast(SST_sobel_stack)
#SST_sobel_stack_sum <- app(SST_sobel_stack, fun=sum)
SST_sobel_stack_sum <- terra::app(SST_sobel_stack,fun=sum, na.rm=T)
writeCDF(SST_sobel_stack_sum, "R Project/Variable_Maps/SST_sobel.nc", overwrite=TRUE, varname="SST_sobel", unit = "°C")

# raster01 = function(r){
#   # get the min max values
#   minmax_r = range(values(r), na.rm=TRUE) 
#   
#   # rescale 
#   return( (r-minmax_r[1]) / (diff(minmax_r)))
# }
# SST_sobel_stack_norm <- raster01(SST_sobel_stack)
# plot(sum(max(SST_sobel_stack_norm))>0.4)

#SST_sobel_stack_std <- scale(SST_sobel_stack, center = T, scale=T)
#SST_fr_sum <- app(SST_sobel_stack_std, sum)

# weekly_s_sobel <- focal(weekly_s_mean[[1]], w =matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3), na.rm=FALSE)
# 
# tapp(weekly_s_mean, )
# weekly_sst_freq <- terra::terrain(weekly_s_mean, v="slope" )
# plot(terra::boundaries(weekly_s_mean, directions=4, classes=TRUE))



#### Import Chl-A ####
# Source: Copernicus 
# Time from Jan 2011 - Dec 2021
print("Start Chl process.")
# c = terra::rast("RAW_Data/aphiwe_layers/Chl/Data/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_CHL_10.02E-39.98E_39.98S-25.02S_2011-01-01-2021-12-31.nc")
# c_nlyr <- nlyr(c)
# ext(c) <- c(10, 40, -40, -25)
# crs(c) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# c_new <- NULL
# 
# # performed this way to keep track of progress
# for(i in 1:c_nlyr){
#   c_sub <- c[[i]]
#   c_sub[c_sub>10] <- 10
#   c_new <- c(c_new,c_sub)
# 
#   progress(i, c_nlyr)
# }
# 
# c_new <- rast(c_new)
# 
# writeCDF(c_new, "RAW_Data/aphiwe_layers/Chl/Data/Chl_a_adjusted.nc", overwrite=TRUE, varname="Chl_a_adj",
#          longname="CHL (Chlorophyll-a concentration - Mean of the binned pixels)", unit="milligram m-3")
c_new = terra::rast("RAW_Data/aphiwe_layers/Chl/Data/Chl_a_adjusted.tif")
crs(c_new) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# plot(c[[100]])
# 
# yearly_c_mean = terra::tapp(c, index = yearly_b, fun = mean)
# yearly_c_mean_new <- ResSampl_func(ext=e, ras=yearly_c_mean, res=new_res)
# yearly_c_min = terra::tapp(c, index = yearly_b, fun = min)
# yearly_c_min_new <- ResSampl_func(ext=e, ras=yearly_c_min, res=new_res)
# yearly_c_max = terra::tapp(c, index = yearly_b, fun = max)
# yearly_c_max_new <- ResSampl_func(ext=e, ras=yearly_c_max, res=new_res)
# yearly_c_var = terra::tapp(c, index = yearly_b, fun = var)
# yearly_c_var_new <- ResSampl_func(ext=e, ras=yearly_c_var, res=new_res)
# yearly_c_range = yearly_c_max-yearly_c_min
# yearly_c_range_new <- ResSampl_func(ext=e, ras=yearly_c_range, res=new_res)

print("Chl mean")
yearly_c_mean = terra::tapp(c_new, index = yearly_i, fun = mean)
#c_mean = terra::mean(yearly_c_mean, na.rm=T)
c_mean_new <- ResSampl_func(ext=e, ras=yearly_c_mean, res=new_res)
writeCDF(c_mean_new, "R Project/Variable_Maps/Chl_mean.nc", overwrite=TRUE, varname="Chl_mean", unit = "mg m-3")

print("Chl min")
yearly_c_min = terra::tapp(c_new, index = yearly_i, fun = min)
c_min = mean(yearly_c_min, na.rm=T)
c_min_new <- ResSampl_func(ext=e, ras=c_min, res=new_res)
writeCDF(c_min_new, "R Project/Variable_Maps/Chl_min.nc", overwrite=TRUE, varname="Chl_min", unit = "mg m-3")

print("Chl max")
yearly_c_max = terra::tapp(c_new, index = yearly_i, fun = max)
c_max = mean(yearly_c_max, na.rm=T)
c_max_new <- ResSampl_func(ext=e, ras=c_max, res=new_res)
writeCDF(c_max_new, "R Project/Variable_Maps/Chl_max.nc", overwrite=TRUE, varname="Chl_max", unit = "mg m-3")

print("Chl var")
yearly_c_var = terra::tapp(c_new, index = yearly_i, fun = sd)^2
c_var = mean(yearly_c_var, na.rm=T)
c_var_new <- ResSampl_func(ext=e, ras=c_var, res=new_res)
writeCDF(c_var_new, "R Project/Variable_Maps/Chl_var.nc", overwrite=TRUE, varname="Chl_var", unit = "mg m-3")

print("Chl range")
c_range = c_max - c_min
c_range_new <- ResSampl_func(ext=e, ras=c_range, res=new_res)
writeCDF(c_range_new, "R Project/Variable_Maps/Chl_range.nc", overwrite=TRUE, varname="Chl_range", unit = "mg m-3")


print("Start Chl Sobel process.")

dts = time(c_new) |> lubridate::ymd()
weekly_i = paste(
  lubridate::isoweek(dts),
  lubridate::year(dts),
  sep = '-') |> factor()

weekly_c_mean = terra::tapp(c_new, index = weekly_i, fun = mean)
weekly_c_mean_new <- ResSampl_func(ext=e, ras=weekly_c_mean, res=new_res)

c_stack_l <- nlyr(weekly_c_mean)
Chl_sobel_stack <- NULL

for(i in c(1:c_stack_l)){
  w = matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  c_sobel_st <- focal(weekly_c_mean_new[[i]], w = w)
  c_vsobel_st <- focal(weekly_c_mean_new[[i]], w = t(w))
  comb_c_sobel_st <- sqrt(c_sobel_st^2 + c_vsobel_st^2)
  Chl_sobel_stack <- c(Chl_sobel_stack,comb_c_sobel_st)
  
  progress(i, c_stack_l)
  
  rm(c_sobel_st)
  rm(c_vsobel_st)
  gc()
}

Chl_sobel_stack <- rast(Chl_sobel_stack)
Chl_sobel_stack_sum <- app(Chl_sobel_stack, fun=sum)
writeCDF(Chl_sobel_stack_sum, "R Project/Variable_Maps/Chl_sobel.nc", overwrite=TRUE, varname="Chl_sobel", unit = "mg m-3")

#Chl_sobel <- calc(Chl_sobel_stack, function(x) sum(x, na.rm=T)/length(x))
#Chl_sobel <- sum(Chl_sobel_stack, na.rm=T)

#chl_freq_stack <- freq(Chl_sobel_stack)

#par(mfrow = c(1,2))


#### Import Sea Surface Height #####
# Source: Copernicus 
# Time from Jan 2011 - Dec 2021

print("Start SSH process.")

h = terra::rast("RAW_Data/aphiwe_layers/SSH/Data/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_sla-ugos-vgos_10.12E-39.88E_39.88S-25.12S_2011-01-01-2021-12-31.nc", subds = "sla")
#h[is.na(h)] <-  -5000
#b <- b - 273.15 
#b[b==0] <- NA
#ext(b) <- c(10, 40, -40, -25)
crs(h) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# plot(h[[100]])
# 
# yearly_h_mean = terra::tapp(h, index = yearly_b, fun = mean)
# yearly_h_mean_new <- ResSampl_func(ext=e, ras=yearly_h_mean, res=new_res)
# yearly_h_min = terra::tapp(h, index = yearly_b, fun = min)
# yearly_h_min_new <- ResSampl_func(ext=e, ras=yearly_h_min, res=new_res)
# yearly_h_max = terra::tapp(h, index = yearly_b, fun = max)
# yearly_h_max_new <- ResSampl_func(ext=e, ras=yearly_h_max, res=new_res)
# yearly_h_var = terra::tapp(h, index = yearly_b, fun = var)
# yearly_h_var_new <- ResSampl_func(ext=e, ras=yearly_h_var, res=new_res)
# yearly_h_range = yearly_h_max-yearly_h_min
# yearly_h_range_new <- ResSampl_func(ext=e, ras=yearly_h_range, res=new_res)

print("Start Okubo-Weiss process.")

dts_h = time(h) |> lubridate::ymd()
weekly_i = paste(
  lubridate::isoweek(dts_h),
  lubridate::year(dts_h),
  sep = '-') |> factor()
yearly_i_h = paste(
  lubridate::year(dts_h),
  sep = '-') |> factor()

weekly_h_mean = terra::tapp(h, index = weekly_i, fun = mean)
weekly_h_mean_new <- ResSampl_func(ext=e, ras=weekly_h_mean, res=new_res)

h_stack_l <- nlyr(weekly_h_mean_new)
Okubo_W_stack <- NULL

for(i in c(1:h_stack_l)){
  
  slope <- terrain(weekly_h_mean_new[[i]], "slope")
  
  # Calculate aspect
  aspect <- terrain(weekly_h_mean_new[[i]], "aspect")
  
  # Calculate Okubo-Weiss parameter
  okubo_weiss <- slope^2 - aspect^2
  
  Okubo_W_stack <- c(Okubo_W_stack,okubo_weiss)
  
  progress(i, h_stack_l)
  
  rm(okubo_weiss)
  rm(slope)
  rm(aspect)
  gc()
  
}

Okubo_W_stack <- rast(Okubo_W_stack)
Okubo_W_stack_mean = terra::mean(Okubo_W_stack, na.rm=T)
Okubo_W_stack_mean = terra::app(Okubo_W_stack, fun=sum)
Okubo_W_stack_mean_new <- ResSampl_func(ext=e, ras=Okubo_W_stack_mean, res=new_res)
writeCDF(Okubo_W_stack_mean_new, "R Project/Variable_Maps/OkuboWeiss.nc", overwrite=TRUE, varname="OkuboWeiss", unit = "m")

# Plot Okubo-Weiss parameter
#plot(Okubo_W_stack[[4]], main="Okubo-Weiss Parameter")



print("SSH mean")
h_mean = terra::mean(h, na.rm=T)
h_mean_new <- ResSampl_func(ext=e, ras=h_mean, res=new_res)
writeCDF(h_mean_new, "R Project/Variable_Maps/SSH_mean.nc", overwrite=TRUE, varname="SSH_mean", unit = "mg m-3")
print("SSH min")
yearly_h_min = terra::tapp(h, index = yearly_i_h, fun = min)
h_min = mean(yearly_h_min, na.rm=T)
h_min_new <- ResSampl_func(ext=e, ras=h_min, res=new_res)
writeCDF(h_min_new, "R Project/Variable_Maps/SSH_min.nc", overwrite=TRUE, varname="SSH_min", unit = "mg m-3")
print("SSH max")
yearly_h_max = terra::tapp(h, index = yearly_i_h, fun = max)
h_max = mean(yearly_h_max, na.rm=T)
h_max_new <- ResSampl_func(ext=e, ras=h_max, res=new_res)
writeCDF(h_max_new, "R Project/Variable_Maps/SSH_max.nc", overwrite=TRUE, varname="SSH_max", unit = "mg m-3")
print("SSH var")
yearly_h_var = terra::tapp(h, index = yearly_i_h, fun = sd)^2
h_var = terra::mean(yearly_h_var, na.rm=T)^2
h_var_new <- ResSampl_func(ext=e, ras=h_var, res=new_res)
writeCDF(h_var_new, "R Project/Variable_Maps/SSH_var.nc", overwrite=TRUE, varname="SSH_var", unit = "mg m-3")
print("SSH range")
h_range = h_max - h_min
h_range_new <- ResSampl_func(ext=e, ras=h_range, res=new_res)
writeCDF(h_range_new, "R Project/Variable_Maps/SSH_range.nc", overwrite=TRUE, varname="SSH_range", unit = "mg m-3")


#### Import Current ####
# Source: Copernicus 
# Time from Jan 2011 - Dec 2021
print("Start Current process.")
# ugos is eastward sea water velocity
hc = terra::rast("RAW_Data/aphiwe_layers/SSH/Data/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_sla-ugos-vgos_10.12E-39.88E_39.88S-25.12S_2011-01-01-2021-12-31.nc", subds = "ugos")
crs(hc) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# print("Horizontal current process")
# hc_mean = terra::mean(hc, na.rm=T)
# hc_mean_new <- ResSampl_func(ext=e, ras=hc_mean, res=new_res)
# hc_min = min(hc, na.rm=T)
# hc_min_new <- ResSampl_func(ext=e, ras=hc_min, res=new_res)
# hc_max = max(hc, na.rm=T)
# hc_max_new <- ResSampl_func(ext=e, ras=hc_max, res=new_res)
# hc_var = terra::stdev(hc, na.rm=T)^2
# hc_var_new <- ResSampl_func(ext=e, ras=hc_var, res=new_res)
# hc_range = hc_max - hc_min
# hc_range_new <- ResSampl_func(ext=e, ras=hc_range, res=new_res)

# vgos is northward sea water velocity
# Source: Copernicus 
# Time from Jan 2011 - Dec 2021

vc = terra::rast("RAW_Data/aphiwe_layers/SSH/Data/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_sla-ugos-vgos_10.12E-39.88E_39.88S-25.12S_2011-01-01-2021-12-31.nc", subds = "vgos")
crs(vc) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# vc_mean = terra::mean(vc, na.rm=T)
# vc_mean_new <- ResSampl_func(ext=e, ras=vc_mean, res=new_res)
# vc_min = min(vc, na.rm=T)
# vc_min_new <- ResSampl_func(ext=e, ras=vc_min, res=new_res)
# vc_max = max(vc, na.rm=T)
# vc_max_new <- ResSampl_func(ext=e, ras=vc_max, res=new_res)
# vc_var = terra::stdev(vc, na.rm=T)^2
# vc_var_new <- ResSampl_func(ext=e, ras=vc_var, res=new_res)
# vc_range = vc_max - vc_min
# vc_range_new <- ResSampl_func(ext=e, ras=vc_range, res=new_res)


curr_l <- nlyr(hc)

curr_stack <- NULL

for(i in c(1:curr_l)){
  curr_stack_single <-sqrt(hc[[i]]^2 + vc[[i]]^2)
  
  curr_stack <- c(curr_stack,curr_stack_single)
  
  progress(i, curr_l)
  
  rm(curr_stack_single)
  gc()
}

curr_stack <- rast(curr_stack)
yearly_h_min = terra::tapp(h, index = yearly_i_h, fun = min)
h_min = mean(yearly_h_min, na.rm=T)

print("Current mean")
curr_mean = terra::mean(curr_stack, na.rm=T)
curr_mean_new <- ResSampl_func(ext=e, ras=curr_mean, res=new_res)
writeCDF(curr_mean_new, "R Project/Variable_Maps/Curr_mean.nc", overwrite=TRUE, varname="Curr_mean", unit = "mg m-3")
print("Current min")
curr_min = min(curr_stack, na.rm=T)
curr_min_new <- ResSampl_func(ext=e, ras=curr_min, res=new_res)
writeCDF(curr_min_new, "R Project/Variable_Maps/Curr_min.nc", overwrite=TRUE, varname="Curr_min", unit = "mg m-3")
print("Current max")
curr_max = max(curr_stack, na.rm=T)
curr_max_new <- ResSampl_func(ext=e, ras=curr_max, res=new_res)
writeCDF(curr_max_new, "R Project/Variable_Maps/Curr_max.nc", overwrite=TRUE, varname="Curr_max", unit = "mg m-3")
print("Current var")
curr_var = terra::stdev(curr_stack, na.rm=T)^2
curr_var_new <- ResSampl_func(ext=e, ras=curr_var, res=new_res)
writeCDF(curr_var_new, "R Project/Variable_Maps/Curr_var.nc", overwrite=TRUE, varname="Curr_var", unit = "mg m-3")
print("Current range")
curr_range = curr_max - curr_min
curr_range_new <- ResSampl_func(ext=e, ras=curr_range, res=new_res)
writeCDF(curr_range_new, "R Project/Variable_Maps/Curr_range.nc", overwrite=TRUE, varname="Curr_range", unit = "mg m-3")




cur_stack_l <- nlyr(curr_stack)
cur_sobel_stack <- NULL

for(i in c(1:cur_stack_l)){
  w = matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  cur_sobel_st <- focal(curr_stack[[i]], w = w, na.rm=TRUE)
  cur_vsobel_st <- focal(curr_stack[[i]], w = t(w), na.rm=TRUE)
  comb_cur_sobel_st <- sqrt(cur_sobel_st^2 + cur_vsobel_st^2)
  cur_sobel_stack <- c(cur_sobel_stack,comb_cur_sobel_st)
  
  progress(i, cur_stack_l)
  
  rm(cur_sobel_st)
  rm(cur_vsobel_st)
  gc()
}

cur_sobel_stack <- rast(cur_sobel_stack)
cur_sobel_stack_sum <- app(cur_sobel_stack, fun=sum, na.rm=T)
ext(cur_sobel_stack_sum) <- c(10, 40, -40, -25)
cur_sobel_stack_sum <- ResSampl_func(ext=e, ras=cur_sobel_stack_sum, res=new_res)
writeCDF(cur_sobel_stack_sum, "R Project/Variable_Maps/Current_sobel.nc", overwrite=TRUE, varname="Curr_sobel", unit = "mg m-3")





#### Import NPP ####
# Source: Oregon State
# Timeframe: Jan 2012 -  Dec 2022  ## Needs to be fixed to include 2011
print("Start NPP process.")
files <- dir(paste("RAW_Data/NPP/Standard VGPM - Monthly/VGPM_Monthly_1222/"), pattern="*.hdf$")
files <- paste("RAW_Data/NPP/Standard VGPM - Monthly/VGPM_Monthly_1222/", files, sep="")
list <- NULL
region <- ext(10, 40, -40, -25)

list <- rast(files)
ext(list) <- c(-180, 180, -90, 90)
crs(list) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
npp_comb <- crop(list,  region)
npp_comb[npp_comb<0] <- NA
#rast_clip[rast_clip == -9999] <- -1
#rast_clip_st <- scale(npp_comb, center = T, scale=T)

print("NPP mean")
npp_mean = terra::mean(npp_comb, na.rm=T)
npp_mean_new <- ResSampl_func(ext=e, ras=npp_mean, res=new_res)
writeCDF(npp_mean_new, "R Project/Variable_Maps/NPP_mean.nc", overwrite=TRUE, varname="NPP_mean", unit = "")
print("NPP min")
npp_min = min(npp_comb, na.rm=T)
npp_min_new <- ResSampl_func(ext=e, ras=npp_min, res=new_res)
writeCDF(npp_min_new, "R Project/Variable_Maps/NPP_min.nc", overwrite=TRUE, varname="NPP_min", unit = "")
print("NPP max")
npp_max = max(npp_comb, na.rm=T)
npp_max_new <- ResSampl_func(ext=e, ras=npp_max, res=new_res)
writeCDF(npp_max_new, "R Project/Variable_Maps/NPP_max.nc", overwrite=TRUE, varname="NPP_max", unit = "")
print("NPP var")
npp_var = terra::stdev(npp_comb, na.rm=T)^2
npp_var_new <- ResSampl_func(ext=e, ras=npp_var, res=new_res)
writeCDF(npp_var_new, "R Project/Variable_Maps/NPP_var.nc", overwrite=TRUE, varname="NPP_var", unit = "")
print("NPP range")
npp_range = npp_max - npp_min
npp_range_new <- ResSampl_func(ext=e, ras=npp_range, res=new_res)
writeCDF(npp_range_new, "R Project/Variable_Maps/NPP_range.nc", overwrite=TRUE, varname="NPP_range", unit = "")

npp_stack_l <- nlyr(npp_comb)
npp_sobel_stack <- NULL

for(i in c(1:npp_stack_l)){
  w = matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  npp_sobel_st <- terra::focal(npp_comb[[i]], w = w, na.rm=TRUE)
  npp_vsobel_st <- terra::focal(npp_comb[[i]], w = t(w), na.rm=TRUE)
  comb_npp_sobel_st <- sqrt(npp_sobel_st^2 + npp_vsobel_st^2)
  npp_sobel_stack <- c(npp_sobel_stack,comb_npp_sobel_st)
  
  progress(i, npp_stack_l)
  
  rm(npp_sobel_st)
  rm(npp_vsobel_st)
  gc()
}

npp_sobel_stack <- rast(npp_sobel_stack)
npp_sobel_stack_sum <- app(npp_sobel_stack, fun=sum)
npp_sobel_stack_sum <- ResSampl_func(ext=e, ras=npp_sobel_stack_sum, res=new_res)
writeCDF(npp_sobel_stack_sum, "R Project/Variable_Maps/NPP_sobel.nc", overwrite=TRUE, varname="NPP_sobel", unit = "mg m-3")


#### Turbidity #### 
# Source: Ocean Data NASA
# Timeframe: 2012 - 2022 ## Needs to be fixed to include 2011
print("Start Turbidity process.")
Tur_files <- dir(paste("RAW_Data/Turbidity_K490_OceanDataNASA/requested_files/"), pattern="*.png$")
Tur_files <- paste("RAW_Data/Turbidity_K490_OceanDataNASA/requested_files/", Tur_files, sep="")
Tur_list <- NULL
region <- ext(10, 40, -40, -25)
Tur_list <- rast(Tur_files)
crs(Tur_list) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

Tur_date <- substr(names(Tur_list), 12,19)
Tur_date_2 <- as.Date(paste(substr(Tur_date, 1,4),
                            substr(Tur_date, 5,6),
                            substr(Tur_date, 7,8), sep="-"), 
                      format="%Y-%m-%d")

n_tur_list <- nlyr(Tur_list)
tur_stack <- NULL

# Function is a loop, otherwise R crashes
for(i in c(1:n_tur_list)){
  Tur_sub <- Tur_list[[i]]
  Tur_sub[Tur_sub==255] <- NA
  tur_stack <- c(tur_stack, Tur_sub)
  
  progress(i, n_tur_list)
  
  rm(Tur_sub)
  gc()

}

tur_stack <- rast(tur_stack)

time(tur_stack) <- Tur_date_2

dts_tur = time(tur_stack) |> lubridate::ymd()
weekly_tur = paste(
  lubridate::isoweek(dts_tur),
  lubridate::year(dts_tur),
  sep = '-') |> factor()

print("Tur mean")
weekly_Tur_mean = terra::tapp(tur_stack, index = weekly_tur, fun = mean, na.rm=T)
Tur_mean <- terra::mean(weekly_Tur_mean, na.rm=T)
ext(Tur_mean) <- c(10, 40, -40, -25)
Tur_mean_new <- ResSampl_func(ext=e, ras=Tur_mean, res=new_res)
writeCDF(Tur_mean_new, "R Project/Variable_Maps/Tur_mean.nc", overwrite=TRUE, varname="Tur_mean", unit = "")
print("Tur min")
weekly_Tur_min = terra::tapp(tur_stack, index = weekly_tur, fun = min, na.rm=T)
Tur_min <- terra::mean(weekly_Tur_min, na.rm=T)
ext(Tur_min) <- c(10, 40, -40, -25)
Tur_min_new <- ResSampl_func(ext=e, ras=Tur_min, res=new_res)
writeCDF(Tur_min_new, "R Project/Variable_Maps/Tur_min.nc", overwrite=TRUE, varname="Tur_min", unit = "")
print("Tur max")
weekly_Tur_max = terra::tapp(tur_stack, index = weekly_tur, fun = max, na.rm=T)
Tur_max <- terra::mean(weekly_Tur_max, na.rm=T)
ext(Tur_max) <- c(10, 40, -40, -25)
Tur_max_new <- ResSampl_func(ext=e, ras=Tur_max, res=new_res)
writeCDF(Tur_max_new, "R Project/Variable_Maps/Tur_max.nc", overwrite=TRUE, varname="Tur_max", unit = "")
print("Tur range")
Tur_range = Tur_max - Tur_min
ext(Tur_range) <- c(10, 40, -40, -25)
Tur_range_new <- ResSampl_func(ext=e, ras= Tur_range, res=new_res)
writeCDF(Tur_range_new, "R Project/Variable_Maps/Tur_range.nc", overwrite=TRUE, varname="Tur_range", unit = "")
print("Tur var")
weekly_Tur_var = terra::tapp(tur_stack, index = weekly_tur, fun = sd, na.rm=T)^2
Tur_var <- terra::mean(weekly_Tur_var, na.rm=T)
ext(Tur_var) <- c(10, 40, -40, -25)
Tur_var_new <- ResSampl_func(ext=e, ras=Tur_var, res=new_res)
writeCDF(Tur_var_new, "R Project/Variable_Maps/Tur_var.nc", overwrite=TRUE, varname="Tur_var", unit = "")



tur_stack_l <- nlyr(weekly_Tur_mean)
tur_sobel_stack <- NULL

for(i in c(1:tur_stack_l)){
  w = matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  tur_sobel_st <- focal(weekly_Tur_mean[[i]], w = w)
  tur_vsobel_st <- focal(weekly_Tur_mean[[i]], w = t(w))
  comb_tur_sobel_st <- sqrt(tur_sobel_st^2 + tur_vsobel_st^2)
  tur_sobel_stack <- c(tur_sobel_stack,comb_tur_sobel_st)
  
  progress(i, tur_stack_l)
  
  rm(tur_sobel_st)
  rm(tur_vsobel_st)
  gc()
}

tur_sobel_stack <- rast(tur_sobel_stack)
tur_sobel_stack_sum <- app(tur_sobel_stack, fun=sum, na.rm=T)
ext(tur_sobel_stack_sum) <- c(10, 40, -40, -25)
tur_sobel_stack_sum <- ResSampl_func(ext=e, ras=tur_sobel_stack_sum, res=new_res)

writeCDF(tur_sobel_stack_sum, "R Project/Variable_Maps/Turbidity_sobel.nc", overwrite=TRUE, varname="Tur_sobel", unit = "mg m-3")




#### Combine datasets #### 
print("Start Combining Dataset process.")
comb <- c(s_mean_new, s_min_new, s_max_new, s_var_new, s_range_new,
          c_mean_new, c_min_new, c_max_new, c_var_new, c_range_new,
          h_mean_new, h_min_new, h_max_new, h_var_new, h_range_new, 
          npp_mean_new, npp_min_new, npp_max_new, npp_var_new, npp_range_new,
          curr_mean_new, curr_min_new, curr_max_new, curr_var_new, curr_range_new,
          Tur_mean_new, Tur_min_new,Tur_max_new, Tur_var_new, Tur_range_new,
          SST_sobel_stack_sum, Chl_sobel_stack_sum, Okubo_W_stack_mean_new,
          # hc_mean_new, hc_min_new, hc_max_new, hc_var_new, hc_range_new, 
          # vc_mean_new, vc_min_new, vc_max_new, vc_var_new, vc_range_new, 
          b_new, sl)
names(comb) <- c("SST_mean", "SST_min", "SST_max", "SST_var", "SST_range",
                 "Chl_mean", "Chl_min", "Chl_max", "Chl_var", "Chl_range",
                 "SSH_mean", "SSH_min", "SSH_max", "SSH_var", "SSH_range", 
                 "NPP_mean", "NPP_min", "NPP_max", "NPP_var", "NPP_range", 
                 "Curr_mean", "Curr_min", "Curr_max", "Curr_var", "Curr_range",
                 "Tur_mean", "Tur_min", "Tur_max", "Tur_var", "Tur_range",
                 "SST_Sobel", "Chl_Sobel",  "Okubo_Weiss_stack",
                 # "HC_mean", "HC_min", "HC_max", "HC_var", "HC_range", 
                 # "VC_mean", "VC_min", "VC_max", "VC_var", "VC_range", 
                 "Bathy", "Slope")

writeCDF(comb, "Combined_raster_data.nc", overwrite=TRUE, varname="comb", 
         longname="Pelagic datasets", unit="m")
print("Process finished - Data saved to Combined_raster_data.nc")
