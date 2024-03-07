library(terra)
library(tidyverse)
library(ggplot2)
library(gganimate)
library(tidyterra)
library(gifski)
library(av)

# Set WD to where you want to have your SST Animation
setwd("~/Documents/SANBI_Pelagic_Ecosystems_PostDoc/")
# Make sure there is a folder within this WD where you can place the individual slides (e.g. "SST_Animated")

# Read in the dataset in RasterStack format 
# Set to a file within your WD
s = terra::rast("RAW_Data/aphiwe_layers/SST/Data/METOFFICE-GLO-SST-L4-REP-OBS-SST_multi-vars_10.02E-39.97E_39.97S-25.02S_2011-01-01-2021-12-31.nc", subds = "analysed_sst")

# Subset to year required 
s_sub <- subset(s, time(s) >= as.Date("2021-01-01"))

# Data is in kelvin, thus needs to be converted to °C
s_sub <- s_sub - 273.15 

# Confirm NAs are NAs
s_sub[is.na(s_sub)] <- NA

# Fix the extent of the data to the required frame
ext(s_sub) <- c(10, 40, -40, -25)

# Specify the CRS
crs(s_sub) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 


# Change names for times so date is used for slide-name
names(s_sub) <- time(s_sub)
s_names <- names(s_sub)

# Determine extreme values (min, max) for ggplot label 
s_min_val <- floor(min(values(s_sub), na.rm = T))
s_max_val <- ceiling(max(values(s_sub), na.rm = T))


# Loop through date name, create graph for each slide, and safe it to your computer 
for(i in seq_len(length(s_names))){
  
  s_layer_name <- s_names[[i]]
  
  s_layer <- s_sub |>
    dplyr::select(
      dplyr::all_of(
        s_layer_name
      )
    )
  
# This makes the graph which will become the individual frames of animation. Use standard ggplot adjustments where necessary.  
  p <- ggplot() + 
    geom_spatraster(data=s_layer, aes()) + 
    scale_fill_distiller(palette = "Spectral", 
                         direction = -1, 
                         limits=c(s_min_val, s_max_val)) + 
    labs(subtitle = s_layer_name,
         fill="SST") + 
    theme_classic() +
    theme(line = element_blank(),
          plot.subtitle = element_text(hjust=0.5))
    
# Create list of file paths to the individual slides. 
# Change "SST_Animated" to the name of the folder where you want the individual slides saved.
  map_name <- file.path(
    getwd(), paste0("SST_Animated/",
                    s_layer_name,
      ".png"
    )
  )
  
# Save the ggplot maps to folder above (as png). You can adjust the size of the image as required.  
  ggsave(
    map_name, p,
    height = 5, width = 7.5 ,
    units = "in", bg = "white"
  )
  
}

# This calls out the paths and file-names of the individual maps. 
map_files <- file.path(
  getwd(), paste0("SST_Animated/",
    s_names,
    ".png"
  )
)

# This takes the individual png files in the folder and turns them into a gif.
gifski::gifski(
  map_files,
  loop = T, 
  delay = .1, # the speed of transition between one frame to another.
  width = 700,
  height = 450,
  gif_file = "SA_SST_2021.gif" 
)

# This turns the pngs into a video (mp4)
av::av_encode_video(map_files, framerate = 20,
                    output = 'test.mp4')
