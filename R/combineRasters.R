# combine both rasters

library(tidyverse)
library(terra)

wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds"
results_folder <- paste0(wd, "/Results/")  
version_folder <- "v9/"
mydate <- 20221017

farmStack <- rast(paste0(results_folder, version_folder, 
                         "BirdOccMod_dOcc_SRsite_farmRaster_300m", mydate, ".tif"))
forestStack <- rast(paste0(results_folder, version_folder, 
                           "BirdOccMod_dOcc_SRsite_forestRaster_300m", mydate, ".tif"))

mergeStack <- merge(forestStack, farmStack) # merge will use the x in the overlaps. ours is the same. 

# `mosaic` will use a function. `vrt` is more efficient. 

plot(mergeStack[[1]])
plot(mergeStack[[2:5]])

(merge(forestStack[[1]], farmStack[[1]]) - merge(farmStack[[1]], forestStack[[1]]) +1 ) %>% plot(col="black")