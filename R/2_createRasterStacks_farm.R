# Farm - make raster stacks for projection
# includes resampling and transformations

library(tidyverse)
library(terra)
library(bestNormalize)

wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds"
results_folder <- paste0(wd, "/Results/")  
version_folder <- "v11/"
stack_folder <- paste0("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_stacks/")
elevation <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Elevation/Elevation_SA_RS_10.tif")
slope <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Slope/Slope_dgr.tif")
OccVarsFarm_scaleModels <- readRDS("~/Desktop/LeuphanaProject/BirdModelling/ETH_birds/Data/nimbleModel_multiOcc_v7_OccVarsFarm_scaleModels.RDS")

# create template if required, or load it ------------------------------------
if(FALSE){
  farm_shp <- vect(paste0(rast_folder, "Shapefile/Farmland_",scenario,".shp"))
  farm_template <- rasterize(farm_shp, sidi1ha)  # default center point only ok here, as is 10m raster
  forest_shp <- vect(paste0(rast_folder, "Shapefile/Forestland_Baseline.shp"))
  forest_template <- rasterize(forest_shp, sidi1ha)
  template <- merge(farm_template, 2*forest_template)
  writeRaster(template, "/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/ETH_template_10m.tif", overwrite=TRUE)
  # also make proportion template for projections at 100m
  farm_template_10 <- farm_template %>% aggregate(fact=10, fun=function(x){sum(x, na.rm=TRUE)/100})
  forest_template_10 <-  forest_template %>% aggregate(fact=10, fun=function(x){sum(x, na.rm=TRUE)/100})
  writeRaster(rast(list(farm_proportion = farm_template_10, 
                        forest_proportion = forest_template_10)),
              "/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/ETH_farmForestProportions_100m.tif", overwrite=TRUE)
  rm(farm_shp, forest_shp, farm_template, forest_template, farm_template_10, forest_template_10)
}
template <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/ETH_template_10m.tif")

# loop over scenarios ----------------------------------------------------------
#scenario <- "Baseline"  # Baseline, BR, CC, CI, FF ## all done
for (scenario in c("Baseline", "BR", "CC", "FF", "CI")){

  rast_folder <- paste0("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/",scenario,"/")
  fl_dis <- rast(paste0(rast_folder, "fldis_",scenario,".tif"))
  sidi1ha <- rast(paste0(rast_folder, "FLsidi1ha_",scenario,".tif"))
  farm_85 <- rast(paste0(rast_folder, "FLtype_",scenario,".tif"))
  
  # resample all rasters, tile, and save in a stack -------------
  # order: # 1"elevation1" 2"elevation2" 3"farm_85" 4"fl_dis" 5 "sidi1ha" 6"slope
  elevation <- resample(elevation, template, method = "bilinear") %>% mask(., template)
  farm_85 <- resample(farm_85, template, method = "mode") %>% mask(., template)
  fl_dis <- resample(fl_dis, template, method = "bilinear") %>% mask(., template)
  sidi1ha <- resample(sidi1ha, template, method = "mode") %>% mask(., template)
  slope <- resample(slope, template, method = "bilinear") %>% mask(., template) 
  
  farm_stack <- rast(list(elevation, farm_85, fl_dis, sidi1ha, slope))
  writeRaster(farm_stack, paste0(stack_folder, "farm_stack_10m_",scenario,".tif"), overwrite=TRUE)
  farm_stack <- rast(paste0(stack_folder, "farm_stack_10m_",scenario,".tif"))
  
# resample to coarser resolution
  farm_stack_100m <- rast(list(
    elevation = farm_stack$elevation %>% aggregate(fact=10, fun=mean, na.rm=TRUE), 
    farm_85 = farm_stack$farm_85 %>% aggregate(fact=10, fun=median, na.rm=TRUE),
    fl_dis = farm_stack$fl_dis %>% aggregate(fact=10, fun=mean, na.rm=TRUE), 
    sidi1ha = farm_stack$sidi1ha %>% aggregate(fact=10, fun=mean, na.rm=TRUE), 
    slope = farm_stack$slope %>% aggregate(fact=10, fun=mean, na.rm=TRUE)
    ))
  writeRaster(farm_stack_100m, paste0(stack_folder, "farm_stack_100m_",scenario,".tif"), overwrite=TRUE)
  farm_stack_100m <- rast(paste0(stack_folder, "farm_stack_100m_",scenario,".tif"))
  
# Transform stacks -------------------------------------------------------------
# loop over resolutions 10m and 100m
  for (i in 1:2){
    
    if(i==1){
      mystack <- farm_stack
      savename <- "farm_pred_stack_10m"
    }

    if(i==2){
      mystack <- farm_stack_100m
      savename <- "farm_pred_stack_100m"
    }

    # extract then transform the raster variables (elevation, fl_dis, sidi1ha)
    values_to_predict <- terra::values(mystack, na.rm=FALSE, dataframe = TRUE) 
    values_to_predict$cellID <- 1:nrow(values_to_predict)
    values_to_predict <- values_to_predict[complete.cases(values_to_predict),]
    
    # transform
    elevPoly <- predict(OccVarsFarm_scaleModels$elevation, newdata = values_to_predict$elevation)
    values_to_predict <- values_to_predict %>% 
      select(-elevation) %>% 
      mutate(elevation1 = predict(OccVarsFarm_scaleModels$elevation1, 
                                 newdata = elevPoly[,1]),
             elevation2 = predict(OccVarsFarm_scaleModels$elevation2, 
                                 newdata = elevPoly[,2]),
             fl_dis = predict(OccVarsFarm_scaleModels$fl_dis,
                             newdata = fl_dis),
             sidi1ha = predict(OccVarsFarm_scaleModels$sidi1ha,
                               newdata = sidi1ha),
             slope = predict(OccVarsFarm_scaleModels$slope,
                             newdata = slope),
             farm_85 = predict(OccVarsFarm_scaleModels$farm_85,
                                 newdata = farm_85)
      ) 
    
    # convert back to raster
    pred_stack <- subset(mystack, c(1,1:5))
    names(pred_stack) <- c("elevation1", "elevation2", "farm_85", "fl_dis", "sidi1ha", "slope")
    pred_stack$elevation1[values_to_predict$cellID] <- values_to_predict$elevation1
    pred_stack$elevation2[values_to_predict$cellID] <- values_to_predict$elevation2
    pred_stack$farm_85[values_to_predict$cellID] <- values_to_predict$farm_85
    pred_stack$fl_dis[values_to_predict$cellID] <- values_to_predict$fl_dis
    pred_stack$sidi1ha[values_to_predict$cellID] <- values_to_predict$sidi1ha
    pred_stack$slope[values_to_predict$cellID] <- values_to_predict$slope

    # save raster stack
    writeRaster(pred_stack, paste0(results_folder, version_folder, 
                                   scenario, "/",
                                   savename, "_", scenario, ".tif"), 
                overwrite = TRUE)
    
  } #end resolution loop
} #end scenario loop  