# Forest - make raster stacks for projection
# includes resampling and transformations

library(tidyverse)
library(terra)
library(bestNormalize)

wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds"
results_folder <- paste0(wd, "/Results/")  
version_folder <- "v11/"
stack_folder <- paste0("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_stacks/")
elevation <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Elevation/Elevation_SA_RS_10.tif")
hli <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/HLI/hli.tif")
#slope <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Slope/Slope_dgr.tif")
OccVarsforest_scaleModels <- readRDS("~/Desktop/LeuphanaProject/BirdModelling/ETH_birds/Data/nimbleModel_multiOcc_v7_OccVarsfrst_scaleModels.RDS")

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
  forest_85 <- rast(paste0(rast_folder, "FRtype_",scenario,".tif"))
  fr_dis <- rast(paste0(rast_folder, "FRdis_",scenario,".tif"))

  # resample all rasters, tile, and save in a stack -------------
  # order: "elevation1" "elevation2" "forest_85"  "fr_dis" "hli" 
  elevation <- resample(elevation, template, method = "bilinear") %>% mask(., template)
  fr_dis <- resample(fr_dis, template, method = "bilinear") %>% mask(., template)
  forest_85 <- resample(forest_85, template, method = "mode") %>% mask(., template)
  hli <- resample(hli, template, method = "bilinear") %>% mask(., template)
  
  forest_stack <- rast(list(elevation=elevation, 
                            forest_85=forest_85, 
                            fr_dis=fr_dis, 
                            hli=hli))
  writeRaster(forest_stack, paste0(stack_folder, "forest_stack_10m_",scenario,".tif"), overwrite=TRUE)
  forest_stack <- rast(paste0(stack_folder, "forest_stack_10m_",scenario,".tif"))
  
  # resample to coarser resolution
  forest_stack_100m <- rast(list(
    elevation = elevation %>% aggregate(fact=10, fun=mean, na.rm=TRUE), 
    forest_85 = forest_85 %>% aggregate(fact=10, fun=median, na.rm=TRUE), 
    fr_dis = fr_dis %>% aggregate(fact=10, fun=mean, na.rm=TRUE), 
    hli = hli %>% aggregate(fact=10, fun=mean, na.rm=TRUE)
  ))
  writeRaster(forest_stack_100m, paste0(stack_folder, "forest_stack_100m_",scenario,".tif"), overwrite=TRUE)
  forest_stack_100m <- rast(paste0(stack_folder, "forest_stack_100m_",scenario,".tif"))

  # Transform stacks -------------------------------------------------------------
  # loop over resolutions 10m and 100m
  for (i in 1:2){
    
    if(i==1){
      mystack <- forest_stack
      savename <- "forest_pred_stack_10m"
    }
    
    if(i==2){
      mystack <- forest_stack_100m
      savename <- "forest_pred_stack_100m"
    }

    # extract then transform the raster variables (elevation, fr_dis, sidi1ha)
    values_to_predict <- terra::values(mystack, na.rm=FALSE, dataframe = TRUE) 
    values_to_predict$cellID <- 1:nrow(values_to_predict)
    values_to_predict <- values_to_predict[complete.cases(values_to_predict),]
    
    # transform "elevation1" "elevation2" "forest_85"  "fr_dis" "hli" 
    elevPoly <- predict(OccVarsforest_scaleModels$elevation, newdata = values_to_predict$elevation)
    values_to_predict <- values_to_predict %>% 
      mutate(elevation1 = predict(OccVarsforest_scaleModels$elevation1, 
                                  newdata = elevPoly[,1]),
             elevation2 = predict(OccVarsforest_scaleModels$elevation2, 
                                  newdata = elevPoly[,2]),
             forest_85 = predict(OccVarsforest_scaleModels$forest_85,
                                 newdata = forest_85),
             fr_dis = predict(OccVarsforest_scaleModels$fr_dis,
                             newdata = fr_dis),
             hli = predict(OccVarsforest_scaleModels$hli,
                             newdata = hli)
      ) 
    
    # convert back to raster
    pred_stack <- subset(mystack, c(1,1:4))
    names(pred_stack) <- c("elevation1", "elevation2", "forest_85", "fr_dis", "hli")
    pred_stack$elevation1[values_to_predict$cellID] <- values_to_predict$elevation1
    pred_stack$elevation2[values_to_predict$cellID] <- values_to_predict$elevation2
    pred_stack$forest_85[values_to_predict$cellID] <- values_to_predict$forest_85
    pred_stack$fr_dis[values_to_predict$cellID] <- values_to_predict$fr_dis
    pred_stack$hli[values_to_predict$cellID] <- values_to_predict$hli
    
    # save and reload raster stack
    writeRaster(pred_stack, paste0(results_folder, version_folder, 
                                   scenario, "/",
                                   savename, "_", scenario, ".tif"), 
                overwrite = TRUE)
    
  } # end resolution loop
} # end scenario loop.

