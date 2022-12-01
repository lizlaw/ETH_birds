# project to raster- forest - 10m  - rep by rep

# libraries --------------------------------------------------------------------
library(tidyverse)
library(terra)
library(bestNormalize)
library(foreach)
library(doParallel)

# data -------------------------------------------------------------------------
wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds"
results_folder <- paste0(wd, "/Results/")  
version_folder <- "v9/"
mydate <- 20221017
samples <- readRDS(paste0(results_folder, version_folder, "BirdOccMod_dOcc_samples_forest_", mydate, ".RDS"))
data_input <- readRDS(paste0(wd, "/Data/Forest_nimbleOccData_v6_", mydate, ".RDS"))
list2env(data_input[[1]], envir = .GlobalEnv); list2env(data_input[[2]], envir = .GlobalEnv); rm(data_input)
OccVarsForest_scaleModels <- readRDS("~/Desktop/LeuphanaProject/BirdModelling/ETH_birds/Data/nimbleModel_multiOcc_v6_OccVarsFrst_scaleModels.RDS")

# out_temp_psi_folder  <- paste0(results_folder, version_folder,"temp_psi_rasters/")
# out_temp_rep_folder  <- paste0(results_folder, version_folder,"temp_rep_rasters/")
out_temp_psi_folder  <- "/Volumes/WCC_2021/temp_psi_rasters/"
out_temp_rep_folder  <- "/Volumes/WCC_2021/temp_rep_rasters/"
out_final_folder <- "forest_summary_rasters/"

# if raster not pre-created ----------------------------------------------------
rast_folder <- "/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Baseline/"

makeStack <-  FALSE
if (makeStack == TRUE){
  
  elevation <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Elevation/Elevation_SA_RS_10.tif")
  frdis <- rast(paste0(rast_folder, "FRdis baseline/FRdis.tif"))
  frtype <- rast(paste0(rast_folder, "FRType baseline/FR_FRType.tif"))
  
  forest_shp <- vect(paste0(rast_folder, "Shapefile/Forestland_Baseline.shp"))
  
  # create template, resample all rasters, tile, and save in a stack -------------
  forest_template <- rasterize(forest_shp, frtype)
  elevation <- resample(elevation, forest_template, method = "bilinear") %>% mask(., forest_template) 
  frdis <- resample(frdis, forest_template, method = "bilinear") %>% mask(., forest_template) 
  frtype <- resample(frtype, forest_template, method = "mode") %>% mask(., forest_template) 
  
  forest_stack <- rast(list(elevation, frdis, frtype, forest_template))
  names(forest_stack) <- c("elevation", "frdis", "foresttype", "forest_template")
  
  writeRaster(forest_stack, paste0(rast_folder, "forest_stack_10m.tif"))
  rm(forest_shp, forest_template, elevation, frdis, forest_stack)
  
} # end make stack

# load the raster stack
forest_stack <- rast(paste0(rast_folder, "forest_stack_10m.tif"))

# set up for projection  ------------------------------------------------------

make_pred_stack <- FALSE
if(make_pred_stack == TRUE){
  # extract then transform the raster variables (elevation, frdis, frtype)
  values_to_predict <- terra::values(forest_stack[[c(1,2,3)]], na.rm=FALSE, dataframe = TRUE) 
  values_to_predict$cellID <- 1:nrow(values_to_predict)
  values_to_predict <- values_to_predict[complete.cases(values_to_predict),]
  
  # transform
  values_to_predict <- values_to_predict %>% 
    mutate(elevation = predict(OccVarsForest_scaleModels$elevation, newdata = elevation),
           frdis = predict(OccVarsForest_scaleModels$fr_dis,
                           newdata = frdis),
           foresttype = predict(OccVarsForest_scaleModels$forest_type,
                                newdata = foresttype)
    ) 
  
  # convert back to raster
  pred_stack <- subset(forest_stack, 1:3)
  pred_stack$elevation[values_to_predict$cellID] <- values_to_predict$elevation
  pred_stack$frdis[values_to_predict$cellID] <- values_to_predict$frdis
  pred_stack$foresttype[values_to_predict$cellID] <- values_to_predict$foresttype
  
  # save and reload raster stack
  writeRaster(pred_stack, paste0(results_folder, version_folder, "forest_pred_stack.tif"), overwrite = TRUE)
  rm(values_to_predict)
}

pred_stack <- rast(paste0(results_folder, version_folder, "forest_pred_stack.tif"))

# rep-wise calculations -------------------------------------------------------

# select draw from each chain
startSamples <- 1001
nSamples <- 80
nChains <- 3
mydraw <- startSamples:(startSamples+nSamples-1)
posterior_matrix <- rbind(samples$chain1[mydraw, ], samples$chain2[mydraw, ], samples$chain3[mydraw, ]) 
nRep <- dim(posterior_matrix)[1]

# run for each rep at a time:

foreach (i = 101:150, #1:nRep,
         .packages = c("terra", "foreach")
         ) %do% {

  # for each species, pull out the intercepts, the betalpsi 1, 2, and 3
  foreach (k = 1:(rareSpp)) %do% {
    cat("Starting: spp", k, "rep", i, date(), "\n")
    out_logitpsi <- posterior_matrix[i, paste0("lpsi", "[", k, "]")]  +
      posterior_matrix[i, paste0("betalpsi", "[", k, ", 1]")] * pred_stack[[1]] + 
      posterior_matrix[i, paste0("betalpsi", "[", k, ", 2]")] * pred_stack[[2]] +
      posterior_matrix[i, paste0("betalpsi", "[", k, ", 3]")] * pred_stack[[3]] 
    
    #out_psi[, k,] <- 
    out_psi <- exp(out_logitpsi)/(exp(out_logitpsi) + 1)
    writeRaster(out_psi, 
                paste0(out_temp_psi_folder, 
                       "psi_sp", k, ".tif"),
                overwrite=TRUE)
    rm(out_logitpsi, out_psi)
    file.remove(list.files(tempdir(), full.names = T, pattern = "^spat"))
  }
  # all the rare species share parameters so can just read-write these
  foreach (k = (rareSpp+1):M) %do% {
    cat("Starting: spp", k, "rep", i, date(), "\n")
    out_psi <- rast(paste0(out_temp_psi_folder, 
                       "psi_sp", k-1,  ".tif"))
    writeRaster(out_psi, 
                paste0(out_temp_psi_folder, 
                       "psi_sp", k, ".tif"),
                overwrite=TRUE)
  }
  
  psi_stack <- rast(paste0(out_temp_psi_folder,
                           "psi_sp", 1:M, ".tif"))
  
  # derived quantities
  cat("Starting: SR_site_all rep", i, date(), "\n")
  SR_site_all <- sum(psi_stack)             # Number of occurring species at each site (sum of psi)
  writeRaster(SR_site_all, 
              paste0(out_temp_rep_folder, 
                     "SR_site_all", "_rep", i, ".tif"),
              overwrite=TRUE)
  rm(SR_site_all)
  cat("Starting: SR_site_fspp rep", i, date(), "\n")
  SR_site_fspp <- sum(psi_stack[[which(fspp==TRUE)]])       # forest dependent species
  writeRaster(SR_site_fspp, 
              paste0(out_temp_rep_folder, 
                     "SR_site_fspp", "_rep", i, ".tif"),
              overwrite=TRUE)
  rm(SR_site_fspp)
  cat("Starting: SR_site_mig rep", i, date(), "\n")
  SR_site_mig <- sum(psi_stack[[which(mig==TRUE)]])          # migrants
  writeRaster(SR_site_mig, 
              paste0(out_temp_rep_folder, 
                     "SR_site_mig", "_rep", i, ".tif"),
              overwrite=TRUE)
  rm(SR_site_mig)
  cat("Starting: SR_site_fnDiet rep", i, date(), "\n")
  SR_site_fnDiet <- sum(psi_stack[[which(fnDiet==TRUE)]])   # fruit/nectar species
  writeRaster(SR_site_fnDiet, 
              paste0(out_temp_rep_folder, 
                     "SR_site_fnDiet", "_rep", i, ".tif"),
              overwrite=TRUE)
  rm(SR_site_fnDiet)
  cat("Starting: SR_site_invDiet rep", i, date(), "\n")
  SR_site_invDiet <- sum(psi_stack[[which(invDiet==TRUE)]]) # invertebrate species
  writeRaster(SR_site_invDiet, 
              paste0(out_temp_rep_folder, 
                     "SR_site_invDiet", "_rep", i, ".tif"),
              overwrite=TRUE)
  rm(SR_site_invDiet)
  
  file.remove(list.files(tempdir(), full.names = T, pattern = "^spat"))
  
} # end for rep i loop 

# Starting: spp 1 rep 2 Fri Nov 18 13:39:10 2022 
# 

# stopCluster()

# then need to combine the reps as required

# session info -----------------------------------------------------------------
date()
sessionInfo() 
