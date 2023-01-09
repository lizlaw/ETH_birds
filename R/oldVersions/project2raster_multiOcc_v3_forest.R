# project to raster- forest - 10m  - cell by cell

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
  ) %>% as.matrix(.)

# cell-wise calculations -------------------------------------------------------
# set up results matrix
predictions <- matrix(data = NA, 
                      nrow = nrow(values_to_predict), ncol = 6, 
                      dimnames = list(values_to_predict[,"cellID"], 
                                      c("SR_site_all",
                                        "SR_site_tobs",
                                        "SR_site_fspp",
                                        "SR_site_mig",
                                        "SR_site_fnDiet",
                                        "SR_site_invDiet")))

# select draw from each chain
startSamples <- 1001
nSamples <- 80
nChains <- 3
mydraw <- startSamples:(startSamples+nSamples-1)
posterior_matrix <- rbind(samples$chain1[mydraw, ], samples$chain2[mydraw, ], samples$chain3[mydraw, ]) 

# set up functions for row-wise prediction
countSet <- round(seq(1, nrow(values_to_predict), length.out = 1000))

tictoc::tic("all rows")
for (i in 74552:nrow(values_to_predict)){
  # counter
  if(i %in% countSet) {
    cat(paste0(which(countSet==i)/10), " % ", date())
  }
  # extract psi for all species, samples
  out_psi <- foreach (k = 1:M, 
                      .combine = 'cbind' 
  ) %do% {
    out_logitpsi <- posterior_matrix[, paste0("lpsi", "[", k, "]")] +
      posterior_matrix[, paste0("betalpsi", "[", k, ", 1]")] * values_to_predict[i,1] + 
      posterior_matrix[, paste0("betalpsi", "[", k, ", 2]")] * values_to_predict[i,2] +
      posterior_matrix[, paste0("betalpsi", "[", k, ", 3]")] * values_to_predict[i,3] 
    
    #out_psi[, k,] <- 
    exp(out_logitpsi)/(exp(out_logitpsi) + 1)
  }
  
  # derived quantities
  predictions[i,"SR_site_all"] <- rowSums(out_psi) %>% mean()             # Number of occurring species at each site (sum of psi)
  predictions[i,"SR_site_tobs"] <- rowSums(out_psi[,1:Nobs]) %>% mean()     # Number of observed species occurring at each site (i.e. not including the non-observed species)
  predictions[i,"SR_site_fspp"] <- rowSums(out_psi[,which(fspp==TRUE)]) %>% mean()         # forest dependent species
  predictions[i,"SR_site_mig"] <- rowSums(out_psi[,which(mig==TRUE)]) %>% mean()           # migrants
  predictions[i,"SR_site_fnDiet"] <- rowSums(out_psi[,which(fnDiet==TRUE)]) %>% mean()        # fruit/nectar species
  predictions[i,"SR_site_invDiet"] <- rowSums(out_psi[,which(invDiet==TRUE)]) %>% mean()      # invertebrate species
  
  rm(out_psi)
  } # end for cell i loop 

tictoc::toc()   

# create convert values to raster function -----------------------------------
insert2raster <- function(x, 
                          name=deparse(substitute(x)), 
                          template = forest_stack[[4]],
                          cells = values_to_predict[,"cellID"]){
  out <- template
  out[cells] <- x
  names(out) <- name
  return(out)
}

# convert to raster stack
SR_site_forest_stack <- rast(list(
  SR_site_all_raster = insert2raster(predictions[,"SR_site_all"], name = "SR_site_all"),
  SR_site_fspp_raster = insert2raster(predictions[,"SR_site_fspp"], name = "SR_site_fspp"),
  SR_site_mig_raster = insert2raster(predictions[,"SR_site_mig"], name = "SR_site_mig"),
  SR_site_fnDiet_raster = insert2raster(predictions[,"SR_site_fnDiet"], name = "SR_site_fnDiet"),
  SR_site_invDiet_raster = insert2raster(predictions[,"SR_site_invDiet"], name = "SR_site_invDiet")
))
  
  # save raster stack to file
  writeRaster(SR_site_forest_stack, 
              paste0(results_folder, 
                     version_folder, 
                     "BirdOccMod_dOcc_SRsite_forestRaster_10m_", 
                     mydate, 
                     ".tif"), 
              overwrite=TRUE)

# session info -----------------------------------------------------------------
date()
sessionInfo() 


