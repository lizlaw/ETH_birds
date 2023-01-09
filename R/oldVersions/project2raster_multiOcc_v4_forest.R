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

# set up functions for sets of row-wise prediction
startPredn <- 1
endPredn <- nrow(values_to_predict)
byPredn <- 5000
countSet <- c(seq(from=startPredn, to=endPredn, by=byPredn), endPredn) %>% unique()
acomb <- function(...) abind::abind(..., along=3)

tictoc::tic("sets of row-wise prediction")

cl <- makeCluster(2)
registerDoParallel(cl)

for (i in 1:length(countSet)){
  # stop if the end
  if(i==length(countSet)) stop("Not an error, just the end :)")
  # else, print counter
  cat("Starting:", round(i/length(countSet+1)*100, 3), "%", date(), "\n")
  
  # set up 
  iCells <- countSet[i]:(countSet[i+1]-ifelse(i==(length(countSet)-1),0,1))
  nCells <- length(iCells)
  array_NA <- array(data = NA, dim = c(nCells, nSamples*nChains, M))  # all species
  myMatrix <- function(x){matrix(data = x, nrow = nCells, ncol = nSamples*nChains, byrow = FALSE)}
  
  # for each species, pull out the intercepts, the betalpsi 1, 2, and 3
  out_psi <- foreach (k = 1:M, 
                      .combine = 'acomb', .multicombine = TRUE, .maxcombine = M
  ) %dopar% {
    
    out_logitpsi <- myMatrix(posterior_matrix[, paste0("lpsi", "[", k, "]")])  +
      myMatrix(posterior_matrix[, paste0("betalpsi", "[", k, ", 1]")]) * values_to_predict[iCells,1] + 
      myMatrix(posterior_matrix[, paste0("betalpsi", "[", k, ", 2]")]) * values_to_predict[iCells,2] +
      myMatrix(posterior_matrix[, paste0("betalpsi", "[", k, ", 3]")]) * values_to_predict[iCells,3] 
    
    #out_psi[, k,] <- 
    exp(out_logitpsi)/(exp(out_logitpsi) + 1)
  }
  
  # derived quantities
  predictions[iCells,"SR_site_all"] <- apply(out_psi, c(1,2), sum) %>% apply(., 1, mean)             # Number of occurring species at each site (sum of psi)
  predictions[iCells,"SR_site_tobs"] <- apply(out_psi[,,1:Nobs],  c(1,2), sum) %>% apply(., 1, mean)     # Number of observed species occurring at each site (i.e. not including the non-observed species)
  predictions[iCells,"SR_site_fspp"] <- apply(out_psi[,,which(fspp==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)          # forest dependent species
  predictions[iCells,"SR_site_mig"]<- apply(out_psi[,,which(mig==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)           # migrants
  predictions[iCells,"SR_site_fnDiet"] <- apply(out_psi[,,which(fnDiet==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)       # fruit/nectar species
  predictions[iCells,"SR_site_invDiet"] <- apply(out_psi[,,which(invDiet==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)     # invertebrate species
  
  rm(out_psi)
  
  write_csv(as.data.frame(predictions[iCells,]), 
            paste0(results_folder, version_folder, "forest_predictions/forest_predictions_", i, ".csv"))
} # end for cell i loop 

stopCluster()
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


