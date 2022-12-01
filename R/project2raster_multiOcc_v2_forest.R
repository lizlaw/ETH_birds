# project to raster- forest - 10m  - using tiles 

# libraries --------------------------------------------------------------------
library(tidyverse)
library(terra)
library(bestNormalize)
library(foreach)
library(doParallel)
library(tictoc)  # timer for code development

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
makeTiles <-  FALSE
if (makeTiles == TRUE){
  
  rast_folder <- "/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Baseline/"
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
  
  tile_template <- aggregate(forest_template, fact=500) ## fact 500 = 270 tiles, fact 1000 = 72 tiles
  tile_template[] <- 1
  # ncell(tile_template)
  forest_tile_stack <-
    makeTiles(forest_stack, 
              tile_template, 
              filename = "/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Baseline/forest_stack_10m_tile_.tif",
              overwrite = TRUE)
  
  rm(forest_shp, forest_template, elevation, frdis, forest_stack, tile_template)
  
} else { 

  forest_tile_stack <- list.files(path = "/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Baseline/",
                                pattern = "forest_stack_10m_tile_",
                                full.names = TRUE)
} # end make tiles 

# set up for projection 1 ------------------------------------------------------
startSamples <- 1001
nSamples <- 70
nChains <- 3
acomb <- function(...) abind::abind(..., along=3)
# select draw from each chain
mydraw <- startSamples:(startSamples+nSamples)
posterior_matrix <- rbind(samples$chain1[mydraw, ], samples$chain2[mydraw, ], samples$chain3[mydraw, ]) 

# for each of the tiles --------------------------------------------------------

for (tile in 1:3){#:length(forest_tile_stack)){
  # tile <- 4
  forest_stack <- rast(forest_tile_stack[tile])
  
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
  
  if (nrow(values_to_predict)==0){
    SR_site_all <- rep(NA, nrow(values_to_predict))             # Number of occurring species at each site (sum of psi)
    SR_site_tobs <- rep(NA, nrow(values_to_predict))      # Number of observed species occurring at each site (i.e. not including the non-observed species)
    SR_site_fspp <- rep(NA, nrow(values_to_predict))          # forest dependent species
    SR_site_mig <- rep(NA, nrow(values_to_predict))            # migrants
    SR_site_fnDiet <- rep(NA, nrow(values_to_predict))        # fruit/nectar species
    SR_site_invDiet <- rep(NA, nrow(values_to_predict))      # invertebrate species
    
  } else {
   
    # set up for the projection 2 
    nCells <- nrow(values_to_predict)
    array_NA <- array(data = NA, dim = c(nCells, nSamples*nChains, M))  # all species
    myMatrix <- function(x){matrix(data = x, nrow = nCells, ncol = nSamples*nChains, byrow = FALSE)}
    
    # for each species, pull out the intercepts, the betalpsi 1, 2, and 3
    tic("method 1")
    cl <- makeCluster(3)
    registerDoParallel(cl)
    
    out_psi <- foreach (k = 1:M, 
                        .combine = 'acomb', .multicombine = TRUE, .maxcombine = M
    ) %dopar% {
      
      out_logitpsi <- myMatrix(posterior_matrix[, paste0("lpsi", "[", k, "]")])  +
        myMatrix(posterior_matrix[, paste0("betalpsi", "[", k, ", 1]")]) * values_to_predict[,1] + 
        myMatrix(posterior_matrix[, paste0("betalpsi", "[", k, ", 2]")]) * values_to_predict[,2] +
        myMatrix(posterior_matrix[, paste0("betalpsi", "[", k, ", 3]")]) * values_to_predict[,3] 
      
      #out_psi[, k,] <- 
      exp(out_logitpsi)/(exp(out_logitpsi) + 1)
    }
    
    stopCluster(cl)
    toc()
    
    # save psi results 
    # saveRDS(out_psi, 
    #         paste0(results_folder, 
    #                version_folder, 
    #                "BirdOccMod_dOcc_psi_forestRaster_10m_", 
    #                tile, "_",
    #                mydate, 
    #                ".RDS"))
    
    # derived quantities
    SR_site_all <- apply(out_psi, c(1,2), sum) %>% apply(., 1, mean)             # Number of occurring species at each site (sum of psi)
    SR_site_tobs <- apply(out_psi[,,1:Nobs],  c(1,2), sum) %>% apply(., 1, mean)     # Number of observed species occurring at each site (i.e. not including the non-observed species)
    SR_site_fspp <- apply(out_psi[,,which(fspp==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)          # forest dependent species
    SR_site_mig <- apply(out_psi[,,which(mig==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)           # migrants
    SR_site_fnDiet <- apply(out_psi[,,which(fnDiet==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)       # fruit/nectar species
    SR_site_invDiet <- apply(out_psi[,,which(invDiet==TRUE)],  c(1,2), sum) %>% apply(., 1, mean)     # invertebrate species
  
    rm(out_psi)
  } # end if else
  
  # create convert values to raster function
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
    SR_site_all_raster = insert2raster(SR_site_all),
    SR_site_fspp_raster = insert2raster(SR_site_fspp),
    SR_site_mig_raster = insert2raster(SR_site_mig),
    SR_site_fnDiet_raster = insert2raster(SR_site_fnDiet),
    SR_site_invDiet_raster = insert2raster(SR_site_invDiet)
    ))
  
  # save raster stack to file
  writeRaster(SR_site_forest_stack, 
              paste0(results_folder, 
                     version_folder, 
                     "BirdOccMod_dOcc_SRsite_forestRaster_10m_", 
                     tile, "_",
                     mydate, 
                     ".tif"), 
              overwrite=TRUE)
  
  rm(SR_site_forest_stack, insert2raster,
     SR_site_all, SR_site_fspp, SR_site_mig, SR_site_fnDiet, SR_site_invDiet)
  
} # end tile loop
 

# merge all the tiles together and save ----------------------------------------

SR_site_forest_stack_merge <- merge(
  sprc(
    map(paste0(results_folder, 
               version_folder, 
               "BirdOccMod_dOcc_SRsite_forestRaster_10m_", 
               1:length(forest_tile_stack), "_",
               mydate, 
               ".tif"),
        rast)
    )
  )

# session info -----------------------------------------------------------------
date()
sessionInfo() 


  