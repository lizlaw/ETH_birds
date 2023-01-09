# project to raster - forest - rep by rep

# libraries --------------------------------------------------------------------
library(tidyverse)
library(terra)
library(bestNormalize)
library(foreach)
library(doParallel)
library(bayestestR)

# data -------------------------------------------------------------------------
wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds"
results_folder <- paste0(wd, "/Results/")  
version_folder <- "v10/"
mydate <- 20221207
samples <- readRDS(paste0(results_folder, version_folder, "BirdOccMod_dOcc_samples_forest_", mydate, "_seed101.RDS"))
data_input <- readRDS(paste0(wd, "/Data/Forest_nimbleOccData_v6_", mydate, ".RDS"))
list2env(data_input[[1]], envir = .GlobalEnv); list2env(data_input[[2]], envir = .GlobalEnv); rm(data_input)

out_temp_psi_folder  <- "/Volumes/WCC_2021/temp_psi_rasters/"
out_temp_rep_folder  <- "/Volumes/WCC_2021/temp_rep_rasters/"

# select draw from each chain
nRep <- 252 # temporary, will be rounded up if not divisible by 3
mydraw <- seq(1, dim(samples$chain1)[1], length.out=ceiling(nRep/3))
posterior_matrix <- rbind(samples$chain1[mydraw, ], samples$chain2[mydraw, ], samples$chain3[mydraw, ]) 
nRep <- dim(posterior_matrix)[1]

rm(samples) # dont need to hold this in memory anymore

# loop over the scenarios and constrainInputs ----------------------------------

for(scenario in c("FF", "BR")){ #"Baseline", "CC", "CI", )){
  for(constrainInputs in c(TRUE, FALSE)){
    
      # load the raster prediction stack ---------------------------------------
      # scenario <- "CC"   #Baseline, CC, CI, FF, BR
      # constrainInputs <- TRUE
      resolution <- "100m"
    
      pred_stack <- rast(paste0(results_folder, version_folder, "/",
                                scenario, "/",
                                "forest_pred_stack_", resolution,"_",scenario,".tif"))
      
      # if constrained, constrain this
      if(constrainInputs){
        for(i in 1:nlyr(pred_stack)){
          pred_stack[[i]] <-  clamp(pred_stack[[i]], lower=min(Xoc[,i]), upper=max(Xoc[,i]), values=TRUE)
        }
      }
      
      # rep-wise calculations --------------------------------------------------
      # run for each rep at a time:
      
      foreach (i = 1:nRep,
               .packages = c("terra", "foreach")
      ) %do% {
      
        # for each species, pull out the intercepts, the betalpsi 1, 2, and 3
        foreach (k = 1:(rareSpp)) %do% {
          cat("spp", k, "rep", i, date(), "\n")
          out_logitpsi <- posterior_matrix[i, paste0("lpsi", "[", k, "]")]  +
            posterior_matrix[i, paste0("betalpsi", "[", k, ", 1]")] * pred_stack[[1]] + 
            posterior_matrix[i, paste0("betalpsi", "[", k, ", 2]")] * pred_stack[[2]] +
            posterior_matrix[i, paste0("betalpsi", "[", k, ", 3]")] * pred_stack[[3]] +
            posterior_matrix[i, paste0("betalpsi", "[", k, ", 4]")] * pred_stack[[4]]
          
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
          cat("spp", k, "rep", i, date(), "\n")
          out_psi <- rast(paste0(out_temp_psi_folder, 
                             "psi_sp", k-1,  ".tif"))
          writeRaster(out_psi, 
                      paste0(out_temp_psi_folder, 
                             "psi_sp", k, ".tif"),
                      overwrite=TRUE)
        }
        
        psi_stack <- rast(paste0(out_temp_psi_folder,
                                 "psi_sp", 1:M, ".tif"))
        
        # derived quantities ---------------------------------------------------
        cat("SR_site_all rep", i, date(), "\n")
        SR_site_all <- sum(psi_stack)             # Number of occurring species at each site (sum of psi)
        writeRaster(SR_site_all, 
                    paste0(out_temp_rep_folder, 
                           "SR_site_all", "_rep", i, ".tif"),
                    overwrite=TRUE)
        rm(SR_site_all)
        cat("SR_site_fspp rep", i, date(), "\n")
        SR_site_fspp <- sum(psi_stack[[which(fspp==TRUE)]])       # forest dependent species
        writeRaster(SR_site_fspp, 
                    paste0(out_temp_rep_folder, 
                           "SR_site_fspp", "_rep", i, ".tif"),
                    overwrite=TRUE)
        rm(SR_site_fspp)
        cat("SR_site_mig rep", i, date(), "\n")
        SR_site_mig <- sum(psi_stack[[which(mig==TRUE)]])          # migrants
        writeRaster(SR_site_mig, 
                    paste0(out_temp_rep_folder, 
                           "SR_site_mig", "_rep", i, ".tif"),
                    overwrite=TRUE)
        rm(SR_site_mig)
        cat("SR_site_fnDiet rep", i, date(), "\n")
        SR_site_fnDiet <- sum(psi_stack[[which(fnDiet==TRUE)]])   # fruit/nectar species
        writeRaster(SR_site_fnDiet, 
                    paste0(out_temp_rep_folder, 
                           "SR_site_fnDiet", "_rep", i, ".tif"),
                    overwrite=TRUE)
        rm(SR_site_fnDiet)
        cat("SR_site_invDiet rep", i, date(), "\n")
        SR_site_invDiet <- sum(psi_stack[[which(invDiet==TRUE)]]) # invertebrate species
        writeRaster(SR_site_invDiet, 
                    paste0(out_temp_rep_folder, 
                           "SR_site_invDiet", "_rep", i, ".tif"),
                    overwrite=TRUE)
        rm(SR_site_invDiet)
        
        file.remove(list.files(tempdir(), full.names = T, pattern = "^spat"))
        
      } # end for rep i loop 
      
      # stopCluster()
      # 252 reps takes about 2.5hrs
      
      # then combine the reps --------------------------------------------------
      
      focalMetricList <- c("SR_site_all", "SR_site_fspp", "SR_site_mig", "SR_site_fnDiet", "SR_site_invDiet")
      
      for(focalMetric in focalMetricList){
        mystack <- rast(list.files(out_temp_rep_folder, pattern = focalMetric, full.names = TRUE))
        writeRaster(mystack, paste0(results_folder, version_folder, 
                                    scenario, "/forest_", 
                                    resolution, "_",
                                    ifelse(constrainInputs, "constrained", "unconstrained"),
                                    "/", focalMetric, ".tif"),
                    overwrite=TRUE)
      }
      
  } # end for constrainInputs loop
} # end scenario loop

# post processing options:
# SR_mean <- mean(mystack)
# SR_median <- median(mystack)
# SR_HDI_low <- app(mystack, fun=function(i) hdi(i)$CI_low) # will give warnings re NA cells is OK.
# SR_HDI_high <- app(mystack, fun=function(i) hdi(i)$CI_high)
# SR_HDI_width <- SR_HDI_high - SR_HDI_low
# 
# SR_mean
# plot(SR_mean)
# density(SR_mean)


# session info -----------------------------------------------------------------
date()
sessionInfo() 
