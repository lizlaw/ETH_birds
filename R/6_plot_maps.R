# merge forest and farm rasters

# libraries --------------------------------------------------------------------
library(tidyverse)
library(terra)
library(rasterVis)

# load and specify data --------------------------------------------------------
wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds"
results_folder <- paste0(wd, "/Results/")  
version_folder <- "v10/"

# combine all the reps for all the metrics, for each scenario ------------------
for(scenario in c("Baseline", "CC", "CI", "FF", "BR")){
  for(constrainInputs in c(TRUE, FALSE)){
    resolution <- "100m"
    
    rep_summary <- function(focalMetric, sumfun=mean){
      forest_reps <- rast(paste0(results_folder, version_folder, scenario, "/forest_", 
                                 resolution, "_",
                                 ifelse(constrainInputs, "constrained", "unconstrained"),
                                 "/", focalMetric, ".tif"))
      farm_reps <- rast(paste0(results_folder, version_folder, scenario, "/farm_", 
                               resolution, "_",
                               ifelse(constrainInputs, "constrained", "unconstrained"),
                               "/", focalMetric, ".tif"))
      
      # get mean value from the reps 
      merge(sumfun(farm_reps), sumfun(forest_reps))
    }
    
    SR_site_stack <- rast(list(
      SR_site_all = rep_summary("SR_site_all", mean),
      SR_site_fspp = rep_summary("SR_site_fspp", mean),
      SR_site_mig = rep_summary("SR_site_mig", mean),
      SR_site_fnDiet = rep_summary("SR_site_fnDiet", mean),
      SR_site_invDiet = rep_summary("SR_site_invDiet", mean)
    ))
    
    writeRaster(SR_site_stack, paste0(results_folder, version_folder, 
                                      scenario, "/meanMetrics_", 
                                      resolution, "_",
                                      ifelse(constrainInputs, "constrained", "unconstrained"),
                                      ".tif"),
                overwrite=TRUE)
  } 
}

# plot, combine by scenario and metric -----------------------------------------------

plot_by_scenario <- function(focalMetric, fillName, constrainInputs, width_adj = c(1,1)){
      Baseline <- rast(paste0(results_folder, version_folder, "Baseline", 
                              "/meanMetrics_", resolution, "_",
                              ifelse(constrainInputs, "constrained", "unconstrained"),
                              ".tif"))[[focalMetric]]
      CC <- rast(paste0(results_folder, version_folder, "CC", 
                        "/meanMetrics_", resolution, "_",
                        ifelse(constrainInputs, "constrained", "unconstrained"),
                        ".tif"))[[focalMetric]]
      CI <- rast(paste0(results_folder, version_folder, "CI", 
                        "/meanMetrics_", resolution, "_",
                        ifelse(constrainInputs, "constrained", "unconstrained"),
                        ".tif"))[[focalMetric]]
      FF <- rast(paste0(results_folder, version_folder, "FF", 
                        "/meanMetrics_", resolution, "_",
                        ifelse(constrainInputs, "constrained", "unconstrained"),
                        ".tif"))[[focalMetric]]
      BR <- rast(paste0(results_folder, version_folder, "BR", 
                        "/meanMetrics_", resolution, "_",
                        ifelse(constrainInputs, "constrained", "unconstrained"),
                        ".tif"))[[focalMetric]]
      
      CC <- resample(CC, Baseline)
      CI <- resample(CI, Baseline)
      FF <- resample(FF, Baseline)
      BR <- resample(BR, Baseline)
      
      SR_site_stack <- rast(list(
        Baseline = Baseline, 
        `Coffee Invest` = CI,
        `Food First` = FF,
        `Cash Crops` = CC,
        `Biosphere Reserve` = BR
      ))
      
      mymin <- min(global(SR_site_stack, min, na.rm=TRUE))
      mymax <- max(global(SR_site_stack, max, na.rm=TRUE))
      
      p1 <- gplot(SR_site_stack[[1]]) + 
        geom_raster(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_viridis_c(limits = c(mymin, mymax), na.value = NA) +
        coord_equal() +
        theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "none")
      
      p2 <- gplot(SR_site_stack[[2:5]]) + 
        geom_raster(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_viridis_c(limits = c(mymin, mymax), na.value = NA) +
        coord_equal() +
        labs(fill=fillName) +
        theme_bw() + theme(axis.text.x = element_text(colour = "white",angle=90, vjust = 0.5), 
                           axis.title.x = element_text(colour = "white"), 
                           axis.text.y = element_blank(),
                           axis.title.y = element_blank())
      
      cowplot::plot_grid(p1, p2, rel_widths = width_adj)
}

# all -----

plot_by_scenario("SR_site_all", 
                 "Species\nrichness\n(all species)",
                 constrainInputs = FALSE,
                 width_adj = c(0.935,1))
ggsave(paste0(results_folder, version_folder,"Plot_SR_site_all_mean_unconstrained.pdf"))

plot_by_scenario("SR_site_all", 
                 "Species\nrichness\n(all species,\nconstrained)",
                 constrainInputs = TRUE,
                 width_adj = c(0.935,1))
ggsave(paste0(results_folder, version_folder,"Plot_SR_site_all_mean_constrained.pdf"))

## fspp -----

plot_by_scenario("SR_site_fspp", 
                 "Species\nrichness\n(forest\nspecies)",
                 constrainInputs = FALSE,
                 width_adj = c(0.985,1))
ggsave(paste0(results_folder, version_folder,"Plot_SR_site_fspp_mean_unconstrained.pdf"))

plot_by_scenario("SR_site_fspp", 
                 "Species\nrichness\n(forest\nspecies,\nconstrained)",
                 constrainInputs = TRUE,
                 width_adj = c(0.93,1))
ggsave(paste0(results_folder, version_folder,"Plot_SR_site_fspp_mean_constrained.pdf"))


## map of forest/farm
forest_reps <- rast(paste0(results_folder, version_folder, scenario, "/forest_", 
                           resolution, "_",
                           ifelse(constrainInputs, "constrained", "unconstrained"),
                           "/", focalMetric, ".tif"))
farm_reps <- rast(paste0(results_folder, version_folder, scenario, "/farm_", 
                         resolution, "_",
                         ifelse(constrainInputs, "constrained", "unconstrained"),
                         "/", focalMetric, ".tif"))

# get mean value from the reps 
FarmForMap <- merge(0*(mean(farm_reps)), (mean(forest_reps)>-Inf))

gplot(FarmForMap) + 
  geom_raster(aes(fill = factor(value, levels = c(0,1), labels = c("Farm", "Forest")))) +
  scale_fill_viridis_d(begin = 0.5, na.value = NA, direction = -1) +
  coord_equal() +
  labs(fill="") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5))

gplot(FarmForMap) + 
  geom_raster(aes(fill = factor(value, levels = c(0,1), labels = c("Farm", "Forest")))) +
  scale_fill_viridis_d(begin = 0.5, na.value = NA, direction = -1) +
  coord_equal() +
  labs(fill="") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

