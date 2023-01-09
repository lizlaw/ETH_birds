##  "bayesBirdModel_multiSPecies_v7 - allData"
##  author: "Liz Law"
##  date: 21 Dec 2022
##  workingconservation@gmail.com

## AIMS ----
#' extract and format data (observation, site variables, and species traits)
#' Script comments out all unnecessary code and comments, but these are retained for clarity.

## Setup ----
library(tidyverse)     # programming and data manipulation
library(bestNormalize) # scaling factors

## data location and load ----
wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/ETH_birds/"
BirdData <- readRDS(paste0(wd, "Data/LandscapeBirds_cleandata.RDS"))
version <- "v7"

## remove extra variables ----
## Remove two sites that only have one of the rounds
# BirdData$ObsvData %>% 
#   group_by(point_id) %>% 
#   summarise(minRound = min(round), maxRound = max(round)) %>% 
#   filter(minRound == maxRound) 
## GBD12GR01   has only round 1
## TBW05TF05   has only round 2
## remove these sites from data 
BirdData[1:3] <- BirdData[1:3] %>% 
  map(function(x){x %>% filter(!(point_id %in% c("GBD12GR01", "TBW05TF05")))})

## Bird presence/absence [site, visit, species] ----
## expand out to presence/absence, correcting the above sites to NA
BirdPA <- BirdData$ObsvData %>% 
  select(point_id, round, species_code) %>% 
  distinct() %>% 
  mutate(PA = 1) 

BirdPA <- expand_grid(point_id = unique(BirdPA$point_id),
                      round = 1:2,
                      species_code = unique(BirdPA$species_code)) %>% 
  left_join(BirdPA) %>% 
  mutate(PA = replace_na(PA, 0)) 

## remove all the species identified to genus only ----
BirdPA <- BirdPA %>% filter(species_code %>% str_detect(pattern = "_sp$", negate = TRUE))

## convert to binary matrix [site, visit, species] ----
y_all <- array(BirdPA %>% arrange(species_code, round, point_id) %>% pull(PA), 
               dim = c(length(unique(BirdPA$point_id)),
                       length(unique(BirdPA$round)),
                       length(unique(BirdPA$species_code))),
               dimnames = list(unique(BirdPA$point_id),
                               unique(BirdPA$round),
                               unique(BirdPA$species_code))
)

## order of species ----
## by the total number of observations per species (rather than the current alphabetical)
oderPA <- apply(y_all, 3, sum) %>% order(decreasing=TRUE)
y_all <- y_all[,,oderPA]

## Occupancy variables [site, covariate] matrix ----
## (projectable) environmental predictors [1:site, covariate] matrix, same order as sites in y
## initial matrix
occ_vars <- BirdData$SiteData %>% 
  arrange(point_id)    

occVL_joint <- c("point_id", "kebele", "habitat", "hab_details", "X", "Y", "elevation", "slope", "twi")
occVL_farm <- c("farm_type", "sidi1ha", "sidi200", "pwv1ha", "pwv500", "fl_dis", "fl_dis85")
occVL_frst <- c("forest_type", "fr_dis", "fr_dis85", "hli", "pwv2km")

## select variables and check for missing data ----
## change name for farm_type and forest_type to farm_85 and forest_85
OccVarsFarm <- occ_vars %>% 
  filter(habitat == "FARM") %>% 
  select(all_of(c(occVL_joint, occVL_farm))) %>% 
  rename(farm_85 = farm_type) 
# summary(OccVarsFarm)
# is.na(OccVarsFarm) %>% sum()  ## all complete

OccVarsFrst <- occ_vars %>% 
  filter(habitat == "FOR") %>% 
  select(all_of(c(occVL_joint, occVL_frst))) %>% 
  rename(forest_85 = forest_type)
# summary(OccVarsFrst)
# is.na(OccVarsFrst) %>% sum()  ## all complete

## all these are complete with data, do not need a data generating process.

## Observation covariates - random variables [site, round, covariate] array ----
obs_vars <- BirdData$SurveyData %>% 
  arrange(point_id, round) %>% 
  left_join(BirdData$SiteData %>% select(point_id, habitat)) %>% 
  mutate(n_observers = str_count(observer, "\\+") + 1) %>% 
  mutate(PR_observer = str_detect(observer, "PR"))

## Select variables. ----
## We might have also liked canopy and understorey but there are too many missing variables.
ObsVarsSelect <- c("habitat",
                   "n_observers", "recording",
                   "start", "date", 
                   "visibility", "cloud")

baseObs <- expand_grid(point_id = unique(BirdPA$point_id), round = unique(BirdPA$round))

ObsVars <- baseObs %>% 
  left_join(obs_vars %>% select(point_id, round, all_of(ObsVarsSelect))) 

## check these
# summary(ObsVars)

## several of these have missing variables
#' rather than estimating these from their forms (e.g. binary, categorical, etc), 
#'  we might be able to make better estimates given the other information. 
#' we take our cues from the sampling patterns and observations of the other points 
#' sampled alongside that point in the known rounds. 
#' Time is drawn from the known round, as this fit within the time schedule 
#' for the point sampling.

# # these two are the missing visits, and thus all information: GBD12GR01 r2, and TBW05TF05 r1.
# ObsVars %>% filter(point_id %>% str_detect("GBD") & round == 1 & habitat == "FARM") 
# ObsVars %>% filter(date == "2016-01-12") 
# add_GBD12GR01_r2 <- list(point_id = "GBD12GR01", round = 2,
#                       habitat = "FARM", n_observers = 1, recording = 1, start = hms::hms(hours = 9),
#                       date = as.Date("2016-01-12"), visibility = 2, cloud = 1) #, understorey = 0, canopy = 0)
# 
# nearTBW <- ObsVars %>% filter(point_id %>% str_detect("TBW") & round == 2 & date == "2016-02-04") %>% pull(point_id)
# ObsVars %>% filter(round == 1 & point_id %in% nearTBW) 
# add_TBW05TF05_r1 <- list(point_id = "TBW05TF05", round = 1,
#                         habitat = "FARM", n_observers = 2, recording = 1, start = hms::hms(hours = 8, minutes = 46),
#                       date = as.Date("2015-12-31"), visibility = 3, cloud = 0) #, understorey = 0, canopy = 0)
# 
## add these into the data 
# which(ObsVars$point_id == "GBD12GR01" & ObsVars$round == 2)
# ObsVars[14,]
# ObsVars[14,] <- add_GBD12GR01_r2
# which(ObsVars$point_id == "TBW05TF05" & ObsVars$round == 1)
# ObsVars[205,]
# ObsVars[205,] <- add_TBW05TF05_r1 

## now we have 3 NA in recording, and 1 NA each in visibility and cloud 
## (ignoring understorey and canopy for now)
# summary(ObsVars)

## recording ----
# ObsVars %>% filter(is.na(recording))
# mydates <- ObsVars %>% filter(is.na(recording)) %>% pull(date)
# ObsVars %>% filter(date %in% mydates)
# ObsVars$recording %>% na.omit %>% table()
## the first 2 seem to have no issues, so ok to assume there are recordings.
## the SGBGRC01  r2 may have had issues - there are two other issues that day. 
## however, the overwhelming trend is for no issues. Recordings likely have more id than no recording, so going with recording == 1 will underestimate (i.e. is conservative)
ObsVars$recording <- replace_na(ObsVars$recording, 1)

## visibility & cloud ----
# ObsVars %>% filter(is.na(visibility))
# mydates <- ObsVars %>% filter(is.na(visibility)) %>% pull(date)
# ObsVars %>% filter(date %in% mydates)
## likely to be 2 for both
ObsVars$visibility <- replace_na(ObsVars$visibility, 2)
ObsVars$cloud <- replace_na(ObsVars$cloud, 2)

## see if we can address understorey and canopy easily ----
# ObsVars %>% filter(is.na(understorey)) %>% print(n=nrow(.))
## not all of these are 0 for farm 
# ObsVars %>% filter(habitat=="FARM") %>% pull(understorey) %>% table
# mySites <- ObsVars %>% filter(is.na(understorey)) %>% pull(point_id)
# ObsVars %>% filter(point_id %in% mySites) %>% print(n=nrow(.))
## this looks promising, but I just don't quite trust these are being measured in the same way or that we can transfer the scores over (due to harvest), so we leave it out for now.

#' now we can transform this into the correct parameterisation array 
#' and centre scale as appropriate

## select all the required data and center-scale ----
farmSites <- ObsVars %>% 
  filter(habitat=="FARM") %>% 
  pull(point_id)

forestSites <- ObsVars %>% 
  filter(habitat=="FOR") %>% 
  pull(point_id)

## FARM ----
## occupancy ~ elevation + farm_85 + fl_dis + sidi1ha + slope
## observation ~ date + poly(start, 2) + visibility

occVL_farm <- c("elevation", "farm_85", "fl_dis", "sidi1ha", "slope")

OccVarsFarm <- OccVarsFarm %>% 
  filter(habitat == "FARM") %>% 
  select(point_id, all_of(occVL_farm)) 

## use bestNormalize (and poly)for transformations
elevMod_farm <- poly(OccVarsFarm$elevation, 2)
OccVarsFarm_scaleModels <- list(
  elevation = elevMod_farm,
  elevation1 = bestNormalize::center_scale(elevMod_farm[,1]),
  elevation2 = bestNormalize::center_scale(elevMod_farm[,2]),
  slope = bestNormalize::sqrt_x(OccVarsFarm$slope, standardize = TRUE),
  fl_dis = bestNormalize::sqrt_x(OccVarsFarm$fl_dis, standardize = TRUE),
  sidi1ha = bestNormalize::center_scale(OccVarsFarm$sidi1ha),
  farm_85 = bestNormalize::no_transform(OccVarsFarm$farm_85) # binary
)
saveRDS(OccVarsFarm_scaleModels, 
        paste0(wd, "Data/nimbleModel_multiOcc_", version, "_OccVarsFarm_scaleModels.RDS"))
  
OccVarsFarm <- OccVarsFarm %>% 
  mutate(
    elevation1 = predict(OccVarsFarm_scaleModels$elevation1),
    elevation2 = predict(OccVarsFarm_scaleModels$elevation2),
    farm_85 = predict(OccVarsFarm_scaleModels$farm_85),  # binary,
    fl_dis = predict(OccVarsFarm_scaleModels$fl_dis),
    sidi1ha = predict(OccVarsFarm_scaleModels$sidi1ha),
    slope = predict(OccVarsFarm_scaleModels$slope)
  ) %>% 
  select(point_id, elevation1, elevation2, farm_85, fl_dis, sidi1ha, slope)

obsVL_farm <- c("start", "date", "visibility")

ObsVarsFarm <- ObsVars %>% 
  filter(habitat == "FARM") %>% 
  select(point_id, round, all_of(obsVL_farm)) %>% 
  mutate(
    start = start %>% as.numeric(),
    start1 = poly(start, 2)[,1] %>% scale() %>% .[,1],
    start2 = poly(start, 2)[,2] %>% scale() %>% .[,1],
    date = date %>% as.numeric() %>% scale() %>% .[,1],
    visibility = (visibility > 1) %>% as.numeric()      # binary
  ) %>% 
  select(point_id, round, date, start1, start2, visibility)

## FOREST ----
## occupancy ~ elevation + forest_85 + fr_dis + hli + slope
## observation ~ date + poly(start,2) + n_observers

occVL_forest <- c("elevation", "forest_85", "fr_dis", "hli", "slope")

OccVarsFrst <- OccVarsFrst %>% 
  filter(habitat == "FOR") %>% 
  select(point_id, all_of(occVL_forest)) 

elevMod_forest <- poly(OccVarsFrst$elevation, 2)
OccVarsFrst_scaleModels <- list(
  elevation = elevMod_forest,
  elevation1 = bestNormalize::center_scale(elevMod_forest[,1]),
  elevation2 = bestNormalize::center_scale(elevMod_forest[,2]),
  forest_85 = bestNormalize::no_transform(OccVarsFrst$forest_85), # binary
  fr_dis = bestNormalize::sqrt_x(OccVarsFrst$fr_dis, standardize = TRUE),
  hli = bestNormalize::center_scale(OccVarsFrst$hli, standardize = TRUE),
  slope = bestNormalize::sqrt_x(OccVarsFrst$slope, standardize = TRUE)
  )
saveRDS(OccVarsFrst_scaleModels, 
        paste0(wd, "Data/nimbleModel_multiOcc_", version, "_OccVarsFrst_scaleModels.RDS"))

OccVarsFrst <- OccVarsFrst %>%
  mutate(
    elevation1 = predict(OccVarsFrst_scaleModels$elevation1),
    elevation2 = predict(OccVarsFrst_scaleModels$elevation2),
    slope = predict(OccVarsFrst_scaleModels$slope),
    hli = predict(OccVarsFrst_scaleModels$hli),
    fr_dis = predict(OccVarsFrst_scaleModels$fr_dis),
    forest_85 = predict(OccVarsFrst_scaleModels$forest_85) # binary
  ) %>% 
  select(point_id, elevation1, elevation2, forest_85, fr_dis, hli, slope)

obsVL_forest <- c("start", "date", "cloud")

ObsVarsFrst <- ObsVars %>% 
  filter(habitat == "FOR") %>% 
  select(point_id, round, all_of(obsVL_forest)) %>% 
  mutate(
    start = start %>% as.numeric(),
    start1 = poly(start, 2)[,1] %>% scale() %>% .[,1],
    start2 = poly(start, 2)[,2] %>% scale() %>% .[,1],
    date = date %>% as.numeric() %>% scale() %>% .[,1],
    cloud = (cloud > 1) %>% as.numeric()      # binary
  ) %>% 
  select(point_id, round, date, start1, start2, cloud)

## species trait groups of interest ----
sppTraits <- left_join(
  BirdPA %>% select(species_code) %>% unique(),
  BirdData$TraitData %>% select(species_code=spcode, 
                                diet_5cat, 
                                for_strat5levels, 
                                size=body_mass_value_elthon,
                                mig=migration_birdlife, # mig2=migration_hbw,
                                dist=endemic_based_on_distribution_from_bird_life_datazone,
                                #dist2=biome_restricted_afrotropical,
                                forDep=forest_dependency_birdlife,
                                status)
)

# table(sppTraits$diet_5cat)
# table(sppTraits$for_strat5levels)
# hist(sppTraits$size %>% as.numeric() %>% log10())
# table(sppTraits$mig)
# table(sppTraits$dist)
# table(sppTraits$forDep)
# table(sppTraits$status)

## common/rare species (can be used to simplify parameters further) ----
sppTraits <- sppTraits %>% 
  add_column(
    sumObserved = apply(y_all,3,sum)
  ) %>% 
  mutate(
    sumObs1 = sumObserved<=1,
    sumObs2 = sumObserved<=2
  )

## forest specialists ----
sppTraits <- sppTraits %>% 
  mutate(
    forSpec = forDep %in% c("high", "medium") #, "low")   # 54 spp
  )

## migrants ----
sppTraits <- sppTraits %>% 
  mutate(
    migYes = mig %in% c("full migrant", "altitudinal migrant") # 30 spp
  )

## diet specialists - fruit&nect ----
sppTraits <- sppTraits %>% 
  mutate(
    FruiNect = diet_5cat %in% c("FruiNect"),  # 18 spp
    Invertebrate = diet_5cat %in% c("Invertebrate") # 73 spp
  )

# #list to the data, constants - name the X variables appropriately
saveRDS(list(version=paste0("allData_", version), 
             y_all=y_all,
             forestSites=forestSites,
             farmSites=farmSites,
             OccVarsFrst=OccVarsFrst,
             ObsVarsFrst=ObsVarsFrst,
             OccVarsFarm=OccVarsFarm,
             ObsVarsFarm=ObsVarsFarm,
             sppTraits=sppTraits),
        paste0(wd, "Data/nimbleModel_multiOcc_", version, "_allData.RDS")
)

