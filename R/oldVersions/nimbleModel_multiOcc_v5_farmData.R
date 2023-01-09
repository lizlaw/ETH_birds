#  "bayesBirdModel_multiSPecies_v5 - farmData"
#  author: "Liz Law"
#  date: "10 Jan 2022
#  workingconservation@gmail.com

# Data development ----
# Setup
library(tidyverse)

# data location and load ----
wd <- "/Users/elaw/Desktop/LeuphanaProject/BirdModelling/Leuphana_Bird_Modelling/"
BirdData <- readRDS(paste0(wd, "Data/LandscapeBirds_cleandata.RDS"))

# Note two sites only have one of the rounds
BirdData$ObsvData %>% 
  group_by(point_id) %>% 
  summarise(minRound = min(round), maxRound = max(round)) %>% 
  filter(minRound == maxRound) 
# GBD12GR01   has only round 1
# TBW05TF05   has only round 2

# remove these sites from data 
BirdData[1:3] <- BirdData[1:3] %>% 
  map(function(x){x %>% filter(!(point_id %in% c("GBD12GR01", "TBW05TF05")))})

# Bird presence/absence [site, visit, species] ----

# expand out to presence/absence, correcting the above sites to NA
BirdPA <- BirdData$ObsvData %>% 
  select(point_id, round, species_code) %>% 
  distinct() %>% 
  mutate(PA = 1) 

BirdPA <- expand_grid(point_id = unique(BirdPA$point_id),
                      round = 1:2,
                      species_code = unique(BirdPA$species_code)) %>% 
  left_join(BirdPA) %>% 
  mutate(PA = replace_na(PA, 0)) %>% 
  mutate(PA = pmap_dbl(list(point_id, round, PA),
                       ~ ifelse( ((..1 == "GBD12GR01") & (..2 == 2)) | ((..1 == "TBW05TF05") & (..2 == 1)),
                                 NA, ..3)))

# convert to binary matrix [site, visit, species]

y_all <- array(BirdPA %>% arrange(species_code, round, point_id) %>% pull(PA), 
               dim = c(length(unique(BirdPA$point_id)),
                       length(unique(BirdPA$round)),
                       length(unique(BirdPA$species_code))),
               dimnames = list(unique(BirdPA$point_id),
                               unique(BirdPA$round),
                               unique(BirdPA$species_code))
)

# Occupancy variables [site, covariate] matrix ---
## (projectable) environmental predictors [1:site, covariate] matrix, same order as sites in y

# initial matrix
occ_vars <- BirdData$SiteData %>% 
  arrange(point_id)    

occVL_joint <- c("point_id", "kebele", "habitat", "hab_details", "X", "Y", "elevation", "slope", "twi")
occVL_farm <- c("farm_type", "sidi1ha", "sidi200", "pwv1ha", "pwv500", "fl_dis", "fl_dis85")
occVL_frst <- c("forest_type", "fr_dis", "fr_dis85", "hli", "pwv2km")

# select variables and check for missing data 

OccVarsFarm <- occ_vars %>% 
  filter(habitat == "FARM") %>% 
  select(all_of(c(occVL_joint, occVL_farm)))
summary(OccVarsFarm)
is.na(OccVarsFarm) %>% sum()  ## all complete

OccVarsFrst <- occ_vars %>% 
  filter(habitat == "FOR") %>% 
  select(all_of(c(occVL_joint, occVL_frst)))
summary(OccVarsFrst)
is.na(OccVarsFrst) %>% sum()  ## all complete

## all these are complete with data, do not need a data generating process.

# Observation covariates - random variables [site, round, covariate] array ----

obs_vars <- BirdData$SurveyData %>% 
  arrange(point_id, round) %>% 
  left_join(BirdData$SiteData %>% select(point_id, habitat)) %>% 
  mutate(n_observers = str_count(observer, "\\+") + 1) %>% 
  mutate(PR_observer = str_detect(observer, "PR"))

# Select variables. We might have also liked canopy and understorey but there are too many missing variables.
ObsVarsSelect <- c("habitat",
                   "n_observers", "recording",
                   "start", "date", 
                   "visibility", "cloud")

baseObs <- expand_grid(point_id = unique(BirdPA$point_id), round = unique(BirdPA$round))

ObsVars <- baseObs %>% 
  left_join(obs_vars %>% select(point_id, round, all_of(ObsVarsSelect))) 

# check these
summary(ObsVars)

# several of these have missing variables
# I propose that, rather than estimating these from their forms (e.g. binary, categorical, etc), we might be able to make better estimates given the other information. 

# these two are the missing visits, and thus all information: GBD12GR01 r2, and TBW05TF05 r1.
# we take our cues from the sampling patterns and observations of the other points sampled alongside that point in the known rounds. Time is drawn from the known round, as this fit within the time schedule for the point sampling.
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
# # add these into the data 
# which(ObsVars$point_id == "GBD12GR01" & ObsVars$round == 2)
# ObsVars[14,]
# ObsVars[14,] <- add_GBD12GR01_r2
# which(ObsVars$point_id == "TBW05TF05" & ObsVars$round == 1)
# ObsVars[205,]
# ObsVars[205,] <- add_TBW05TF05_r1 

# now we have 3 NA in recording, and 1 NA each in visibility and cloud (ignoring understorey and canopy for now)
summary(ObsVars)

# recording
ObsVars %>% filter(is.na(recording))
mydates <- ObsVars %>% filter(is.na(recording)) %>% pull(date)
ObsVars %>% filter(date %in% mydates)
ObsVars$recording %>% na.omit %>% table()
# the first 2 seem to have no issues, so ok to assume there are recordings.
# the SGBGRC01  r2 may have had issues - there are two other issues that day. 
# however, the overwhelming trend is for no issues. Recordings likely have more id than no recording, so going with recording == 1 will underestimate (i.e. is conservative)
ObsVars$recording <- replace_na(ObsVars$recording, 1)

# visibility & cloud
ObsVars %>% filter(is.na(visibility))
mydates <- ObsVars %>% filter(is.na(visibility)) %>% pull(date)
ObsVars %>% filter(date %in% mydates)
# likely to be 2 for both
ObsVars$visibility <- replace_na(ObsVars$visibility, 2)
ObsVars$cloud <- replace_na(ObsVars$cloud, 2)

# see if we can address understorey and canopy easily
# ObsVars %>% filter(is.na(understorey)) %>% print(n=nrow(.))
# not all of these are 0 for farm 
# ObsVars %>% filter(habitat=="FARM") %>% pull(understorey) %>% table
# mySites <- ObsVars %>% filter(is.na(understorey)) %>% pull(point_id)
# ObsVars %>% filter(point_id %in% mySites) %>% print(n=nrow(.))
# this looks promising, but I just don't quite trust these are being measured in the same way or that we can transfer the scores over (due to harvest), so we leave it out for now.

# now we can transform this into the correct parameterisation array and centre scale as appropriate

# select all the required data and center-scale

farmSites <- ObsVars %>% 
  filter(habitat=="FARM") %>% 
  pull(point_id)

## FARM
## occupancy ~ elevation + slope + farm_type + fl_dis + sidi200
## observation ~ poly(start, 2) + date + visibility

occVL_farm <- c("elevation", "slope", "farm_type", "fl_dis", "sidi1ha")

OccVarsFarm <- OccVarsFarm %>% 
  filter(habitat == "FARM") %>% 
  select(point_id, all_of(occVL_farm)) %>% 
  mutate(
    elevation = elevation %>% scale() %>% .[,1],
    slope = sqrt(slope) %>% scale() %>% .[,1],
    fl_dis = sqrt(fl_dis) %>% scale() %>% .[,1],
    sidi1ha = sidi1ha %>% scale() %>% .[,1],
    farm_type = farm_type                   # binary
  ) %>% 
  select(point_id, elevation, fl_dis, sidi1ha, slope, farm_type)

obsVL_farm <- c("start", "date", "visibility")

ObsVarsFarm <- ObsVars %>% 
  filter(habitat == "FARM") %>% 
  select(point_id, round, all_of(obsVL_farm)) %>% 
  mutate(
    start = start %>% as.numeric() %>% scale() %>% .[,1],
    start1 = poly(start, 2)[,1],
    start2 = poly(start, 2)[,2],
    date = date %>% as.numeric() %>% scale() %>% .[,1],
    visibility = (visibility > 1) %>% as.numeric()      # binary
  ) %>% 
  select(point_id, round, date, start1, start2, visibility)

# add data inits params ----

# y should be binary matrix [site, visit, species], and is defined above for all sites, we need to define it here for just farm or forest
# nSite, nVisits are fixed (get dim of y)
#  M is total number of species (including, if assumed, unobserved species) 
# 159 species observed in study. 881 species in Ethiopia total, so we might need to increase this. 

# list to the data, constants - name the X variables appropriately
y <- y_all[(dimnames(y_all)[[1]] %in% farmSites), , ]       # [site, visit, species]
nSite <- dim(y)[1]
nVisits <- dim(y)[2]
Nobs <- dim(y)[3]
Nadd <- 17  # based on individual rarefaction curves for this habitat type.
M <- Nobs + Nadd
Y <- abind::abind(y, array(0, dim = c(dim(y)[1:2], Nadd)))  ## add the 'missing' species - not seen so = 0, 'absence'

# additional to deal with missing site visits
# v <- Y[,,1]; v[] <- ifelse(is.na(v[]), 0, 1)  ## if site is visited or not in that round
# Y[is.na(Y)] <- 0   # modify the NA to now be zero observed. 

z <- array(apply(Y, c(1,3), sum, na.rm = TRUE), dim=dim(Y)[c(1,3)]); z[] <- ifelse(z > 0, 1, NA) ## we know that if we observe them in a survey they are there, so can partially fill this
w <- c(rep(1, Nobs), rep(NA, Nadd)) ## similarly if we observe the species in our study they are there

Xoc <- OccVarsFarm %>% select(-point_id) %>% as.matrix() %>% array(dim = dim(.))   # [site, vars]
Xob <- array(NA, dim = c(nSite, nVisits, dim(ObsVarsFarm)[2]-2)) # [site, visit, vars]
Xob[,,1] <- ObsVarsFarm %>% arrange(round, point_id) %>% pull(date)
Xob[,,2] <- ObsVarsFarm %>% arrange(round, point_id) %>% pull(start1)
Xob[,,3] <- ObsVarsFarm %>% arrange(round, point_id) %>% pull(start2)
Xob[,,4] <- ObsVarsFarm %>% arrange(round, point_id) %>% pull(visibility)

# constain to best functioning variables
Xoc <- Xoc[,1:3]
Xob <- Xob[,,1:3]

npsi <- dim(Xoc)[2]
np <- dim(Xob)[3]

data <- list(Y = Y, 
             # v = v,
             z = z,
             w = w,
             Xoc = Xoc,
             Xob = Xob)  
constants <- list(nSite = nSite, 
                  nVisits = nVisits, 
                  Nobs = Nobs,
                  Nadd =  Nadd, 
                  M = M,
                  npsi = npsi,
                  np = np)

# Initial values - match up with the number of observed vars used in the model
omega <- Nobs/M
wst <- rep(1, M)                   # Simply set everybody at occurring
zst <- array(1, dim = c(nSite, M)) # ditto
lpsi <- rnorm(n = M)
betalpsi <- array(rnorm(n = M * npsi), dim = c(M, npsi))
lp <- rnorm(n = M)
betalp <- array(rnorm(n = M * np), dim = c(M, np))

inits <- list(
  omega = omega,
  z = zst, 
  w = wst, 
  lpsi = lpsi, 
  betalpsi = betalpsi, 
  lp = rnorm(n = M), 
  betalp = betalp #,
  #Yb = Y
)

# Set 1
params1 <- c("mu.lpsi", "sd.lpsi", "mu.betalpsi", "sd.betalpsi", 
             "mu.lp", "sd.lp", "mu.betalp", "sd.betalp", 
             "omega", "SR_site_all", "SR_site_obs", "n0", "Ntotal")

