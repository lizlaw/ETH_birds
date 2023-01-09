# climate data

library(tidyverse)
library(terra)

study_shp <- vect("/Volumes/WCC_2021/LeuphanaData/Study Area/StudyAreaWGSUTMz37.shp")
wc_temp <- rast(list.files("/Volumes/WCC_2021/LeuphanaData/WorldClim/Historical/wc2.1_30s_tavg", 
                                pattern=".tif", full.names = TRUE))
wc_precip <- rast(list.files("/Volumes/WCC_2021/LeuphanaData/WorldClim/Historical/wc2.1_30s_prec", 
                       pattern=".tif", full.names = TRUE))

study_shp <- project(study_shp, wc_temp)
wc_mean_temp <- crop(wc_temp, study_shp, mask=TRUE) %>% mean()
wc_ann_precip <- crop(wc_precip, study_shp, mask=TRUE) %>% sum()

wc_monthly_precip <- crop(wc_precip, study_shp, mask=TRUE) %>% 
  global(mean, na.rm=TRUE) %>% 
  as_tibble() %>% 
  add_column(month=1:12, .before=1) %>% 
  plot()

## how does this correlate with elevation?
elev <- rast("/Users/elaw/Desktop/LeuphanaProject/ETH_SpatialData/data_10m/Elevation/Elevation_SA_RS_10.tif")[[1]]
elev_1km <- project(elev, wc_mean_temp, method="bilinear")
plot(elev_1km)

mystack <- rast(list(elev_1km, wc_mean_temp, wc_ann_precip))
myvalues <- values(mystack, na.rm=TRUE) %>% as_tibble()
names(myvalues) <- c("elevation", "temp", "precip")

# nrow(myvalues) = 3326

ggplot(myvalues, aes(x=elevation, y=temp)) +
  geom_point() + 
  geom_smooth(method="lm")

ggplot(myvalues, aes(x=elevation, y=precip)) +
  geom_point()

## temperature - elevation line
myvalues <- myvalues %>% mutate(elev1000 = elevation/1000)
m1 <- lm(elevation ~ temp, data=myvalues)
summary(m1)

## every degree increase in temp = approx -224.8253 m in elevation

## what do the climate models predict for this region?
bc_future_list <- map(1:12,
                      ~rast(list.files("/Volumes/WCC_2021/LeuphanaData/WorldClim/Future/2021_2040/", 
                                         pattern=".tif", full.names = TRUE)[.x])
)

bc_future_temp <- map(bc_future_list, function(x) crop(x[[1]], study_shp, mask=TRUE))
bc_future_prec <- map(bc_future_list, function(x) crop(x[[12]], study_shp, mask=TRUE))

temp_stack <- rast(bc_future_temp)
prec_stack <- rast(bc_future_prec)

layernames <- list.files("/Volumes/WCC_2021/LeuphanaData/WorldClim/Future/2021_2040/", 
                                pattern=".tif", full.names = FALSE) %>% 
  str_replace(
    pattern = "wc2.1_30s_bioc_", replacement = "") %>% 
  str_replace(
    pattern = "_2021-2040.tif", replacement = ""
  )

names(temp_stack) <- layernames
names(prec_stack) <- layernames

# which models have the lowest/highest/median prediction values?

temp_vals <- temp_stack %>% global(mean, na.rm=TRUE) %>% 
  as_tibble(rownames = "GCM") %>% 
  add_column(median = temp_stack %>% global(median, na.rm=TRUE)) %>% 
  arrange(mean)

prec_vals <- prec_stack %>% global(mean, na.rm=TRUE) %>% 
  as_tibble(rownames = "GCM") %>% 
  add_column(median = prec_stack %>% global(median, na.rm=TRUE)) %>% 
  arrange(mean)

temp_vals
prec_vals

wc_mean_temp %>% global(mean, na.rm=TRUE)
wc_mean_temp %>% global(median, na.rm=TRUE)
wc_ann_precip %>% global(mean, na.rm=TRUE)
wc_ann_precip %>% global(median, na.rm=TRUE)

temp_change <-  temp_stack - wc_mean_temp
temp_change %>% global(mean, na.rm=TRUE) %>% 
  as_tibble(rownames = "GCM")
plot(temp_change) # interesting that while some have gradients (N-S, NE-SW) there is usually a corresponding reverse

prec_change <-  prec_stack - wc_ann_precip
prec_change %>% global(mean, na.rm=TRUE) %>% 
  as_tibble(rownames = "GCM")
plot(prec_change) # similarly, stronger patterns than temp here, but equally variable between GCMs. 

climate_means <- full_join(temp_vals, prec_vals, by="GCM")[,c(1,2,4)]
names(climate_means) <- c("GCM", "mean_temp", "mean_prec")
climate_means <- climate_means %>% add_row(
  GCM = "Current", 
  mean_temp = wc_mean_temp %>% global(mean, na.rm=TRUE) %>% pull(),
  mean_prec = wc_ann_precip %>% global(mean, na.rm=TRUE) %>% pull(),
)
library(ggrepel)
ggplot(climate_means, aes(x=mean_temp, y=mean_prec, color=GCM, label=GCM)) +
  geom_point() +
  geom_text_repel(size=2) +
  theme_bw()

c(0.85144, 1.49884) * 224.8

# is this change trending across altitude? -- No, doesnt appear so.
evalues <- values(elev_1km) %>% as_tibble()
cvalues <- values(temp_change) %>% as_tibble()
myvalues <- bind_cols(evalues, cvalues) %>% na.omit()

ggplot(myvalues, aes(x=Elevation_SA_RS_10, y=CanESM5_ssp245)) +
  geom_point() + 
  geom_smooth(method="lm")

### summary ------------------
#' Temperature mean is currently 17.54192. 
#' Temperature mean is forecast to range from 18.39336 (MIROC6_ssp245) 
#' to 19.04076 (CNRM-CM6-1_ssp245). 
#' This is a 0.85144 - 1.49884 change in the mean temperature.
#' On a cell-by cell basis, mean change in temperature ranges from 
#' 0.965 (MPI-ESM1-2-LR_ssp370) to 1.31 (CanESM5_ssp245)
#' 
#' This mean change in temp, based on 1 degree = 224.8m, gives a 
#' 191.4037  - 336.9392 m change in altitude. 
#' 
#' Ann precip is currently mean 1599.857, median 1616.5.
#' Ann precip is forecast to range from 1602.638 (MPI-ESM1-2-LR_ssp245) to 
#' to 1952.039 (CanESM5_ssp370). 
#' This is a 2.780535 - 352.181936 change in annual precipitation (0-22%)
#' On a cell-by cell basis, mean change in annual precip ranges from 
#' 7.14 (MPI-ESM1-2-LR_ssp370) to 352.18 (CanESM5_ssp370).
#' 
#' Patterns of climate change are quite variable, more distinct for precip than temp. 
#' 

