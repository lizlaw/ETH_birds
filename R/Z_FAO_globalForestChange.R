# global forest change 
# https://fra-data.fao.org/WO/fra2020/forestAreaChange/


library(tidyverse)

# forest change in 1000 ha/year
fac <- read_csv("Data/Data_raw/fra2020-forestAreaChange.csv", skip = 1) 

names(fac) <- c("country", 
                paste0("Expansion_", c("1990_2000", "2000_2010", "2010_2015", "2015_2020")),
                paste0("Deforestation_", c("1990_2000", "2000_2010", "2010_2015", "2015_2020"))
)

fac <- fac %>% mutate(
  Expansion_2010_2020 = map2_dbl(Expansion_2010_2015, Expansion_2015_2020, mean),
  Deforestation_2010_2020 = map2_dbl(Deforestation_2010_2015, Deforestation_2015_2020, mean)
  )

fac <- fac %>% mutate(
  across(starts_with("Exp"), ~-.x)
)
                       
fac <- fac %>% mutate(
  Change_1990_2000 = Deforestation_1990_2000 + Expansion_1990_2000,
  Change_2000_2010 = Deforestation_2000_2010 + Expansion_2000_2010,
  Change_2010_2020 = Deforestation_2010_2020 + Expansion_2010_2020
)

FAC <- fac %>% select(country, starts_with("Change"))

# Ethiopia
fac %>% filter(country=="Ethiopia") %>% select(starts_with("Change"))
