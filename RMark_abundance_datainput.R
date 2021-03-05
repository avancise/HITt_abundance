library(sf)
library(sp)
library(magrittr)
library(tidyverse)
library(rgdal)
library(lubridate)
library(reshape2)
library(zoo)

####Generate data input files for POPAN abundance estimation####
####            Written by Amy M. Van Cise                  ####
################################################################

#set file location and name of data file to be used
wd <- "C:/Users/Amy/Google Drive/05 CRC/01 Tursiops abundance"
photoIDfile <- "TursiopsPhotoIDJuly2020v6.csv"

setwd(wd)

#levels used to order tables and figures consistently throughout the project
subarealevels <- c("KA","KB","KC","KD","OA","OB","MA","MB","HA","HB")
arealevels <- c("Kauai/Niihau", "Oahu", "Maui Nui", "Hawaii")

#read GPS effort data
# transect.data <- read.csv("Effort data with on off effort_Nov2019.csv", as.is = TRUE) %>% 
#   mutate(date=as.POSIXct(strptime(Date_, format="%m/%d/%Y"))) %>% 
#   mutate(Island = replace(Island, Island %in% c("Lanai", "Maui"), "Maui.Nui")) %>% 
#   mutate(Island = replace(Island, Island == "Kauai", "Kauai.Niihau")) %>%
#   mutate(year = year(date)) %>% 
#   st_as_sf(coords = c("Long_dd","Lat_dd"), crs = (4135))

#create capture dataframe
capdat <-
  read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  mutate(date=as.POSIXct(strptime(Date..dd.MMM.yyyy., format="%d-%b-%Y"))) %>% 
  mutate(Depth=as.numeric(Depth)) %>% 
  mutate(year = year(date)) %>% 
  filter(Area != "Unknown" & Area != "Offshore" & Area != "NWHI") %>% 
  filter(as.numeric(Depth) < 3000 | is.na(Depth)) %>% 
  filter(Distinctiveness..4..very..3...average..2...slightly..1...not. >= 3 & 
           Best.photo.quality..4...excellent..3.good..2.fair..1.poor. >= 3) %>% 
  filter(Partially.entered.matched.encounter...Complete.Incomplete. == "complete") %>% 
  filter(Island != "Kaula Rock") %>% 
  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>%
  mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>%
  filter(!(Area == "Hawaii" & abs(Long) < 155.82)) %>%
  filter(year > 1999 & year <= 2018) %>%
  mutate(Area = replace(Area, Area %in% c("4-island","4-island ","4-Island","Maui"), "Maui Nui")) %>%
  mutate(detect = 1) %>% 
  dplyr::rename(Sight.type = Original..within.year..between.year..final.Note.between.island.only.used.for.between.area.matches) %>% 
  mutate(Sight.type = replace(Sight.type, Sight.type %in% c("Original", "Original ", "original"), "Original")) %>% 
  mutate(Sight.type = replace(Sight.type, Sight.type %in% c("within-day", "Within-day "), "Within-day")) %>% 
  mutate(Sight.type = replace(Sight.type, Sight.type == "within-year", "Within-year"))

#count number of individuals seen only once
onex.dat <- capdat %>% group_by(ID..) %>% filter(n() == 1) %>% ungroup() %>% tally() %>% .$n

#filter data for subarea assignment
locdat <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  mutate(date=as.POSIXct(strptime(Date..dd.MMM.yyyy., format="%d-%b-%Y"))) %>% 
  mutate(Depth=as.numeric(Depth)) %>% 
  mutate(year = year(date)) %>% 
  mutate(Area = replace(Area, Area %in% c("4-island","4-island ","4-Island","Maui"), "Maui Nui")) %>%
  filter(Area != "Unknown" & Area != "Offshore" & Area != "NWHI") %>% 
  filter(as.numeric(Depth) < 3000 | is.na(Depth)) %>%
  filter(Island != "Kaula Rock") %>% 
  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
  filter(!is.na(Long)) %>%
  mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>% 
  filter(!(Area == "Hawaii" & abs(Long) < 155.82)) %>%
  group_by(ID..) %>% 
  arrange(year) %>% 
  slice(1) %>% 
  ungroup()

#get coastline data
setwd("C:/Users/Amy/Google Drive/05 CRC/01 Tursiops abundance/Raw data/HI coastlines shapefile")
coast <- readOGR(".", layer = "Coastline")
coast <- st_as_sf(coast)# make sf object
coastr <- st_transform(coast, crs = 32604)# transform to projection of data
#plot(coastr["agency"], main=NULL)

setwd("C:/Users/Amy/Google Drive/05 CRC/01 Tursiops abundance/Raw data/Hawaii_DepthContours")
depth <- readOGR(".", layer = "Ocean_Depth")
depth <- st_as_sf(depth)# make sf object
depthr <- st_transform(depth, crs = 32604)# transform to projection of data
depthr <- depthr %>% filter(depth == 1000)
#plot(depthr["agency"], main=NULL)

#mapping visuals to manually generate subarea designations
# subdat <- capdat %>% filter(Area == "Hawaii") %>%  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>%
#   filter(!is.na(Long)) %>%
#   mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>%
#   select(ID..,Lat,Long, year)
# #ggplot to visualize thresholds
# ggplot(data = coastr)+
#   geom_sf(fill = "#00695C", alpha=0.5) +
#   coord_sf(crs = st_crs(4135)) +
#   xlim(-156.3,-154.6) +
#   ylim(18.9,20.5) +
#   geom_point(data=subdat, aes(x=Long, y=Lat, color=year)) +
#   geom_hline(yintercept = 19.71) +
#   geom_vline(xintercept = -155.82) +
#   theme_classic(base_size = 18) +
#   theme(legend.position = "none")

#assign indivudals to subarea by first sighting
subdat <- locdat %>% 
  mutate(subarea = {ifelse(.$Area == "Kauai/Niihau",
                           ifelse(abs(.$Long) < 159.5, "KD",
                                  ifelse(.$Lat < 22 & abs(.$Long) < 159.94, "KC", 
                                         ifelse(.$Lat > 22 & abs(.$Long) < 159.94,"KB","KA"))),
                          ifelse(.$Area == "Hawaii",
                                 ifelse(.$Lat <= 19.71, "HB", "HA"),
                          ifelse(.$Area == "Maui Nui",
                                 ifelse(.$Lat >= 21 & abs(.$Long) > 156.8 | abs(.$Long) >= 157.5, "MA","MB"),
                          ifelse(.$Area == "Oahu",
                                 ifelse(.$Lat <= 21.5, "OB","OA"),"NA"))))}) %>% 
  dplyr::select(ID..,subarea) 

#create mapping dataframe
mapdat <- locdat %>% left_join(subdat, by = "ID..")

#add subarea assignment to capdat
capdat <- capdat %>% 
  left_join(subdat, by = "ID..") 

setwd(wd)