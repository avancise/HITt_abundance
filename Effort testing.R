library(RMark)
library(lubridate)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(maps)
library(sp)

####Generate summary metrics of field effort####
####     Written by Amy M. Van Cise         ####
################################################

#count individuals seen only once
uniqueinds <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  group_by(ID..) %>% 
  filter(n() == 1)

#count total individuals in dataset
numind <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  filter(!duplicated(ID..)) %>% 
  mutate(date=as.POSIXct(strptime(Date..dd.MMM.yyyy., format="%d-%b-%Y"))) %>% 
  mutate(Depth=as.numeric(Depth)) %>% 
  mutate(year = year(date)) %>%
  filter(year > 1999 & year <= 2018)

#count total number of observations
numobs <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  mutate(date=as.POSIXct(strptime(Date..dd.MMM.yyyy., format="%d-%b-%Y"))) %>% 
  mutate(Depth=as.numeric(Depth)) %>% 
  mutate(year = year(date)) %>%
  filter(year > 1999 & year <= 2018)

#count number of "pelagic" individuals removed from dataset
numpel <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  filter(Area == "Offshore" | as.numeric(Depth) > 3000 | Island == "Kaula Rock") %>% 
  group_by(ID..) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  tally()

#find ind. encountered at Kaula Rock
kaula <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  filter(Island == "Kaula Rock") %>% 
  .$ID..

kaula.inds <- read.csv(photoIDfile, stringsAsFactors = FALSE) %>% 
  filter(ID.. %in% kaula)

#calculate median ind. sightings
#capdat dataset comes from RMark_abundance_datainput.R
med.sight <- capdat %>% 
  group_by(ID..) %>% 
  summarize(numsight = n()) %>% 
  ungroup() %>% 
  summarize(median(numsight))

#calculate max ind. sightings
max.sight <- capdat %>% 
  group_by(ID..) %>% 
  summarize(numsight = n()) %>% 
  ungroup() %>% 
  summarize(max(numsight))

#calculate median individual sightings per year
med.sight.yr <- capdat %>% 
  group_by(Area,year,ID..) %>% 
  summarize(numsight = n()) %>% 
  ungroup() %>% 
  group_by(Area,year) %>% 
  summarize(median(numsight))

#summarize sightings by type (original, within year, between year) within each stock
capdatsummary <- capdat %>% 
  group_by(Area,year) %>% 
  dplyr::count(Sight.type) %>% 
  spread(Sight.type, n) %>% 
  dplyr::select(Original, "Within-year", "Between-year") %>% 
  ungroup() %>%
  replace_na(list(Original = 0, "Within-year" = 0, "Between-year" = 0)) %>% 
  mutate(Total.Sight = rowSums(.[3:5])) %>%
  group_by(Area) %>%
  mutate(Cum.Original = cumsum(Original)) %>%
  mutate(Cum.Sight = cumsum(Total.Sight)) %>% 
  ungroup()

#summarize sightings by type (original, within year, between year) within each subarea
capdatsubsum <- capdat %>% 
  filter(!is.na(subarea)) %>% 
  group_by(subarea,year) %>% 
  dplyr::count(Sight.type) %>% 
  spread(Sight.type, n) %>% 
  dplyr::select(Original, "Within-year", "Between-year") %>% 
  ungroup() %>%
  replace_na(list(Original = 0, "Within-year" = 0, "Between-year" = 0)) %>% 
  mutate(Total.Sight = rowSums(.[3:5])) %>%
  group_by(subarea) %>%
  mutate(Cum.Original = cumsum(Original)) %>%
  mutate(Cum.Sight = cumsum(Total.Sight)) %>% 
  ungroup()

#calculate number of individuals without subarea designations
nogps <- capdat[is.na(capdat$subarea),]
num.nogps <- length(unique(nogps$ID..))
numsight.nogps <- table(nogps$ID..)
num1sight.nogps <- length(which(numsight.nogps == 1)) #num ind without subarea and seen 1x

#calculate number of ind. in each subarea
sub.ind <- capdat %>% 
  filter(!is.na(subarea)) %>% 
  group_by(Area,subarea) %>% 
  summarize(Num.individuals = n_distinct(ID..)) %>% 
  ungroup() %>% 
  rename(Subarea = subarea) %>% 
  rename(Stock = Area)

#calculate number of encounters in each subarea by source of encounter
sub.enc <- capdat %>% 
  filter(!is.na(subarea)) %>% 
  distinct(Group, .keep_all = TRUE) %>% 
  group_by(Area,subarea) %>% 
  summarize(Span.of.Years = paste(min(year),"-",max(year), sep=""),
            CRC.encounters = sum(Source == "RWB"), 
            PWF.encounters = sum(Source == "PWF Line Transect" | Source == "PWF Platform of Opportunity" |
                                Source == "PWF Donation" | Source == "PWF" | Source == "PWF "),
            Other.encounters = sum(!(Source %in% c("RWB", "PWF Line Transect", "PWF Platform of Opportunity", "PWF Donation", "PWF", "PWF ")))) %>% 
  ungroup() %>% 
  mutate(Total.encounters = rowSums(.[4:6])) %>% 
  mutate(Area = factor(.$Area, levels = arealevels)) %>% 
  arrange(Area) %>% 
  rename(Subarea = subarea) %>% 
  rename(Stock = Area) %>% 
  left_join(sub.ind, by = c("Stock","Subarea")) %>% #merge with sub.ind 
  rbind(., data.frame(Stock= "Total", Subarea = "", Span.of.Years = "",
                      Num.individuals = sum(.$Num.individuals),
                      CRC.encounters = sum(.$CRC.encounters, na.rm=T), 
                      PWF.encounters = sum(.$PWF.encounters, na.rm=T), 
                      Other.encounters = sum(.$Other.encounters, na.rm=T),
                      Total.encounters = sum(.$Total.encounters, na.rm=T, check.names=FALSE))) %>% 
  rename("Num individuals" = Num.individuals, "Span of Years" = Span.of.Years, "CRC encounters" = CRC.encounters,
         "PWF encounters" = PWF.encounters, "Other encounters" = Other.encounters, "Total encounters" = Total.encounters) %>% 
  select(Stock, Subarea, "Span of Years", "Num individuals", everything())

#calculate number of encounters in each subarea and year by source of encounter
sub.enc.year <- capdat %>% 
  filter(!is.na(subarea)) %>% 
  distinct(Group, .keep_all = TRUE) %>% 
  group_by(Area,subarea,year) %>% 
  summarize(CRC.encounters = sum(Source == "RWB"), 
            PWF.encounters = sum(Source == "PWF Line Transect" | Source == "PWF Platform of Opportunity" |
                                   Source == "PWF Donation" | Source == "PWF" | Source == "PWF "),
            Other.encounters = sum(!(Source %in% c("RWB", "PWF Line Transect", "PWF Platform of Opportunity", "PWF Donation", "PWF", "PWF ")))) %>% 
  ungroup() %>% 
  #mutate(Total.encounters = rowSums(.[4:6])) %>% 
  mutate(Area = factor(.$Area, levels = arealevels)) %>% 
  arrange(Area) %>% 
  rename(Subarea = subarea) %>% 
  rename(Stock = Area) %>% 
  left_join(sub.ind, by = c("Stock","Subarea")) %>% #merge with sub.ind 
  # rbind(., data.frame(Stock= "Total", Subarea = "", year = "",
  #                     Num.individuals = sum(.$Num.individuals),
  #                     CRC.encounters = sum(.$CRC.encounters, na.rm=T), 
  #                     PWF.encounters = sum(.$PWF.encounters, na.rm=T), 
  #                     Other.encounters = sum(.$Other.encounters, na.rm=T),
  #                     Total.encounters = sum(.$Total.encounters, na.rm=T, check.names=FALSE))) %>% 
  rename("Num individuals" = Num.individuals, "CRC encounters" = CRC.encounters,
         "PWF encounters" = PWF.encounters, "Other encounters" = Other.encounters) %>% 
  select(Stock, Subarea, year, "Num individuals", everything()) %>% 
  pivot_longer(c("CRC encounters","PWF encounters","Other encounters"), names_to = "enc_type", values_to = "num_enc")

my.pal = pal_futurama()(12)[c(2:4)]
pdf(file = "Supplemental Figure S1.sightings-by-year.pdf")
ggplot(data = sub.enc.year, aes(x=year, y = num_enc, fill = enc_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Stock)+
  theme_classic(base_size = 14) +
  scale_fill_manual(values = my.pal) +
  ylab("Number of encounters") +
  labs(fill = "Encounter type")
dev.off()
  

#count number of inds seen in both Maui Nui and Oahu
MnOinds <- capdat %>% 
  filter(Area == "Maui Nui" & subarea %in% c("OB","OA") | Area == "Oahu" & subarea %in% c("MB","MA")) %>% 
  .$ID.. %>% 
  n_distinct()
  
#PWF encounters
PWFenc <- capdat %>% 
  filter(Source == "PWF Line Transect" | Source == "PWF Platform of Opportunity" |
           Source == "PWF Donation" | Source == "PWF" | Source == "PWF ")

#find kaula rock encounters
kaularock <- capdat %>% 
  filter(ID.. %in% capdat[which(Island == "Kaula Rock"),1])

#summarize by number of contributions per source
enc.cont <- capdat %>% 
  group_by(Source) %>% 
  summarize(Contributions = n()) %>% 
  arrange(desc(Contributions))

#summarize number of sightings in oahu by year
oahueffort <- capdat %>% 
  filter(Area == "Oahu") %>% 
  group_by(year) %>% 
  summarize(Sightings = n(), Individuals = length(unique(ID..)))

#number of sightings on Hawaii windward side by year
hawaiiwindeffort <- capdat %>% 
  filter(Area == "Hawaii") %>% 
  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
  filter(!is.na(Long)) %>%
  mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>%
  filter(abs(Long) < 155.8) %>% 
  group_by(year) %>% 
  summarize(Sightings = n(), Individuals = length(unique(ID..)))

hawaiiwind <- capdat %>% 
  filter(Area == "Hawaii") %>% 
  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
  filter(!is.na(Long)) %>%
  mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>%
  filter(abs(Long) < 155.8) 

hawaiiwindInd <- unique(hawaiiwind$ID..)

hawaiiwindIndresight <- capdat %>% 
  filter(ID.. %in% hawaiiwindInd) %>%
  # mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
  # filter(!is.na(Long)) %>%
  # mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>%
  group_by(ID..) %>% 
  filter(n() >1) %>% 
  n_groups()

# ggplot(data = coastr)+
#     geom_sf(fill = "#00695C", alpha=0.5) +
#     coord_sf(crs = st_crs(4135)) +
#     xlim(-156.3,-154.6) +
#     ylim(18.9,20.5) +
#     geom_point(data=hawaiiwindIndresight, aes(x=Long, y=Lat, color=ID..), position = "jitter")

#map geographic coverage of sightings by year
# oahucoverage <- capdat %>% 
#   filter(Area == "Oahu") %>% 
#   mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
#   filter(!is.na(Long)) %>%
#   mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)})
# 
# oahuplots <- list()
# oahuyears <- unique(oahucoverage$year)
# 
# for (i in 1:length(oahuyears)){
#   oahuyear <- oahucoverage %>% filter(year == oahuyears[i])
#   oahuplots[[i]] <- ggplot() +
#     geom_sf(data = coastr, fill = "#00695C", alpha = 0.5) +
#     geom_sf(data = depthr, fill = "#FFFFFF") +
#     coord_sf(crs = sf::st_crs(4135)) +
#     geom_point(data=oahuyear, aes(x=Long, y=Lat)) +
#     theme_classic(base_size = 18) +
#     #coord_map(xlim = c(-159, -157),ylim = c(20, 52))
#     xlim(-158.6,-157.5)+
#     ylim(21,22) +
#     xlab("") +
#     ylab("") +
#     annotate(geom="text", x=-157.6, y=21.97, label=oahuyears[i],size=8) +
#     theme(axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none")+
#     theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))
# }

#do.call("grid.arrange", c(oahuplots, ncol=2))