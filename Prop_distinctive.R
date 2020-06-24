library(dplyr)
library(magrittr)
library(FSA)

####Calculate proportion distinctive animals by stock and subarea####
####            Written by Amy M. Van Cise                       ####
#####################################################################

#get data
propdat <- read.csv(photoIDfile, stringsAsFactors = FALSE)
propdat$date <- as.POSIXct(strptime(propdat$Date..dd.MMM.yyyy., format="%d-%b-%Y"))

#data cleanup
propdat <- propdat %>% 
  mutate(year = year(date)) %>%
  mutate(Area = ifelse(Area == "4-island" | Area == "4-island " | Area == "Maui" | 
                        Area == "4-Island", "Maui Nui", Area)) %>% 
  filter(Area != "NWHI" & Area != "Offshore" & Area != "Unknown") %>% 
  filter(Source == "RWB")

#find encounters with > 4 distinctive, high quality IDs
encs <- propdat %>% 
  subset(Distinctiveness..4..very..3...average..2...slightly..1...not. >=3 &
           Best.photo.quality..4...excellent..3.good..2.fair..1.poor. >= 3) %>% 
  group_by(Group) %>% 
  filter(n() > 4) %>% 
  "$" (Group)

#subset propdat by encounters with >4 distinctive, high quality IDs
propdat <- propdat %>% 
  filter(Group %in% encs) 
rm(encs)

##Proportion of highly distinctive individuals
#count number of highly distinctive animals by encounter
propdist <- propdat %>% 
  filter(Best.photo.quality..4...excellent..3.good..2.fair..1.poor. >= 3) %>%
  group_by(Area,year,Group) %>% 
  summarize(prop34 = sum(Distinctiveness..4..very..3...average..2...slightly..1...not. >=3)/n())

#mean and variance in number of highly distinctive individuals by encounter
prop<- propdist %>% 
  summarize(mean = mean(prop34), sd = sd(prop34), se = se(prop34), var = var(prop34)) %>% 
  filter(year > 1999 & year < 2019)

#ANOVA test of significant differences in the proportion of high quality photos by year and stock
propdisttest <- aov(prop34~Area + year+year:Area,data=propdist)
disttestsum <- summary(propdisttest)

##Proportion of high quality photos
#count number of high quality photos by encounter
goodphotoprop <- propdat %>%
  group_by(Area,year,Group) %>% 
  summarise(propgood = sum(Best.photo.quality..4...excellent..3.good..2.fair..1.poor. >= 3)/n()) 

#mean and variance in number of good photos by encounter
goodphotosummary <- goodphotoprop %>% 
  summarize(mean = mean(propgood), sd = sd(propgood), se = se(propgood), var = var(propgood))

#ANOVA test of significant differences in the proportion of high quality photos by year and stock
photoqualitytest <- aov(propgood~year + Area, data = goodphotoprop)
photoqualsum <- summary(photoqualitytest)
