library(magrittr)
library(dplyr)

####Calcuate resights within the same subarea####
####     Written by Amy M. Van Cise          ####
#################################################

# capdat dataset comes from RMark_abundance_datainput.R

#assign subarea for each sighting in the dataset
sightlocdat <- capdat %>% 
  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
  filter(!is.na(Long)) %>%
  mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>%
  mutate(sight.subarea = {ifelse(.$Area == "Kauai/Niihau",
                                 ifelse(abs(.$Long) < 159.5, "KA",
                                        ifelse(.$Lat < 22 & abs(.$Long) < 159.94, "KB", 
                                               ifelse(.$Lat > 22 & abs(.$Long) < 159.94,"KC","KD"))),
                           ifelse(.$Area == "Hawaii",
                                  ifelse(.$Lat <= 19.73, "HB", "HA"),
                                  ifelse(.$Area == "Maui Nui",
                                         ifelse(.$Lat >= 21 & abs(.$Long) > 156.8 | abs(.$Long) >= 157.5, "MA","MB"),
                                         ifelse(.$Area == "Oahu",
                                                ifelse(.$Lat <= 21.5, "OB","OA"),"NA"))))})

#does encounter subarea match first encounter subarea?
sub.sight.prob <- sightlocdat %>% 
  group_by(ID..) %>% 
  mutate(insubarea = ifelse(sight.subarea == subarea, TRUE, FALSE)) %>% 
  slice(2:n()) %>% 
  summarize(prob.subarea=sum(insubarea==TRUE)/n()) %>% 
  mutate(prob.subarea = round(prob.subarea,2))

#mean probability of being found in the same subarea as the first encounter
meanprob.sub <- mean(sub.sight.prob$prob.subarea)

#dataset of encounters where individuals were found in a different subarea
sub.sight.diff <- sightlocdat %>% 
  group_by(ID..) %>% 
  mutate(insubarea = ifelse(sight.subarea == subarea, TRUE, FALSE)) %>% 
  slice(2:n()) %>% 
  filter(insubarea==FALSE)