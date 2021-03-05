library(fasterize)
library(raster)
library(sf)
library(magrittr)
library(tidyverse)
library(geosphere)
library(hrbrthemes)
library(gridExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(rgeos)
library(ggpattern)

####    Test whetherindivdiual inter-annual movement is      ####
####  significantly less than stock inter-annual movement    ####
####             adapted by Amy M. Van Cise                  ####
#### from code by K. Alexandra Curtis (alex.curtis@noaa.gov) ####
####              shared on March 6, 2020                    ####
#################################################################

# create dataframe of individual locations
locdat <- capdat[,c("ID..","Long","Lat", "Area", "year")] %>% 
  mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
  filter(!is.na(Long)) %>%
  mutate(Long = {ifelse(.$Long > 0, (.$Long*-1), .$Long)}) %>% 
  group_by(Area,ID..) %>% 
  mutate(ny=n_distinct(year)) %>% 
  ungroup()

# function to calculate distance between points
dist.ind.mean <- function(x) {
  # filter multiyear individuals (n=76)
  my.iyl <- x %>% 
    filter(ny>1)                                    
  # calculate mean distance among within-individual sightings 
  d.ind <- my.iyl %>% 
    tidyr::nest(-ID..) %>% 
    mutate(dists = purrr::map(data, ~ geosphere::distm(data.frame(lon=.x$Long,lat=.x$Lat))/1000),
           xdist = purrr::map(dists, function(x) mean(x[upper.tri(x)])))   # mean per individual
    return(mean(unlist(d.ind$xdist)))   
}

# create looping vector and empty data objects
subset <- c("Kauai/Niihau", "Oahu", "Maui Nui", "Hawaii")
plot.list <- list()
wilcox.p <- vector()

# permute data 1000x, calculate distance between points
for (j in 1:length(subset)){
truth <- as.numeric(rep(NA, 1000))
test <- as.numeric(rep(NA, 1000))
set.seed(1020)
for (i in 1:length(test)) {
  all.iyl <- locdat[which(locdat$Area==subset[j]),] %>% 
    group_by(ID.., year) %>% sample_n(1) %>% ungroup() %>% arrange(year, ID..) 
  truth[i] <- dist.ind.mean(all.iyl)
  perm <- all.iyl %>% 
    group_by(year) %>% 
    mutate(r=sample(1:n()), ID..=ID..[r], ny=ny[r]) %>% 
    ungroup() %>% dplyr::select(-r)
  test[i] <- dist.ind.mean(perm)
}
rm(all.iyl, perm, i)

# test whether individual movement is significantly less than stock movement
permtest <- wilcox.test(truth,test, alternative = "less")
wilcox.p <- c(wilcox.p, permtest$p.value)

# combine permutation data for plotting
sight.dist <- rbind(cbind.data.frame(dist=test,type="test"),cbind.data.frame(dist=truth,type="truth"))

#color plot
plot.list[[j]] <- ggplot(sight.dist, aes(x=dist, fill=type)) +
  geom_histogram(alpha = 0.3, position = "identity", bins = 150) +
  scale_fill_manual(values = c("#404080", "#69b3a2"), labels = c("Simulated data", "Real data")) +
  theme_classic(18) +
  labs(fill="") +
  xlab("") +
  ylab("") +
  annotate("text", label = paste(ifelse(subset[j] == "Oahu", "O\u02BBahu",
                    ifelse(subset[j] == "Kauai/Niihau", "Kaua\u02BBi/Ni\u02BBihau",
                           ifelse(subset[j] == "Hawaii", "Hawai\u02BBi", subset[j]))),
             ", ",
             ifelse(permtest$p.value < 0.0001, "p < 0.0001",
                    paste("p = ", format(permtest$p.value, digits = 2), sep = "")), sep = ""), 
           x = Inf, y = Inf, hjust = 1, vjust = 1, size = 7) +
  coord_cartesian(ylim=c(0, 300)) +
  xlim(0,30) +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"), legend.text=element_text(size=20))

#bw plot
# plot.list[[j]] <- ggplot(sight.dist, aes(x=dist, fill=type)) +
#      geom_histogram(position = "identity", bins = 150) +
#      scale_fill_manual(values = c("black", "grey50"), labels = c("Simulated data", "Real data")) +
#   theme_bw(18) +
#   theme(axis.title.x = element_text(hjust = 1)) +
#   labs(fill="") +
#   annotate("text", label = paste(ifelse(subset[j] == "Oahu", "O\u02BBahu", 
#                     ifelse(subset[j] == "Kauai/Niihau", "Kaua\u02BBi/Ni\u02BBihau",
#                            ifelse(subset[j] == "Hawaii", "Hawai\u02BBi", subset[j]))),
#              ", ",
#              ifelse(permtest$p.value < 0.0001, "p < 0.0001", 
#                     paste("p = ", format(permtest$p.value, digits = 2), sep = "")), sep = ""), x = Inf, y = Inf, hjust = 1, vjust = 1) +
#   ylab("") +
#   xlab("") +
#   coord_cartesian(ylim=c(0, 300)) +
#   xlim(5,30)

}