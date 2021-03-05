library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(magrittr)
library(gridExtra)
library(ggpattern)
library(dplyr)
 
####   Test whether sampling is significantly non-random     ####
####             adapted by Amy M. Van Cise                  ####
#### from code by K. Alexandra Curtis (alex.curtis@noaa.gov) ####
####              shared on March 6, 2020                    ####
#################################################################

##function to get legend from ggplot list
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#### Test for non-random sampling within Areas###

# create looping vector and empty data objects
subset <- c("Kauai/Niihau", "Oahu", "Maui Nui", "Hawaii")
stockwilcox.p <- vector()
stockplot.list <- list()

for (j in 1:length(subset)){
  set.capdat <- capdat[which(capdat$Area == subset[j]),]
  # Define parameter values
  n.occasions <- length(unique(set.capdat$year))                   # Number of capture occasions, defined by the data
  n.original <- set.capdat %>% group_by(year) %>% filter(Sight.type == "Original") %>% 
    summarise(n = n()) %>% summarise(mean(n)) %>% .$"mean(n)" %>% round(.,digits = 0)
  marked <- rep(n.original, ifelse(n.occasions > 1, n.occasions-1, 1))   # Annual number of newly marked individuals, defined by the data
  phi <- rep(0.9, ifelse(n.occasions > 1, n.occasions-1, 1)) # survival (from Stolen and Barlow 2003)
  p <- rep(0.4, ifelse(n.occasions > 1, n.occasions-1, 1))  # capture probability (from Kery 2012)

  # Define matrices with survival and recapture probabilities
  PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
  P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

  # Define function to simulate a capture-history (CH) matrix
  simul.cjs <- function(PHI, P, marked){
    n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
 
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
        # Bernoulli trial: does individual survive occasion?
        sur <- rbinom(1, 1, PHI[i,t-1])
        if (sur==0) break		# If dead, move to next individual 
        # Bernoulli trial: is individual recaptured? 
        rp <- rbinom(1, 1, P[i,t-1])
        if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
  }
  
# Execute function
CH <- simul.cjs(PHI, P, marked)
mean(rowSums(CH))
#hist(rowSums(CH))

randsample <- cbind.data.frame(mean=rowSums(CH), type="random.sample")

# real data for annual resights
annual.resight <- capdat[which(capdat$Area == subset[j]),] %>% group_by(ID..,year) %>% tally() %>% summarize(mean = n()) %>% mutate(type = "annual.sight")

# test for differences from random
yeartest <- wilcox.test(annual.resight$mean,randsample$mean)
stockwilcox.p <- c(stockwilcox.p, yeartest$p.value)

# ggplot

#uniform font sizes
geom.text.size = 7
theme.size = (14/5) * geom.text.size


sight.hist <- rbind(annual.resight[c("mean","type")], randsample)
stockplot.list[[j]] <- ggplot(sight.hist, aes(mean, fill = type)) +
  geom_density(color="#e9ecef", alpha = 0.3, position = "identity") +
  scale_fill_manual(values=c("#69b3a2", "#404080"), labels = c("Real data", "Simulated data")) +
  theme_classic(theme.size) +
  labs(fill="")+
  xlab("") +
  ylab("") +
  annotate("text", label = paste(ifelse(subset[j] == "Oahu", "O\u02BBahu",
                    ifelse(subset[j] == "Kauai/Niihau", "Kaua\u02BBi/Ni\u02BBihau",
                           ifelse(subset[j] == "Hawaii", "Hawai\u02BBi", subset[j]))),
                           ", ",
             ifelse(yeartest$p.value < 0.0001, "p < 0.0001",
                    paste("p = ", format(yeartest$p.value, digits = 2), sep = "")), sep = ""), 
           x = Inf, y = Inf, hjust = 1, vjust = 1, size = geom.text.size) +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"), legend.text=element_text(size=20))

# same plot, in black and white 
# stockplot.list[[j]] <- ggplot(sight.hist, aes(mean)) +
#   geom_density_pattern(aes(pattern=as.factor(type)), alpha = 0.6, position = "identity", pattern_density = 0.6) +
#   scale_pattern_discrete(name = "", labels = c("Real dataset", "Simulated dataset")) +
#   theme_bw(18) +
#   theme(axis.title.x = element_text(hjust = 1)) +
#   labs(fill="")+
#   annotate("text", label = paste(ifelse(subset[j] == "Oahu", "O\u02BBahu", 
#                     ifelse(subset[j] == "Kauai/Niihau", "Kaua\u02BBi/Ni\u02BBihau",
#                            ifelse(subset[j] == "Hawaii", "Hawai\u02BBi", subset[j]))),
#              ", ",
#              ifelse(yeartest$p.value < 0.0001, "p < 0.0001", 
#                     paste("p = ", format(yeartest$p.value, digits = 2), sep = "")), sep = ""), x = Inf, y = Inf, hjust = 1, vjust = 1) +
#   ylab("") +
#   xlab("") +
#   theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"))

}


#### Test for non-random sampling within subreas###

# create looping vector and empty data objects
subarea <- c("KA","KB","KC","KD","OA","OB","MA","MB","HA","HB")
subwilcox.p <- vector()
subplot.list <- list()

for (j in 1:length(subarea)){
  set.capdat <- capdat[which(capdat$subarea == subarea[j]),] %>% 
      filter(!is.na(subarea))
    # Define parameter values
    n.occasions <- length(unique(set.capdat$year))                   # Number of capture occasions, defined by the data
    n.original <- set.capdat %>% group_by(year) %>% filter(Sight.type == "Original") %>% 
      summarise(n = n()) %>% summarise(mean(n)) %>% .$"mean(n)" %>%  round(.,digits = 0)
    marked <- rep(n.original, ifelse(n.occasions > 1, n.occasions-1, 1))   # Annual number of newly marked individuals, defined by the data
    phi <- rep(0.9, ifelse(n.occasions > 1, n.occasions-1, 1)) # survival (from Stolen and Barlow 2003)
    p <- rep(0.4, ifelse(n.occasions > 1, n.occasions-1, 1))# capture probability (from Kery 2012)
    
    # Define matrices with survival and recapture probabilities
    PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
    P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))
    
    # Define function to simulate a capture-history (CH) matrix
    simul.cjs <- function(PHI, P, marked){
      n.occasions <- dim(PHI)[2] + 1
      CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
      # Define a vector with the occasion of marking
      mark.occ <- rep(1:length(marked), marked[1:length(marked)])
      # Fill the CH matrix
      
      for (i in 1:sum(marked)){
        CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
        if (mark.occ[i]==n.occasions) next
        for (t in (mark.occ[i]+1):n.occasions){
          # Bernoulli trial: does individual survive occasion?
          sur <- rbinom(1, 1, PHI[i,t-1])
          if (sur==0) break		# If dead, move to next individual 
          # Bernoulli trial: is individual recaptured? 
          rp <- rbinom(1, 1, P[i,t-1])
          if (rp==1) CH[i,t] <- 1
        } #t
      } #i
      return(CH)
    }
    
    # Execute function
    CH <- simul.cjs(PHI, P, marked)
    mean(rowSums(CH))
    #hist(rowSums(CH))
    
    randsample <- cbind.data.frame(mean=rowSums(CH), type="random.sample")
    
    #real data for annual resights
    annual.resight <- capdat[which(capdat$subarea == subarea[j]),] %>% group_by(ID..,year) %>% 
      tally() %>% summarize(mean = n()) %>% mutate(type = "annual.sight")
    
    #test for differences from random
    yeartest <- wilcox.test(annual.resight$mean,randsample$mean)
    subwilcox.p <- c(subwilcox.p, yeartest$p.value)
    
    #ggplot
    sight.hist <- rbind(annual.resight[c("mean","type")], randsample)
    subplot.list[[j]] <- ggplot(sight.hist, aes(mean, fill = type)) +
      geom_density(color="#e9ecef", alpha = 0.4, position = "identity") +
      scale_fill_manual(values=c("#69b3a2", "#404080"), labels = c("Real data", "Simulated data")) +
      theme_classic(theme.size) +
      labs(fill="") +
      xlab("") +
      ylab("") +
      annotate("text", label = paste(subarea[j], ", ",
                 ifelse(yeartest$p.value < 0.0001, "p < 0.0001",
                        paste("p = ", format(yeartest$p.value, digits = 2), sep = "")), sep = ""), 
               x = Inf, y = Inf, hjust = 1, vjust = 1, size = geom.text.size) +
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"), legend.text=element_text(size=20))
    
 # same plot, in black and white 
    # subplot.list[[j]] <- ggplot(sight.hist, aes(mean)) +
    #   geom_density_pattern(aes(pattern=as.factor(type)), alpha = 0.6, position = "identity", pattern_density = 0.6) +
    #   theme_bw(18) +
    #   theme(axis.title.x = element_text(hjust = 1)) +
    #   labs(fill="")+
    #   annotate("text", label = paste(subarea[j], ", ",
    #              ifelse(yeartest$p.value < 0.0001, "p < 0.0001", 
    #                     paste("p = ", format(yeartest$p.value, digits = 2), sep = "")), sep = ""), x = Inf, y = Inf, hjust = 1, vjust = 1) +
    #   scale_pattern_discrete(name = "", labels = c("Real dataset", "Simulated dataset")) +
    #   theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"))

}
