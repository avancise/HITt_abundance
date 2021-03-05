library(RMark)
library(lubridate)
library(reshape2)
library(magrittr)
library(tidyverse)
library(FSA)
library(dplyr)
library(mvtnorm)

####POPAN abundance estimation and parameter summary tables####
####            Written by Amy M. Van Cise                 ####
###############################################################

#must install MARK separately to run this code

#source code to generatue POPAN input and calculated proportion distinctive animals
setwd("C:/Users/Amy/Google Drive/05 CRC/01 Tursiops abundance")
source("RMark_abundance_datainput.R")
source("Prop_distinctive.R")

#Designate stock subsetting
#subset <- unique(capdat$Area)
subset <- c("Kauai/Niihau", "Oahu", "Maui Nui", "Hawaii")

#create empty dataframes
all_total_abundance <- data.frame()
all_model.results <- data.frame(p = numeric(), pent = numeric(), npar = numeric(), AICc = numeric(), 
                                Delta = numeric(), weight = numeric(), Deviance = numeric(), Stock = character(), chat = numeric())
all_model_phi <- data.frame()
all_model_p <- data.frame()
all_pd <- data.frame()
ests <- list()
sigma <- list()

#for loop to estimate abundance for each stock
for (i in 1:length(subset)){
  #subset data
  set.capdat <- capdat[which(capdat$Area == subset[i]),] %>% 
    filter(!is.na(subarea)) %>% 
    mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
    #mutate(subarea = {ifelse(.$Area == "Maui Nui" & .$subarea == "OB", ifelse(.$Lat >= 21 & abs(.$Long) > 156.8 | abs(.$Long) >= 157.5, "MA","MB"),subarea)}) 
    filter(case_when(Area == "Maui Nui" ~ subarea != "OB", 
                     Area == "Oahu" ~ !grepl("M",subarea),
                     Area == "Hawaii" ~ subarea != ".",
                    Area == "Kauai/Niihau" ~ subarea != "."))
  time.int <- diff(sort(unique(set.capdat$year)))
  
  #format data for RMark ch with year as the sampling period
  tempdat <- set.capdat %>% 
    melt(id.var=c("ID..","year"), measure.var="detect") %>%
    dcast(ID.. ~ year) %>% 
    `rownames<-`(.[,1]) %>% 
    dplyr::select(-ID..) %>% 
    replace(., .>0, 1) 
  
  #create capt hist object
  set.capt.hist <- capHistConvert(tempdat,in.type="individual",out.type="RMark")
  set.capt.hist <- merge(set.capt.hist, unique(set.capdat[,c("ID..","Area","subarea")]), by.x=0, by.y="ID..")
  
  
  #initialize model variables, set time intervals and factor covariates
  obs.data.proc <- process.data(set.capt.hist, model="POPAN", begin.time=min(set.capdat$year),
                                time.intervals = time.int, groups = "subarea")
  
  #test model goodness of fit
  gof<-release.gof(obs.data.proc)
  chat<-gof$Chi.square[3]/gof$df[3]
  
  
  #create PIM for MARK input
  Tt.set.ddl<-make.design.data(obs.data.proc)
  
  #incorporate subarea effort into design data
  #Tt.set.ddl$p <- merge_design.covariates(Tt.set.ddl$p, effort, bygroup = TRUE, bytime = TRUE)
  
  #Function for running Phi, p, pent, and N models
  Tt.all.models<-function() {
    Phi.dot=list(formula=~1)
    
    p.dot=list(formula=~1)
    p.time=list(formula=~time)
    p.subarea=list(formula=~subarea)
    p.timesubarea=list(formula=~subarea*time)
    
    pent.dot=list(formula=~1)
    pent.time=list(formula=~time)
    pent.subarea=list(formula=~subarea)
    pent.timesubarea=list(formula=~subarea*time)
    
    N.subarea=list(formula=~subarea)
    
    cml=create.model.list(model="POPAN")
    results=mark.wrapper(cml,data=obs.data.proc,ddl=Tt.set.ddl, delete = TRUE, output = FALSE)
    
    return(results)
  }
  
  #run model
  Tt.set.results <- Tt.all.models()
  
  #return model AIC, adjust for chat
  if (chat > 1) {Tt.set.results <- adjust.chat(chat=chat,Tt.set.results)}
  Tt.set.results.table <- Tt.set.results %>%
    model.table() %>% 
    dplyr::select(p, pent, npar, contains("AICc"), contains("Delta"), weight, contains("Deviance")) %>% 
    mutate(Area=subset[i]) %>% 
    filter(weight > 0.0001) %>%
    mutate(chat = chat)
  all_model.results <- rbind(all_model.results, setNames(Tt.set.results.table, names(all_model.results)))
  
  #return model estimates of annual N, corrected for pd
  pd <- mean(prop[grep(substr(subset[i],1,4),prop$Area),"mean"]$mean) #calculated within each Area in Prop.distinctive.R
  pdvar <- mean(prop[grep(substr(subset[i],1,4),prop$Area),"var"]$var,na.rm=TRUE)
  all_pd <- rbind(all_pd,c(subset[i],pd,pdvar), stringsAsFactors=FALSE)
  colnames(all_pd) <- c('Stock','$\\theta$', '$\\theta_{var}$')

  Tt.set.est <- popan.derived(obs.data.proc,Tt.set.results)
  Tt.est.corr <- Tt.set.est$Nbyocc %>%
    mutate(N_i = .$N/as.numeric(pd),
           year = colnames(tempdat),
           Area = subset[i]) %>%
    #estimate variance in N_i using delta method (from ALB 2018)
    mutate(Varcond = (se^2/pd^2)+(pdvar*(N^2/pd^4))) %>%
    mutate(SEcond = sqrt(Varcond)) %>%
    mutate(CVcond = (SEcond/N_i)) %>%
    mutate(C_Burn = exp(1.96*sqrt(log(1+CVcond^2)))) %>%
    mutate(LCL = round(N_i/C_Burn,digits=0)) %>%
    mutate(UCL =  round(N_i*C_Burn,digits=0)) %>%
    dplyr::select(-C_Burn) %>%
    rename(Stock = Area)

  all_total_abundance <- rbind(all_total_abundance,Tt.est.corr)
  
  #save abundance estimates and sigma for later analysis of trends
  ests[[i]] <- Tt.set.est$Nbyocc$N
  sigma[[i]] <- Tt.set.est$Nbyocc.vcv
  
  #return model-averaged phi
  MA_phi<-model.average(Tt.set.results,"Phi", drop=TRUE)[1,] %>%
    rownames_to_column() %>% 
    dplyr::rename(Parameter = rowname) %>% 
    mutate(Stock=subset[i]) %>% 
    dplyr::select(Stock, estimate, se)
  all_model_phi <- rbind(all_model_phi,MA_phi)
  
  #return model-averaged p
  MA_p<-model.average(Tt.set.results, "p", drop=TRUE) %>% 
    rownames_to_column() %>% 
    dplyr::rename(Parameter = rowname) %>% 
    dplyr::select(Parameter, estimate, se) %>% 
    mutate(Stock = subset[i])
  all_model_p <- rbind(all_model_p, MA_p)
  
}

write.csv(all_model_p, paste("model.average.p.",Sys.Date(),".csv"))

all_model.results$weight <- ifelse(all_model.results$weight >= 0.0001, round(all_model.results$weight, digits = 4), "< 0.0001")
all_pd$'$\\theta$' <- as.numeric(all_pd$'$\\theta$')
all_pd$'$\\theta_{var}$' <- as.numeric(all_pd$'$\\theta_{var}$')

#add number of sightings and individuals 
annual_sight <- capdat %>% 
  group_by(Area, year) %>% 
  summarize(Sightings = n(), Individuals = length(unique(ID..)))

all_total_abundance <- all_total_abundance %>% 
  mutate(year = as.numeric(year)) %>% 
  left_join(annual_sight, by = c("Stock" = "Area","year"))

#calculate annual growth
ann.growth <- all_total_abundance %>%
  mutate(year = as.numeric(year)) %>% 
  group_by(Stock) %>% 
  arrange(year) %>%
  mutate(Diff_year = year - lag(year),  # Difference in time (just in case there are gaps)
         Diff_growth = N_i - lag(N_i), # Difference in route between years
         Rate_percent = (Diff_growth / Diff_year)/N_i * 100)  # growth rate in percent 
  
mean_growth <- ann.growth %>%
  summarize(Mean_growth = mean(Rate_percent, na.rm = TRUE))

#calculate CI of slope for each stock
nrep <- 10000
stock_slopes <- matrix(0, 4, nrep)
pop_slope <- rep(0,nrep)
stock_intercept <- matrix(0,4,nrep)
pop_intercept <- rep(0,nrep)

for (j in 1:nrep){
  N <- as.data.frame(seq(2000,2018,by=1))
  colnames(N) <- "years"
  for (i in 1:length(subset)){
    years <- all_total_abundance %>% filter(Stock == subset[i]) %>% .$year
    N_istock <- rmvnorm(1,ests[[i]],sigma[[i]])[1,]
    stock_slopes[i,j] <- coefficients(lm(N_istock ~ years,2))[2]
    stock_intercept[i,j] <- coefficients(lm(N_istock ~ years,2))[1]
    N_temp <- cbind.data.frame(N_istock,years)
    stock <- subset[i]
    N <- N %>% left_join(N_temp, by = "years") %>% rename(!!stock := N_istock)
  }
  Ntot <- N %>% drop_na() %>% mutate(Total = rowSums(.[2:5]))
  pop_slope[j] <- coef(lm(Ntot$Total ~ Ntot$years))[2]
  pop_intercept[j] <- coef(lm(Ntot$Total ~ Ntot$years))[1]
}

rownames(stock_slopes) <- subset
pop_slope_ci <- quantile(pop_slope,probs = c(0.025,0.975))
pop_intercept_ci <- quantile(pop_intercept, probs = c(0.025, 0.975))
stock_slopes_ci <- matrix(0,4,2)
stock_intercept_ci <- matrix(0,4,2)
for (i in 1:4){
  stock_slopes_ci[i,1:2] <- quantile(stock_slopes[i,],probs = c(0.025,0.975))
  stock_intercept_ci[i,1:2] <- quantile(stock_intercept[i,],probs = c(0.025,0.975))
}
rownames(stock_slopes_ci) <- subset
rownames(stock_intercept_ci) <- subset
