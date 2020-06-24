library(RMark)
library(lubridate)
library(reshape2)
library(magrittr)
library(tidyverse)
library(FSA)
library(dplyr)

####POPAN abundance estimation and parameter summary tables####
####            Written by Amy M. Van Cise                 ####
###############################################################

#source code to generatue POPAN input and calculated proportion distinctive animals
setwd("C:/Users/Amy/Google Drive/00 CRC/01 Tursiops abundance")
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

#for loop to estimate abundance for each stock
for (i in 1:length(subset)){
  #subset data
  set.capdat <- capdat[which(capdat$Area == subset[i]),] %>% 
    filter(!is.na(subarea)) %>% 
    mutate(Long=as.numeric(Long), Lat = as.numeric(Lat)) %>% 
    mutate(subarea = ifelse(.$Area == "Maui Nui" & 
                    .$subarea == "OB", ifelse(.$Lat >= 21 & abs(.$Long) > 156.8 | abs(.$Long) >= 157.5, "MA","MB"),subarea))
  time.int <- diff(sort(unique(set.capdat$year)))
  #source("subarea effort estimation.R")
    
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
    #p.effort=list(formula=~prop.effort)
    #p.timeeffort=list(formula=~prop.effort*time)
  
    
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
    filter(weight > 0) %>%
    mutate(chat = chat)
  all_model.results <- rbind(all_model.results, setNames(Tt.set.results.table, names(all_model.results)))
  
  #return model estimates of annual N, corrected for pd
  pd <- mean(prop[grep(substr(subset[i],1,4),prop$Area),"mean"]$mean) #calculated within each Area in Prop.distinctive.R
  pdvar <- mean(prop[grep(substr(subset[i],1,4),prop$Area),"var"]$var,na.rm=TRUE)
  all_pd <- rbind(all_pd,c(subset[i],pd,pdvar), stringsAsFactors=FALSE)
  colnames(all_pd) <- c('Stock','$\\theta$', '$\\theta_{var}$')

  Tt.set.est <- popan.derived(obs.data.proc,Tt.set.results)
  Tt.est.corr <- Tt.set.est$N %>%
    mutate(N_i = .$N/as.numeric(pd), 
           year = rep(colnames(tempdat),length(unique(set.capdat$subarea))), 
           Area = subset[i]) %>%
    #estimate variance in N_i using delta method (from ALB 2018)
    mutate(Varcond = (se^2/pd^2)+(pdvar*(N^2/pd^4))) %>% 
    mutate(SEcond = sqrt(Varcond)) %>% 
    mutate(CVcond = (SEcond/N_i)) %>% 
    mutate(C_Burn = exp(1.96*sqrt(log(1+CVcond^2)))) %>% 
    mutate(LCL = round(N_i/C_Burn,digits=0)) %>%
    mutate(UCL =  round(N_i*C_Burn,digits=0)) %>%
    group_by(Area, year) %>% 
    summarize(N_i = round(sum(N_i),digits=0), se = sum(SEcond), CV = sum(CVcond), C_Burn = sum(C_Burn), LCL = sum(LCL), UCL = sum(UCL)) %>% 
    ungroup() %>% 
    dplyr::select(-C_Burn) %>% 
    rename(Stock = Area)
  all_total_abundance <- rbind(all_total_abundance,Tt.est.corr)
  
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