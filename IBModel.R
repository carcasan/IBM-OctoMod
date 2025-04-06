###----------------------------------------------------------###
# Project: Individual-Based Model of gorgonian A. elisabethae ##
# to investigate demographic impact of harvesting in          ##
# the Bahamas                                                 ##
# Build Population model                                      ##
# Developer: Carolina Castro-Sanguino & Howard R. Lasker      ##
###----------------------------------------------------------###


##----------------------------------------
##Packages/libraries needed to run the code
##----------------------------------------
library(RGeode)
library(abind)
library(gamlss.dist)
library(tidyverse)
library(magrittr)
library(fitdistrplus)
library(gamlss)
library(abind)
library(truncnorm)
library(purrr)
library(here)
library(data.table)
library(parallel)
library(purrr)
library(DT)
library(here)
library(tictoc)

required_packages <- .packages()
#save(required_packages, file="mypackages.Rdata")

#source("Set_Parameters.R")## 
load("Set_Parameters.RData")
load("Metapopulations.Rdata")
source("Recruitment.R")
source("Run_Main.R")

##-------------------------------------
##Function to Model Pop Dynamics of AE
##-------------------------------------

  AEmod <- function(simul,recdat){
  simul=Simul  
  recdat=recdat 

  ##Choose Mortality rates  
  tokill=Mortality[c(1,M)]## M for specific scenario
  names(tokill)[2]<-"pk"
  tokill$class<-c(0,5,10,20,30,40,50) # create size classes for total mortality
  partMort$class <- c(0, 10, 20, 30, 40, 50) # create size classes for partial mortality

  
  harvestpop=list()
  
  # Pre-allocate output data tables
  Yield <- data.table::data.table(yield = integer(0), simul = integer(0), year = integer(0))
  countharvested <- data.table::data.table(count = integer(0), simul = integer(0), year = integer(0))
  countsize <- data.table::data.table(class = integer(0), N = integer(0), simul = integer(0), year = integer(0), scenario = NA)
  countage <- data.table::data.table(Age = integer(0), N = integer(0), simul = integer(0), year = integer(0), scenario = NA)
  ageheight <- data.table::data.table(Age = integer(0), height = integer(0), simul = integer(0), year = integer(0), scenario = NA)
  countall <- data.table::data.table(N = integer(0), simul = integer(0), year = integer(0), recdat = NA, scenario = NA)
  
    
    for (s in 1:length(scenarios)) {
      scenario <- scenarios[s]
      print(paste(scenario,recdat))
      
      for (i in 1:simul) {
          
          if (scenario=="nonharvested"){
            
            newpop=data.table::data.table(repPop[[i]]) 
            
            newRpop=data.table::data.table()
            
            for(t in 2:length(years)){
              
              print(paste(t,i,s))
              
              if (t==firstharvest-1){
                harvestpop[[i]]=newpop##save before it grows to avoid mismatch in pop state between scenarios
              }
              
              # Colonies grow  and age
              newpop[, Height :=Height+growth_rate(Age)]
              newpop[, Age :=Age+1]        
              
              # Define size classes for total mortality
              classes <- c(0, 5, 10, 20, 30, 40, 50)
              
              newpop$class <- findInterval(newpop$Height, classes)
              
              IDtokill <- integer(0)  # Initialize empty vector to store IDs to kill
              
              for(nk in 1:nrow(tokill)) {
                # Find indices in the current size class
                idx <- which(newpop$class == nk)
                
                # Determine how many to kill based on the mortality probability (pk)
                n_to_kill <- round(length(idx) * tokill$pk[nk])
                
                # Randomly sample to kill (if there are any in the class)
                if (n_to_kill > 0) {
                  sampled_colonies <- sample(idx, n_to_kill)
                  IDtokill <- c(IDtokill, sampled_colonies)
                }
              }
              
              newpop[, id := .I] 
              newpop = newpop[!id %in% IDtokill] # remove killed colonies
              newpop = droplevels(newpop)
              
              
              # Define size classes for partial mortality
              classes <- c(0, 10, 20, 30, 40, 50)
              
              newpop$class <- findInterval(newpop$Height, classes)
              
              IDtoshrink=integer(0) 
              for(m in 1:nrow(partMort)) {
                idx <- newpop$ColonyNo[newpop$class == m] # ID colony
                
                n_to_shrink <- round(length(idx) * partMort$part.mort[m])  # Number to shrink
                if (n_to_shrink > 0) {
                  sampled_colonies <- sample(idx, n_to_shrink)
                  IDtoshrink <- c(IDtoshrink, sampled_colonies)
                }
              }
              
              
              shrink<-sample(Partial_mortality$GR_cm_year, length(IDtoshrink),replace=TRUE)
              
              if (length(IDtoshrink) > 0) {
                n <- max(length(IDtoshrink), length(newpop$ColonyNo))
                length(IDtoshrink) <- n  # Extend length if necessary
                IDtoshrink[is.na(IDtoshrink)] <- 0  # Replace NAs with 0
                
                length(shrink) <- n  # Extend shrinkage vector if necessary
                shrink[is.na(shrink)] <- 0  # Replace NAs with 0
                
                # shrink selected colonies
                toshrink <- data.frame(IDtoshrink, shrink)
                newpop$Height[newpop$ColonyNo %in% toshrink$IDtoshrink] <- 
                newpop$Height[newpop$ColonyNo %in% toshrink$IDtoshrink] + toshrink$shrink[toshrink$IDtoshrink > 0]
              }
              
              newpop[Height < 0, Height := 1]# shrinking does not lead to mortality
              
              newpop=newpop[, .(ColonyNo, Age, Height)]
              
              AA = newpop[Height >= adultsize, sum(Height^2)] # Adult density threshold
              
              ##density-dependent recruitment
              rec.current = exp(-7.96) * AA ^ 1.19#Based on log-log fitted model
              ##Retention: Proportion that successfully recruit from local supply
              rec.current <- ifelse(recdat == "High", rec.current * 0.6, rec.current * 0.3)
              
              rec.current=ifelse(rec.current<0,0,rec.current)##Set to 0 if negative
              
              rec.neighb = f_rec.neighbor(neighbors = neighbors,recdat=recdat) #External recruits
              ##Total supply
              rec = ceiling(sum(c(rec.current, rec.neighb)))##retention+ migration
              
              ##Succesful recruits
              settlers=settlers # considers all recruit
              rec=ceiling(rec*settlers)
              
              newRpop=new.dp.recruits(rec)## 1year-olds
              
              nmax=max(newpop$ColonyNo)
              
              newRpop[, ColonyNo := nmax+ColonyNo] # Update pop count
              
              newpop=rbind(newpop,newRpop) # Add recruits to population
              
              newpop=newpop[!Age >MaxAge] # Senescence
              
              ## Save pop density outputs
              temp_dt1 <- newpop[, .(N = .N)]
              temp_dt1 <- temp_dt1[, `:=` (simul=i,year=t,recdat=recdat,scenario=scenario)]
              countall <- rbindlist(list(countall, temp_dt1), use.names = TRUE, fill = TRUE)
              
              # Define size classes
              classes <- c(0, 5, 10, 20, 30, 40, 50)
              newpop[, class := findInterval(Height, classes)]
              
              
              ## Save size structure 
              temp_dt2 <- newpop[, .N, by= .(class)]
              temp_dt2 <- temp_dt2[, `:=` (simul=i,year=t,scenario=scenario)]
              countsize=rbindlist(list(countsize, temp_dt2), use.names = TRUE, fill = TRUE)
              
              ## Save Age structure
              temp_dt3 <- newpop[, .N, by= .(Age)]
              temp_dt3 <- temp_dt3[, `:=` (simul=i,year=t,scenario=scenario)]
              countage=rbindlist(list(countage, temp_dt3), use.names = TRUE, fill = TRUE)
              
              ## Save average height per Age              
              temp_dt4 <- newpop[, .(height = mean(Height)), by = Age]
              temp_dt4 <- temp_dt4[, `:=` (simul=i,year=t,scenario=scenario)]
              ageheight=rbindlist(list(ageheight, temp_dt4), use.names = TRUE, fill = TRUE)
              
              
              t=t+1
            } 
            harvestpop = harvestpop[sapply(harvestpop, function(x) !is.null(x) && dim(x)[1] > 0)]
            
          }
          
          else if (scenario=="harvested"){
            
            yield=NA
            
            th=firstharvest-1
            
            newpop<-harvestpop[[i]] #Harvest starts from unharvested pop at year th
            newpop[, harvyr := -1 ] #Track which colonies are harvested
            newpop[, fate := 0 ]
  
           
             newRpop=data.table::data.table()
            
            for (t in th:length(years)){
              IDtoharvest=numeric()
              
              if(t>lastharvest+2){
                break
              }
              
              print(paste(t,i,s))
              
              # Colonies grow  and age
              newpop<-newpop[, Height :=Height+growth_rate(Age)]
              newpop<-newpop[, Age :=Age+1]        
              
              newpop[harvyr %in% 0:3, Height := Height * CG]##compensatory growth for 4 years after harvested
              
              
              # Define size classes
              classes <- c(0, 5, 10, 20, 30, 40, 50)
              
              newpop$class <- findInterval(newpop$Height, classes)
              
              IDtokill <- integer(0)  # Initialize empty vector to store IDs to kill
              
              for(nk in 1:nrow(tokill)) {
                idx <- which(newpop$class == nk)
                
                n_to_kill <- round(length(idx) * tokill$pk[nk])
                
                if (n_to_kill > 0) {
                  sampled_colony <- sample(idx, n_to_kill)
                  IDtokill <- c(IDtokill, sampled_colony)
                }
              }
               
              newpop[, id := .I] 
              newpop = newpop[!id %in% IDtokill]
              newpop = droplevels(newpop)
              
              # Define size classes
              classes <- c(0, 10, 20, 30, 40, 50)
              
              newpop$class <- findInterval(newpop$Height, classes)
              
              IDtoshrink=integer(0) 
              for(m in 1:nrow(partMort)) {
                idx <- newpop$ColonyNo[newpop$class == m] 
 
                n_to_shrink <- round(length(idx) * partMort$part.mort[m])
                if (n_to_shrink > 0) {
                  sampled_colonies <- sample(idx, n_to_shrink)
                  IDtoshrink <- c(IDtoshrink, sampled_colonies)
                }
              }
              
              shrink<-sample(Partial_mortality$GR_cm_year, length(IDtoshrink),replace=TRUE)
              
              if (length(IDtoshrink) > 0) {
                n <- max(length(IDtoshrink), length(newpop$ColonyNo))
                length(IDtoshrink) <- n  
                IDtoshrink[is.na(IDtoshrink)] <- 0  
                
                length(shrink) <- n  
                shrink[is.na(shrink)] <- 0
                
                toshrink <- data.frame(IDtoshrink, shrink)
                newpop$Height[newpop$ColonyNo %in% toshrink$IDtoshrink] <- 
                  newpop$Height[newpop$ColonyNo %in% toshrink$IDtoshrink] + toshrink$shrink[toshrink$IDtoshrink > 0]
              }
              
              newpop[Height < 0, Height := 1]
              
              newpop<-newpop[, .(ColonyNo, Age, Height, harvyr,fate)]
              
              AA = newpop[Height >= adultsize, sum(Height^2)]
              
              rec.current = exp(-7.96) * AA ^ 1.19
              rec.current <- ifelse(recdat == "High", rec.current * 0.6, rec.current * 0.3)
              rec.current=ifelse(rec.current<0,0,rec.current)
              
              rec.neighb = f_rec.neighbor(neighbors=neighbors,recdat=recdat) 
              rec = ceiling(sum(c(rec.current, rec.neighb)))
              
              settlers=settlers
              
              rec=rec*settlers
              
              newRpop=new.dp.recruits(rec)
              
              nmax=max(newpop$ColonyNo)
              
              newRpop[, ColonyNo := nmax+ColonyNo]
              
              newRpop[, harvyr := -1]
              newRpop[, fate := 0]
              
              newpop=rbind(newpop,newRpop)
              
              newpop=newpop[!Age >MaxAge]
              
              
              if(t %in% harvestyears){
                adultdens <- newpop[Height >= adultsize, .N]
                if(adultdens >= minadultpop){
                  dens=newpop[Height >harv.colony, .N]
                  
                  toharvest=efficiency 
                  
                  colhv <- newpop[Height > harv.colony, ColonyNo]
                  
                  if (length(colhv)>0) {
                    ##Biomass available before harvest
                    avail.biomass<-newpop[Height > harv.colony, Height]
                    ibiomass=sum((0.002*avail.biomass^2.075))
                    
                    samp<-colhv[sample(round(1:length(colhv)*toharvest),replace=TRUE)]
                    IDtoharvest=c(IDtoharvest,samp)
                  }
                  
                  trackcol= list(length(IDtoharvest),i, t)
                  # Check for missing values in x
                  if (any(is.na(trackcol))) {
                    warning("Missing values found in x! Skipping this iteration.")
                    next  # Skip this iteration if x contains NA
                  }
                  countharvested=rbind(countharvested, trackcol,fill=TRUE)
                  
                  if(is_empty(IDtoharvest)==FALSE){
                    cut<-base::sample(rs.height, length(IDtoharvest),replace=TRUE)
                    
                    for(e in 1: length(IDtoharvest)){
                      newpop$Height[newpop$ColonyNo == IDtoharvest[e]]=cut[e]
                    }
                  }
                  
                  
                  newpop[ColonyNo %in% IDtoharvest, fate := 1]
                  
                  ##Biomass left after harvest
                  resid.biomass<-newpop[fate == 1, Height]
                  
                  fbiomass=sum(0.002*(resid.biomass^2.075))
                  
                  yield=ibiomass-fbiomass 
                  temp_dt1 <- list(yield=yield,simul=i,year=t)  
                  Yield<- rbindlist(list(Yield, temp_dt1), use.names = TRUE, fill = TRUE)
                  
                }
              }
              
              ## Track harvested colonies for compernsatory growth
              newpop[harvyr == 0, harvyr := 1]  # Set h = 1 where h == 0
              newpop[ColonyNo %in% IDtoharvest & t %in% harvestyears, harvyr := 0]  # Set h = 0 when conditions are met
              newpop[harvyr > 0, harvyr := harvyr + 1]  # Increment h by 1 where h > 0
              newpop[harvyr < 0, harvyr := -1]  # Set h = -1 where h < 0
              
              temp_dt1 <- newpop[, .(N = .N)]
              temp_dt1 <- temp_dt1[, `:=` (simul=i,year=t,recdat=recdat,scenario=scenario)]
              countall <- rbindlist(list(countall, temp_dt1), use.names = TRUE, fill = TRUE)
              
              # Define size classes
              classes <- c(0, 5, 10, 20, 30, 40, 50)
              newpop[, class := findInterval(Height, classes)]
              
              ## Save outputs
              temp_dt2 <- newpop[, .N, by= .(class)]
              temp_dt2 <- temp_dt2[, `:=` (simul=i,year=t,scenario=scenario)]
              
              countsize=rbindlist(list(countsize, temp_dt2), use.names = TRUE, fill = TRUE)
              
              temp_dt3 <- newpop[, .N, by= .(Age)]
              temp_dt3 <- temp_dt3[, `:=` (simul=i,year=t,scenario=scenario)]
              
              
              countage=rbindlist(list(countage, temp_dt3), use.names = TRUE, fill = TRUE)
              
              
              temp_dt4 <- newpop[, .(height = mean(Height)), by = Age]
              temp_dt4 <- temp_dt4[, `:=` (simul=i,year=t,scenario=scenario)]
              
              
              ageheight=rbindlist(list(ageheight, temp_dt4), use.names = TRUE, fill = TRUE)
              
              t=t+1
            }
            
            countall <- countall[!is.na(scenario)]
            countage<-countage[!is.na(scenario)]
            countsize<-countsize[!is.na(scenario)]
            ageheight<-ageheight[!is.na(scenario)]
            countharvested<-countharvested[simul>0]
            Yield<-Yield[yield>0]  
          }
        }
      }
    
      outputs<-return(list(countall = countall, countsize = countsize, 
                           countage = countage,ageheight=ageheight,
                           countharvested=countharvested,
                           Yield=Yield))
      return(outputs)
}
 
#save.image(file="model_data.Rdata")
  