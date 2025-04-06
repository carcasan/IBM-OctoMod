###----------------------------------------------------------###
# Project: Individual-Based Model of gorgonian A. elisabethae ##
# to investigate impact of harvesting in the Bahamas          ##
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
library(shiny)
library(purrr)
library(DT)
library(here)

##-------------------------------NOTE:-------------------------------------------## 
## The following code was built to:                                              ##  
## Create "neighboring" Adult populations based on observed data                 ##
##-------------------------------------------------------------------------------## 

##Load saved VBcurve estimates
load("VB_CHparams.RData")

## wrap pop model in function
f_mod.neighbor <- function(A,recdat){
  ##Create Adult population based on obs across all sites 
  A=rDPO(10, mu = 9.29, sigma = 5.65)## Highly variable
  
  A=sample(A[A>0], 1, replace = TRUE)#Define neighbors density ensuring min 1 adult
  
  Agepop=rPO(A,mu=9.2)#Assign age
  
  ##Sample within the CI of the predicted VBgrowth Curve to give varying size per age
  LCI <- UCI <- LPI <-UPI <-numeric(length(Agepop))
  for (i in 1:length(Agepop)) {
    pv <- ests[,"Linf"]*(1-exp(-ests[,"K"]*(Agepop[i]-ests[,"t0"])))##predict for each "Adult" colony
    LCI[i] <- quantile(pv,0.025)#LowerCI
    UCI[i] <- quantile(pv,0.975)##UpperCI
    LPI[i] <- quantile(pv-bootTypical$rse,0.025)
    UPI[i] <- quantile(pv+bootTypical$rse,0.975)
  }
  
  Height=numeric(length(Agepop))
  Height=runif(length(Height), min=LPI, max=UPI)##Assign height
  
  Npop=data.table::data.table(ColonyNo=1:A,
                  Age=0, 
                  Height=0)
  
  Npop$Height=Height
  Npop$Age=Agepop
  
  ## Based on Harvest efficiency data
  ## Update neighbors according to recruitment scenario
  ## `high efficiency`=	0.54## proportion of colonies harvested--> low recruitment
  ## `low efficiency`=	0.23-->High recruitment
  
  L= round(A*0.54)
  H= round(A*0.23)
  
  ##Define recdat
  recdat <- "Highly variable"
  
  if(A>7){
  A <- switch(recdat,
              "Highly variable" = A,
              "Low" = L,
              "High" = H)

  }
  ##Update density
  Npop=Npop[1:A,]
  return(Npop)
}

#neighbors = 6
#recdat=recdat

f_rec.neighbor <- function(neighbors,recdat){
  NRec=list()
  for (i in 1:neighbors){
    NRec[[i]]<-f_mod.neighbor(A,recdat=recdat)
  }
  
  AA<-purrr::map(NRec, filter, Height>=20)##Select reproductive
  AA<-purrr::map(AA, mutate, Aa= sum(Height^2))
  AA<-purrr::map(AA, mutate, Rec= exp(-7.96)*Aa^1.19)##Based on log-log fitted model
  
  AA<-purrr::map(AA, mutate, Rec= 0.05*Rec)##Only 5% is exported
  AA=bind_rows(AA)
  
  Rec= ceiling(sum(unique(AA$Rec)))
  
  return(Rec)
}


#save.image(file = "Metapopulations.Rdata")

## Use for modelling focal population 

## --------------------------END --------------------------------##