rm(list=ls())

library(fitdistrplus)
library(gamlss.dist)
library(gamlss)
library(abind)
library(truncnorm)
library(magrittr)

##Load vB estimates
load("VB_CHparams.RData")

##-------------------------------------------------------------------------
##Function to generate the initial population of 1yr-old recruits and adults of different sizes
##-------------------------------------------------------------------------

isEmpty <- function(x) {##Detects if random number is zero to avoid errors
  return(length(x)==0)
}

# set.seed(1)


  inipop=function(nPop){
    ## Determine No. of recruits based on 2012 data
    R=rDPO(1, mu=137, sigma=116)##Data  extrapolated to 15m2
    ##Adult population based on 2012
    A=rPO(1, mu=15.75)##based on CH data set
    
    ##Generate random Ages based on obs Age-rings in CH only for >1yr colonies
    Agepop=rPO(A,mu=9.2)
    Hr=runif(R, 2,5)##Variate recruits size based on obs (3-5cm for ~13-14mo old)
    
    ##Sample within the CI of the predicted VBgrowth Curve to give varying size per age
    LCI <- UCI <- LPI <-UPI <-numeric(length(Agepop))
    for (i in 1:length(Agepop)) {
      pv <- ests[,"Linf"]*(1-exp(-ests[,"K"]*(Agepop[i]-ests[,"t0"])))##predict for each "Adult" colony
      LCI[i] <- quantile(pv,0.025)#LowerCI
      UCI[i] <- quantile(pv,0.975)##UpperCI
      ##approximate prediction bounds to fitted line
      ##by +/- the residual SE from each bootstrap model
      LPI[i] <- quantile(pv-bootTypical$rse,0.025)
      UPI[i] <- quantile(pv+bootTypical$rse,0.975)
    }
    
    Height=numeric(length(Agepop))
    for (i in 1:length(Height)) {
      Height[i]=runif(1, min=LPI[i], max=UPI[i])
      Height[i][isEmpty(Height[i])]= 0##Random number from within the CI of predicted lengths
    }
    
    summary(Height)
    
    ##-------------
    ##Create Data
    ##-------------
    
    ##Total Number of Individuals per 15m2 
    nPop=R+A
    
    PEpop=data.frame(ColonyNo=1:nPop,
                     state="A", #Start population with 1yr-old Recruits
                     Age=1, ##Recruits survived the first year
                     Height=0,##set a min size for recruits (3-5 cm for 13-14 month old--CR paper
                     stringsAsFactors = FALSE)
    PEpop$state[1:R]="R"
    PEpop$Height[PEpop$state=="A"]=Height
    PEpop$Age[PEpop$state=="A"]=Agepop
    if(R>0){PEpop$Height[PEpop$state=="R"]=Hr}# avoid negative preds when using LPI
    PEpop%<>%dplyr::select(-state)
    return(PEpop)
  }



repPop=list()
for (i in 1:50){
repPop[[i]]<-inipop(nPop)
}

##-----------------------------------------------------
##Store outputs and save to use for initialization in tool
##All simulations start with same population  in year 1
##-----------------------------------------------------

#save(repPop, file="Initialization.RData")

