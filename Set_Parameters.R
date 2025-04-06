## GATHER PARAMETERS DATA

library(abind)
library(gamlss.dist)
library(tidyverse)
library(magrittr)
library(fitdistrplus)
library(gamlss.dist)
library(gamlss)
library(abind)
library(truncnorm)
library(profvis)
library(here)

##-----------------------------------------
## Datasets
##-----------------------------------------
 
# ##Load VBcurve estimates and growth rate equations
load("VB_CHparams.RData")
 
##Probability of partial mortality (size-based)
partMort=data.frame(size=c("0-10","10-20","20-30","30-40","40-50",">50"),
                           part.mort=c(0.382, 0.238, 0.598, 0.612, 0.296,0.296))#Adjusted to 1year 

### Length of partial mortality
Partial_mortality <- read.csv("Partial_mortality.csv")

## Mortality rates
Mortality<-readxl::read_excel("Parameters.xlsx",sheet="Mortality")

##Load saved Initial Population (50 simuls)
load("Initialization.RData")

##Save data to run Mod

#save.image("Set_Lparameters.RData")
