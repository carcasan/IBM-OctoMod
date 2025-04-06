
##-----------------------------------------------------
# ##Set VB growth curve from CH data to estimate Height at Age
##-----------------------------------------------------
library(simpleboot)
library(car)  # for Boot()
library(FSA)  #Fr VB function
library(nlstools)
library(tidyverse)
library(readxl)
library(ggpmisc)
library(here)

##Download GR data
obsgrowth <- read_csv("CrossHarbour_growthdata.csv")


## vB formula
vbgf <- function(L_inf, k, age, t_0 = 0) {
  L_t <- L_inf * (1 - exp(-k * (age - t_0)))  
  return(L_t)
}

##Predict growth Rate age-based
growth_rate <- function(a) {
  return(k * L_inf * exp(-k * (a - t_0)))
}
#size-based
growth_rate2<-function(h){
  return((L_inf-h)*k)
}


set.seed(100)

##--------------------------------------------------------
##USE CURRENTLY FOR MODEL

vb<-vbFuns(param="Typical")# Synonymous of 'BevertonHolt' parameterization of the vonB function

### get starting values
f.starts<-vbStarts(Height ~ Agerings, data=obsgrowth,methLinf="Walford")

##update as per Goffredo
f.starts[1]= 81.7
f.starts[2]=0.041


##Fit model with starting values
fit1<-nls(Height ~ vb(Agerings,Linf,K,t0), data=obsgrowth, start=f.starts,
control = nls.control(maxiter = 2000, tol = 1e-6))
summary(fit1,correlation = TRUE)
AIC(fit1)
#residPlot(fit1)
hist(residuals(fit1),main="")## approx. normal 

## calculate with boostrapping for nls models
bootTypical <- nlsBoot(fit1,niter=2000)##Kept only that converged

ests <- bootTypical$coefboot
#plot(bootTypical)

#confint(bootTypical,plot=TRUE)
L_inf <- bootTypical$bootCI[1]
k <- bootTypical$bootCI[2]
t_0 <-bootTypical$bootCI[3]

# Create age vector (max 55 as per model)
Age <- seq(0, 55, by = 0.1)

# Predict height at age
`Height (cm)`<-vbgf(L_inf, k, Age + 1, t_0)

preds=data.frame(Age=Age,Height=`Height (cm)`)

preds$Height[preds$Age==55]## max length

# Get confidence intervals
LCI <- UCI <- LPI <- UPI <- numeric(length(Age))
for (i in 1:length(Age)) {
  pv <- ests[,"Linf"]*(1-exp(-ests[,"K"]*(Age[i]-ests[,"t0"])))
  LCI[i] <- quantile(pv,0.025)## CI
  UCI[i] <- quantile(pv,0.975)
  LPI[i] <- quantile(pv-bootTypical$rse,0.025)## prediction bounds (CI - se)
  UPI[i] <- quantile(pv+bootTypical$rse,0.975)## (CI + se)
}

#`Growth rate (cm year)`<-growth_rate(Age)
`Growth rate (cm year)`<-growth_rate2(`Height (cm)`)


preds$Growth<-`Growth rate (cm year)`


##-------------------------
##PLOT Supplementary Figure
##-------------------------

Linf="Linf="
K="K="
t0="t0="

#tiff("FigureS1.tiff", width = 2000, height = 1200, res = 300)

par(mfrow = c(1, 2))

plot(`Height (cm)`~Age,col="grey", type="l",cex=1)
points(obsgrowth$Height~obsgrowth$Agerings,pch=19,cex=0.4)
text(x=45,30, paste(Linf,round(L_inf,2)),cex=0.7)
text(x=45,20, paste(K,round(k,2)),cex=0.7)
text(x=45,10, paste(t0,round(t_0,2)),cex=0.7)
text(x=18,80, labels="predicted", col="grey")
text(x=18,5, labels="observed", col="black")
mtext("a)", side = 3, line = 0, at = -5, cex = 1.5)  

#lines(UCI~Age,type="l",col="blue",lwd=2,lty=2)
#lines(LCI~Age,type="l",col="blue",lwd=2,lty=2)
lines(UPI~Age,type="l",col="grey",lwd=2,lty=2)
lines(LPI~Age,type="l",col="grey",lwd=2,lty=2)

#plot(Age,`Growth rate (cm year)`, col="grey",type="l")
plot(`Height (cm)`,`Growth rate (cm year)`, col="grey", type="l",cex=1)
mtext("b)", side = 3, line = 0, at = -5, cex = 1.5)  

dev.off()

##-----------------------------------------------------
##Store outputs and save to use for initialization
##-----------------------------------------------------

# save(vbgf, L_inf, k, t_0, ests,bootTypical, growth_rate, growth_rate2,file="VB_CHparams.RData")


