## PLOTS sensitivity outputs

library(here)
library(abind)
library(tidyverse)
library(magrittr)
library(fitdistrplus)
library(gamlss.dist)
library(gamlss)
library(truncnorm)
library(RGeode)
library(ggpubr)
library(readr)
library(DT)
library(scales)
library(ggridges)
library(viridis)
library(wesanderson)
library(ggstar)
library(data.table)
library(shiny)
library(gridExtra)
library(patchwork)

##-----------------------------------##
## Load Outputs
##-----------------------------------##

mypath=paste(here(),"Sensitivity", sep="/")
all_files <- list.files(path=mypath)

outputs <- "ModelOutputs"  # keep space

results <- grep(outputs, all_files, value = TRUE)
results=results[-2]#remove main
results=paste(mypath,results,sep="/")

run_names=NA
# Loop through each RData file and load it into R
for (i in seq_along(results)) {
  file_name <- basename(results[i])  # Get the file name 
  list_name <- strsplit(file_name, "_")[[1]][1]
  run_names[i]<-list_name
  load(results[i])
  assign(list_name, model_results)
 
}

#-------------------------------------------------

density_runs=data.frame()

for (run in run_names) {
  current <- get(run, envir = .GlobalEnv)
    df1 <- as.data.frame(current[[1]][1])
    df2 <- as.data.frame(current[[2]][1])
    
    df <- bind_rows(df1, df2)
    df$Parameter <- run 
    density_runs<-bind_rows(density_runs,df)
    
    rm(df1,df2,df,current)
    
  }

names(density_runs)<-c("N","simul","year","recdat","scenario","Parameter")

simul=unique(last(density_runs$simul))

density_runs%<>%
  filter(year==last(year-1))%>%## compare outputs of the last simulated year
  group_by(Parameter, scenario, recdat)%>%
  summarise(density=mean(N/15))

##----------------------------------------------
## Set main as the reference and calculate changes 
##----------------------------------------------

density_change <- density_runs %>%
  group_by(scenario, recdat) %>%
  mutate(
    ref_value = density[Parameter == "main"],  # Reference model for relative change
    relative_change = (density - ref_value) / ref_value, 
  ) %>%
  ungroup() %>%
  dplyr::select(-ref_value) 

head(density_change)


##---------------------------
## Sensitivity Plot
##---------------------------

cols <-  wes_palette(6, 2, type = c("discrete"))

density_change$scenario <- factor(density_change$scenario)
density_change$scenario <- relevel(density_change$scenario, ref = "nonharvested")
levels(density_change$scenario)[1]<-"unharvested"


##Order Parameters
density_change$Parameter <- factor(density_change$Parameter)
density_change%<>%filter(!Parameter=="main")%>%droplevels()

density_change$effect <- str_extract(density_change$Parameter, "(?i)Minus|Plus")
density_change$Parameter<- as.factor(str_remove(density_change$Parameter, "Minus|Plus"))

print(levels(density_change$Parameter))

levels(density_change$Parameter)<-c("Adult Density", "Adult Size", "Regrowth", "Efficiency", "Frequency", 
                                    "Harvest size", "Growth rate", "Mortality", "Max age", "Neighbors", "Recruit mortality")

save(density_change, file="sensitivy to density. Rdata")

# Create the plot
print(levels(density_change$Parameter))

##position of text
label_unh<-levels(density_change$Parameter)[5]
unh<-levels(density_change$Parameter)[4]

label_harv<-levels(density_change$Parameter)[3]
harv<-levels(density_change$Parameter)[2]


density_change%>%
ggplot(aes(y =Parameter, x = log(relative_change+1), fill = scenario)) +
   geom_bar(aes(ifelse(effect=="Minus",log(relative_change+1),0)),
                   stat = "identity", position = "dodge",
                   alpha=0.3)+
  scale_fill_manual(values = rev(cols))+  
  geom_bar(aes(ifelse(effect=="Plus",log(relative_change+1),0)),
           stat = "identity", position = "dodge",alpha=1)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    y = "Parameter",
    x = "Relative Change (log+1)",
    fill = "Scenario"
  ) +
  geom_text(data=subset(density_change,recdat=="High"),aes(y = label_unh, x=1,label=scenario[1]), size = 2.6, color = cols[1])+
  geom_text(data=subset(density_change,recdat=="High"),aes(y = label_harv, x=1,label=scenario[2]), size = 2.6, color = cols[2])+
  geom_text(data=subset(density_change,recdat=="High" & scenario==scenario[1] & effect=="Minus"),aes(y = unh, x=0.7,label="-20%"), size = 2.6, color = cols[1], alpha=0.12)+
  geom_text(data=subset(density_change,recdat=="High" & scenario==scenario[1] & effect=="Minus"),aes(y = unh, x=1.3,label="+20%"), size = 2.6, color = cols[1], alpha=1)+
  geom_text(data=subset(density_change,recdat=="High" & scenario==scenario[1] & effect=="Minus"),aes(y = harv, x=0.7,label="-20%"), size = 2.6, color = cols[2], alpha=0.12)+
  geom_text(data=subset(density_change,recdat=="High" & scenario==scenario[1] & effect=="Minus"),aes(y = harv, x=1.3,label="+20%"), size = 2.6, color = cols[2], alpha=1)+
  
  facet_wrap(~recdat)+
  theme(legend.position = "none", panel.grid = element_blank(),panel.background = element_rect(fill="grey99"))
  
#ggsave(file=paste("Sensitivity/","Fig6.Sensitivity.png", sep=""), width = 5, height = 6, dpi=400)


