###----------------------------------------------------------###
# Project: Individual-Based Model of gorgonian A. elisabethae ##
# to investigate impact of harvesting in the Bahamas          ##
# Developer: Carolina Castro-Sanguino                         ##
###----------------------------------------------------------###

##------------------
## install Libraries
##------------------

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


myplots<-"Plots/"

##-----------------------------------##
## Load Outputs
##-----------------------------------##

all_files <- list.files()

outputs <- "main "  # keep space


# Filter files that contain the string
results <- grep(outputs, all_files, value = TRUE)


# Loop through each RData file and load it into R
for (i in seq_along(results)) {
  load(results[i])
  }
# ##Load the parameters used in the simulation  
for (col_name in colnames(this_run)) {
  value <- this_run[[col_name]]
  assign(col_name, value, envir = .GlobalEnv)
}

scenarios=strsplit(this_run$scenarios, "_")[[1]]
harvestyears=seq(firstharvest,lastharvest,frequency)


##---------------------------##
## Create figures for paper
##---------------------------##

# Check if the folder exists, and if not, create it

folder <- paste(myplots,Sys.Date(), sep="")  # Modify this path

if (!dir.exists(folder)) {
  dir.create(folder)
  print(paste("Folder created at:", folder))
} else {
  print(paste("Folder already exists at:", folder))
}

##-------------------------------    
## Plot Simulations-- pop density
##-------------------------------
scenarios[1]<- "unharvested"

cols <-  wes_palette(6, 2, type = c("discrete"))

plot_env <- new.env()## save all plots here

# Function to create and combine plots for each scenario
pop_plots <- function(results, names) {
  
  for (i in seq_along(model_results)){
    df <- model_results[[i]][[1]]
    simul=max(df$simul)
    recdat=unique(df$recdat)
    
    df$scenario <- factor(df$scenario)
    df$scenario <- relevel(df$scenario, ref = "nonharvested")
    
    ##start from non-harvested state
    df$year[df$scenario=="harvested"] <- df$year[df$scenario=="harvested"]-1 # Years from first harvest
    df%<>%filter(year > firstharvest-1) ##remove burning period (low rec likely taking longer to equilibrium?)
    
    levels(df$scenario)[1]<-"unharvested"
    
    df$year <- df$year - firstharvest # Years from first harvest
    
    df%<>%
      group_by(scenario, year)%>%
      summarise(mean=mean(N/15), sd=sd(N/15))%>%
      mutate(se=(sd/sqrt(simul)))%>%
      mutate(hci=(mean)+(1.96*se), lci=mean-(1.96*se))
    

    p<-df%>%ggplot(aes(x=year, y=mean, fill=scenario))+
      scale_color_manual(values = rev(cols))+
      scale_fill_manual(values = rev(cols))+
      geom_line(alpha=0.7)+xlab("Years from first harvest")+
      ylab(expression("Colonies per m" ^2))+
      scale_y_continuous(limits = c(0,6)) +
      geom_ribbon(aes(ymin = lci, ymax = hci, fill=scenario),alpha = 0.7) +
      theme(legend.position = "none", panel.grid = element_blank(),panel.background = element_rect(fill="grey98"))+
      theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())+
      theme(axis.text = element_text(size = 9))+
      theme(axis.title = element_text(size = 9))+
      ggtitle(paste(recdat,"recruitment", sep=" "))
    
    
    # Add annotation only to one plot (i == 1)
    if (i == 2) {
      p <- p +
        geom_text(aes(y = 6 , x=8,label=scenarios[1]),size = 3, color = rev(cols)[1])+
        geom_text(aes(y = 6 , x=23,label=scenarios[2]),size = 3, color = rev(cols)[2])
    }
    
    assign(paste0("pop_", recdat), p, envir = plot_env)
  }
}

pop_plots(model_results, recdat) ## population density

##---------------------------------------------------
## Size-structure of populations 
##---------------------------------------------------

sfmean=data.frame(scenario=NA, classes=NA, smean=NA, hci=NA, lci=NA)

sizefreq_plots <- function(results, names) {
  
  for (i in seq_along(model_results)){
    df <- model_results[[i]][[2]]  # size-freq data
    simul=max(df$simul)
    recdat=unique(model_results[[i]][[1]]$recdat)
    
    df$scenario <- factor(df$scenario)
    df$scenario <- relevel(df$scenario, ref = "nonharvested")
    levels(df$scenario)[1]<-"unharvested"
    
    df$year[df$scenario=="harvested"] <- df$year[df$scenario=="harvested"]-1 # Years from first harvest
    
    df%<>% ## prepare data for plotting
      filter(year>firstharvest)
      
    
    df$classes=NA
    df$classes[df$class==1]<-"0-5"
    df$classes[df$class==2]<-"5-10"
    df$classes[df$class==3]<-"10-20"
    df$classes[df$class==4]<-"20-30"
    df$classes[df$class==5]<-"30-40"
    df$classes[df$class==6]<-"40-50"
    df$classes[df$class==7]<-">50"
    
    df$classes=factor(df$classes, levels=c("0-5","5-10","10-20","20-30","30-40","40-50",">50"))
    
    
    tot=df%>%group_by(scenario, simul, year)%>%mutate(n=sum(N))
    tot$relative<-tot$N/tot$n
    tot%<>%ungroup()%>%dplyr::select(scenario, classes,relative)
    
    p<-tot%>%filter(!classes=="0-5")%>%
      group_by(scenario,classes)%>%
      ggplot(aes(y = relative,
                   x = classes, 
                   fill = scenario)) +
      geom_bar(stat = "identity", position =position_dodge(width = 0.6))+
      xlab("Size class (cm)") +ylab("Relative abundance")+
      scale_color_manual(values = rev(cols))+
      scale_fill_manual(values = rev(cols))+
      scale_y_continuous(limits = c(0,0.5)) +
      theme(legend.position = "none")+
      theme(axis.text = element_text(size = 8))+
      theme(axis.title = element_text(size = 9))+
      ggtitle(paste(recdat,"recruitment", sep=" "))
      
      # Add annotation only to the first plot (i == 1)
      if (i == 1) {
        p<- p +
          geom_text(aes(y = max(relative), x=classes[4],label=scenarios[1]),size = 4, color = rev(cols)[1])+
          geom_text(aes(y = max(relative)-0.1, x=classes[4],label=scenarios[2]),size = 4, color = rev(cols)[2])
      }
    
    assign(paste0("proportional_sfd",recdat), p, envir = plot_env)
  }
}

sizefreq_plots(model_results, recdat)## Relative size-distribution no recruits

##----------------
##Recruit density
##----------------

recruitment_plots <- function(results, names) {
  for (i in seq_along(model_results)){
    df <- model_results[[i]][[3]]  # age-freq data
    recdat=unique(model_results[[i]][[1]]$recdat)
    
    df$scenario <- factor(df$scenario)
    df$scenario <- relevel(df$scenario, ref = "nonharvested")
    
    levels(df$scenario)[1]<-"unharvested"
    
    df$year[df$scenario=="harvested"] <- df$year[df$scenario=="harvested"]-1 # Years from first harvest
    
    df%<>%filter(year > firstharvest-2) ##remove burning period
    df$year <- df$year - firstharvest # Years from first harvest
    
    
    df1=df%>% ## prepare data for plotting
      filter(Age==1)%>%
      group_by(scenario, year)%>%
      summarise(mean=mean(N/15), sd=sd(N/15))%>%##across simulations
      mutate(se=(sd/sqrt(simul)))%>%
      mutate(hci=(mean)+(1.96*se), lci=mean-(1.96*se))
    
    
    p1<-df1%>%ggplot(aes(x=year, y=mean, fill=scenario))+
      scale_color_manual(values = rev(cols))+
      scale_fill_manual(values = rev(cols))+
      geom_line(alpha=0.7)+xlab("Years from first harvest")+
      ylab(expression("Recruits per m" ^2))+
      geom_ribbon(aes(ymin = lci, ymax = hci, fill=scenario),alpha = 0.7) +
      theme(legend.position = "none", panel.grid = element_blank(),panel.background = element_rect(fill="grey98"))+
      theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())+
      theme(axis.text = element_text(size = 9))+
      theme(axis.title = element_text(size = 9))

        assign(paste0("Recruitment_",recdat), p1, envir = plot_env)
  }
  }


recruitment_plots(model_results, recdat) ## Recruits


##--------------------------------------------------
## Shift in density and size of  Adult 
##--------------------------------------------------

adult_plots <- function(results, names) {
  
  for (i in seq_along(model_results)){
    df <- model_results[[i]][[4]]
    recdat=unique(model_results[[i]][[1]]$recdat)
    
    df$scenario <- factor(df$scenario)
    df$scenario <- relevel(df$scenario, ref = "nonharvested")
    levels(df$scenario)[1]<-"unharvested"
    
    df$year[df$scenario=="harvested"] <- df$year[df$scenario=="harvested"]-1 # Years from first harvest
    
    df1=df%>% filter(year >firstharvest-1)%>%
      filter(Age>1)%>%##Only adults
      group_by(scenario, simul)%>%
      summarise(mheight=mean(height), sd=sd(height), n=length(unique(year)))
    
    df1%<>%##across simulations
      mutate(se=(sd/sqrt(n)))%>%
      mutate(hci=(mheight)+(1.96*se), lci=mheight-(1.96*se))%>%ungroup()
    
    df2<-df1%>%dplyr::select(-simul)%>%group_by(scenario)%>%
      summarise(mean=mean(mheight), hci=mean(hci),lci=mean(lci), .groups = "drop")
    
    
    
    p1<-df1%>%ggplot(aes(x = mheight, fill = scenario, color=scenario)) +
      geom_density(alpha=0.7, bw=1)+
      xlim(min(df1$mheight)-5,max(df1$mheight)+5)+
      geom_segment(aes(x = 20, xend = 20, y = 0, yend = 0.4), 
                   color = "black", linetype = "dashed", linewidth = 1) +
      ylab("Probability Density")+xlab("Adults mean Height (cm)")+
      scale_color_manual(values = rev(cols))+
      scale_fill_manual(values = rev(cols))+
      theme(legend.position = "none",  
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      geom_text(aes(y = 0.35, x=25,label=scenarios[1]),size = 4, color = rev(cols)[1])+
      geom_text(aes(y = 0.35, x=15,label=scenarios[2]),size = 4, color = rev(cols)[2])+
      ggtitle(paste(recdat,"recruitment", sep=" "))
    
    
    p3<-p1+
      geom_point(df1,mapping=aes(x = mheight, y=-0.01,fill = scenario, color=scenario))+
      geom_point(df2, mapping=aes(x=mean, y=-0.01),color="black", size=3)
      
    assign(paste0("MeanAdulHeight_",recdat), p3, envir = plot_env)
    
  }
}

adult_plots(model_results, recdat) # Probability density of adult size

##---------------------------
##--Yield harvested
##---------------------------

cols2 <-  wes_palette(5, 15, type = c("continuous"))

yield_plots <- function(results, names) {
  
  for (i in seq_along(model_results)){
    df <- model_results[[i]][[5]]  ## Harvested Individuals
    recdat=unique(model_results[[i]][[1]]$recdat)
    
    df$year <- df$year - firstharvest # Years from first harvest
    
    p1<-df%>%ggplot(aes(x=factor(year), y=count, fill=factor(year)))+geom_boxplot()+
      ylab(expression("Yield (Ind) per 15m" ^2))+
      xlab("Years from first harvest")+
      ylim(0,16)+
      scale_fill_manual(values = cols2)+
      theme(axis.text = element_text(size = 8))+
      theme(axis.title = element_text(size = 8))+
      theme(legend.position = "none",  
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      ggtitle(paste(recdat,"recruitment", sep=" "))+
      theme(plot.title = element_text(hjust = 0.5, vjust = -5, size = 10, color = "black"))
        
        ## Track proportion of populations that are harvested each year
        df1=df
        
        df%<>%
        complete(simul,year,fill=list(count=0))%>%ungroup()
        df$harvest=1
        df$harvest[df$count==0]=0
        
        ##Proportion of harvested populations each year
        df%<>%group_by(year,harvest)%>%summarise(reps=n())%>%mutate(prop=reps/sum(reps))
        
        h<-max(df$year)
        m<-median(df$year)
        
        p2<-ggplot(df, aes(x = factor(year), y = prop, fill = factor(harvest))) +
          geom_bar(stat = "identity") +
          scale_fill_manual(values = rev(cols))+ 
          xlab("Years from first harvest") +ylab("Proportion of populations")+
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size = 8),
            axis.title = element_text(size = 8))+
          annotate("label",fill = rev(cols)[1],label=scenarios[1], x=6, y=0.95,size = 3, color = rev(cols)[2])+
          annotate("label",fill = rev(cols)[2],label=scenarios[2], x=4, y=0.1,size = 3, color = rev(cols)[1])
            
    assign(paste0("Yield_",recdat), p1, envir = plot_env)
    assign(paste0("Harvested_",recdat), p2, envir = plot_env)
    
    df2 <- model_results[[i]][[6]]  ## Yield= Harvested Biomass
    df2$year <- df2$year - firstharvest # Years from first harvest
    
    ##Effort/yield
    df2%<>%left_join(df1)
    names(df2)=c("yield", "simul","year","colonies")
    
    ##Complete cases
    df2%<>%
      complete(simul,year,fill=list(yield=0,colonies=0))%>%ungroup()
    df2$harvest=1
    df2$harvest[df2$yield==0]=0
    
    ##Averaged across harvested populations
    df.harvested<-df2%>%filter(harvest==1)%>%
      group_by(year)%>%
      summarise(CPUhA=mean(yield), harv.Colonies=mean(colonies))%>%ungroup()
    df.harvested$meanPercapita=df.harvested$CPUhA/df.harvested$harv.Colonies ##per harvested area (prop)
   
     df.visited<-df2%>%
      group_by(year)%>%
      summarise(CPUvA=mean(yield), meanColonies=mean(colonies))%>%ungroup()
    
df2=df.harvested%>%left_join(df.visited)
    
    p3<-df2%>%
      ggplot(aes(x=year))+
      geom_line(aes(y=CPUhA),colour=cols[1], size=1)+
      ylab(expression("Biomass (g) per 15m" ^2))+
      geom_line(aes(y=CPUvA),colour=cols[1], size=1, linetype="dotted", alpha=0.6)+
      xlab(expression("year from first harvest"))+
      theme(axis.text = element_text(size = 8))+
      theme(axis.title = element_text(size = 8))+
      theme(legend.position = "none",  
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      theme(plot.title = element_text(hjust = 0.5, vjust = -5, size = 10, color = "black"))
    
    # Add annotation only to the first plot (i == 1)
    if (i == 1) {
      ymax=max(df2$meanPercapita)
      va=mean(df2$CPUvA)*0.5
      ha=mean(df2$CPUhA)+1
      
      p3 <- p3 +
      #annotate("text",label="per Capita", x=8, y=ymax,size = 4, color = cols[1], alpha=0.6)+
      annotate("text",label="harvested", x=20, y=ha,size = 3, color = cols[1],fontface = "bold")+
      annotate("text",label="visited", x=10, y=va,size = 3, color = cols[1], alpha=0.6)
    }

    assign(paste0("Biomass_",recdat), p3, envir = plot_env)
    
    ##perCapita
    p4<-df2%>%
      ggplot(aes(x=year))+
      geom_line(aes(y=meanPercapita),colour=cols[1], size=1, alpha=0.6)+
      xlab(expression("year from first harvest"))+
      ylab("Biomass per Capita")+
      expand_limits(y=c(2,max(df2$meanPercapita)))+
      theme(axis.text = element_text(size = 8))+
      theme(axis.title = element_text(size = 8))+
      theme(legend.position = "none",  
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      theme(plot.title = element_text(hjust = 0.5, vjust = -5, size = 10, color = "black"))
    
    assign(paste0("ColonyBiomass_",recdat), p4, envir = plot_env)
    
  }
}

yield_plots(model_results,recdat)## Yield


##---------------------------
## Final Plots 
##---------------------------

## Fig2. Validation

##Load Observed SFD SS 

SFD<-readxl::read_excel("Parameters.xlsx",sheet="SFD_SS")
SFD%<>%group_by(classes)%>%filter(!classes=="0-5")%>%
  summarise(N=mean(N2))
##relative
SFD$n=sum(SFD$N)
SFD$relative<-SFD$N/SFD$n
SFD$scenario<-"observed"
SFD%<>%dplyr::select(scenario,classes,relative)

SFD$classes=factor(SFD$classes, levels=c("5-10","10-20","20-30","30-40","40-50",">50"))


obs<-SFD%>%filter(!classes=="0-5")%>%
  ggplot(aes(y = relative,
             x = classes)) +
  geom_bar(stat = "identity", position =position_dodge(width = 0.6))+
  xlab("Size class (cm)") +ylab("Relative abundance")+
  scale_y_continuous(limits = c(0,0.5)) +
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 9))+
  ggtitle("Observed")

FigSFD<-ggarrange(print(plot_env$proportional_sfdHigh),
                  print(plot_env$proportional_sfdLow)+ylab(""),
                  obs+ylab(""), ncol=3)



## Extract Yield Predictions

high_yield=plot_env$Yield_High$data%>%group_by(year)%>%
  summarise(yield=mean(count), min=min(count), max=max(count))%>%
  mutate(scenario=plot_env$Yield_High$labels$title)

low_yield=plot_env$Yield_Low$data%>%group_by(year)%>%
  summarise(yield=mean(count), min=min(count), max=max(count))%>%
  mutate(scenario=plot_env$Yield_Low$labels$title)


predyield=rbind(high_yield,low_yield)

##Observed Yield data
observed.harvest=data_frame(year=10, yield=6.23, sd=5.9, scenario="observed")

reccols=wes_palette(5, 2, type = c("discrete"))

y1<-predyield%>%group_by(scenario)%>%summarise(yield=mean(yield),min=mean(min), max=mean(max))%>%
  ggplot(aes(x=scenario, y=yield), color=reccols)+
  geom_pointrange(aes(ymin=min, ymax=max),color=reccols)+
  geom_point(observed.harvest, mapping=aes(y=yield, x=scenario),color="black", size=3)+
  geom_segment(observed.harvest, mapping=aes(y=yield-sd, yend=yield+sd),color="black")+
  xlab("")+
  ylab(expression("Harvested colonies per 15m"^2))+
  theme(axis.text = element_text(size = 9))+
  theme(axis.title = element_text(size = 9))+
  theme(legend.position = "none",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = -5, size = 13, color = "black"))

Validation<-ggarrange(FigSFD/y1, align="none")

#ggsave(paste(folder,"/","Fig2_validation.png", sep=""), Validation, width = 8.5, height = 5, units = "in", dpi=300)


## Fig3. Pop dynamics

Fig_pop<-ggarrange(print(plot_env$pop_High),print(plot_env$pop_Low),
                   print(plot_env$Recruitment_High),print(plot_env$Recruitment_Low))

#ggsave(paste(folder,"/",Sys.Date(),"Fig3_popdynamics.png", sep=""), Fig_pop, width = 5, height = 4.5, dpi=400)      


## Fig4. Size probability distribution 

Fig_Asize<-ggarrange(print(plot_env$MeanAdulHeight_High),print(plot_env$MeanAdulHeight_Low)+ylab(""))

#ggsave(paste(folder,"/",Sys.Date(),"Fig4_adultsize.png", sep=""), Fig_Asize, width = 5, height = 4, units = "in", dpi = 300)


## Fig5. Yield


Yield<-ggarrange(print(plot_env$Yield_High),print(plot_env$Yield_Low)+ylab(""),
                 print(plot_env$Harvested_High),print(plot_env$Harvested_Low)+ylab(""),
                 print(plot_env$Biomass_High)+ylim(0,40),print(plot_env$Biomass_Low)+ylab("")+ylim(0,40), 
                 print(plot_env$ColonyBiomass_High)+ylim(1,6),print(plot_env$ColonyBiomass_Low)+ylab("")+ylim(1,6),
                 nrow = 4, ncol = 2, labels=c("a","","b","","c","","d"))


#ggsave(paste(folder,"/",Sys.Date(),"Fig5_Yield.png", sep=""), Yield, width = 5, height = 8, units = "in", dpi=400)


#save(plot_env, file = paste(folder,"/",Sys.Date(),"plots.RData", sep=""))
