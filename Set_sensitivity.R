##------------------------------
## Set parameter for  sensitivity
##------------------------------

##List of scenarios to run
Model_run <- read.csv("Simulation_runs.csv")

#scn=2# First for sensitivity
##Loop through scenarios
this_run=Model_run[scn,]##remove Default model

scn.name=this_run$Run

##Load the parameters used in the simulation
for (col_name in colnames(this_run)) {
  value <- this_run[[col_name]]
  assign(col_name, value, envir = .GlobalEnv)
}

Simul=as.numeric(simul)
scenarios=strsplit(scenarios, "_")[[1]]
Recdata=strsplit(Recdata, "_")[[1]]

lastharvest=firstharvest+totalharvest*3
harvestyears=seq(firstharvest,lastharvest,frequency)
years=c(1:(lastharvest+2))
rangeobs1=c(2,14)
rs.height=RGeode::rgammatr(length(years), 0.15, 1/40,range = rangeobs1)

#save.image(file="sensitivity_scenario.Rdata")
