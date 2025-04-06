
##----------------
## RUN scenarios
##----------------

##List of scenarios to run
Simulation_runs <- read.csv("Simulation_runs.csv")

Model_run=Simulation_runs[1,]##Default model


##Load parameters 
for (col_name in colnames(Model_run)) {
  value <- Model_run[[col_name]]
  assign(col_name, value, envir = .GlobalEnv)
}

Simul=as.numeric(simul)
scenarios=strsplit(scenarios, "_")[[1]]
Recdata=strsplit(Recdata, "_")[[1]]

lastharvest=firstharvest+totalharvest*3
harvestyears=seq(firstharvest,lastharvest,frequency)
years=c(1:(lastharvest+2))

##set residual heights
rangeobs1=c(2,14)
rs.height=RGeode::rgammatr(length(years), 0.15, 1/40,range = rangeobs1)

