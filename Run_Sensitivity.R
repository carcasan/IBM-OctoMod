##------------------
### RUN SENSITIVITY
##------------------

source("IBModel.R")
simul_runs <- ("Simulation_runs.csv")
Model_run<-read.csv(simul_runs)
#Model_run<-Model_run[-1,]# not main model

# ##-----------------------------------------------------
# # Run models in parallel for each recruitment condition
# ##-----------------------------------------------------

library(parallel)

model<-paste(here(),"model_data.Rdata", sep="/")
params<-paste(here(),"sensitivity_scenario.Rdata", sep="/")

# Define the function to run the simulation in parallel
run_sensitivity <- function() {
  # Set up parallel workers (use 4 cores, adjust as necessary)
  num_cores <- 4
  cl <- makeCluster(num_cores)
  # Export the .RData file to the workers
  model<-paste(here(),"model_data.Rdata", sep="/")
  simul_runs <- ("~/Documents/Bahamas_OctoMod/Paper/Scripts&Data/Simulation_runs.csv")
  
  clusterExport(cl, c("model","scn","Model_run","simul_runs","params"))  # Make sure the file is available to workers
  
  # Load the .RData file on each worker
  clusterEvalQ(cl, {
    load(model)  # Load the .RData file into each worker's environment
    
    Model_run<-read.csv(simul_runs)
    
    # Set the CRAN mirror (important for package installation)
    options(repos = c(CRAN = "https://cran.rstudio.com/"))  # Set a CRAN mirror
    
    # Load the required packages on each worker
    for (pkg in required_packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cran.rstudio.com/")  # Set mirror here
        library(pkg, character.only = TRUE)  # Load the package
      }
    }
    source(paste(here(),"Run_sensitivity.R", sep="/"))

  })
  
  # Use parLapply to run the loop in parallel for each level in Rec
  results <- parLapply(cl, 1:length(Recdata), function(r) {
    recdat <- Recdata[r]  # Get the specific level from Rec
    model <- AEmod(simul, recdat)  # Run the model for the current level and simulation
    return(model)  # Return the result for this iteration
  })
  
  return(results)
  
  # Stop the cluster
  stopCluster(cl)
  
}

#-----------------------------
# Run the parallel simulations
#-----------------------------
Model_run$Run

scn=2#c(2:length(Model_run)) Choose scenario option

source("Set_sensitivity.R")

model_results <- run_parallel_simulations()




