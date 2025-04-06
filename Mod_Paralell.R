source("IBModel.R") ## The model and associated functions
load("mypackages.Rdata") ## required packages
source("Run_Main.R") ## parameters for default model

##---------------------------------  
# # Parallel setup for simulations
##---------------------------------
library(parallel)

model<-paste(here(),"model_data.Rdata", sep="/") ## the path for output from "IBModel.R"

# Function to run simulations in parallel
run_parallel_simulations <- function() {
  # Set up parallel workers
  num_cores <- 4
  cl <- makeCluster(num_cores)
  # Export the .RData file to the workers
  model<-paste(here(),"model_data.Rdata", sep="/")
  params<-paste(here(),"Run_Main.R", sep="/")
  
  clusterExport(cl, c("model"))  # Make sure the file is available to workers
  
  # Load the .RData file on each worker
  clusterEvalQ(cl, {
    load(model)
    # Set the CRAN mirror (important for package installation)
    options(repos = c(CRAN = "https://cran.rstudio.com/"))  # Set a CRAN mirror
    
    # Load the required packages on each worker
    for (pkg in required_packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cran.rstudio.com/")  # Set mirror here
        library(pkg, character.only = TRUE)  # Load the package
      }
    }
    source(paste(here(),"Run_Main.R", sep="/"))
  })
  
  # Run each recruitment scenario
  results <- parLapply(cl, 1:length(Recdata), function(r) {
    recdat <- Recdata[r]  # Get the specific level from Rec
    model <- AEmod(simul, recdat)  # Run model for the sencenario and simulation
    return(model)
  })
  
  return(results)
  
  # Stop the cluster
  stopCluster(cl)
  
}

# Run the parallel simulations
model_results <- run_parallel_simulations()

##Save outputs
#save(model_results, Simulation_runs, file=paste(Run,"_ModelOutputs",Sys.Date(),".Rdata"))

 