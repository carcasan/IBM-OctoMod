These scripts built the Individual-based model of Antillogorgia elisabethae for modelling population dynamics and harvest impacts in the Bahamas


1.	1_VBCurve_estimates.R fits the growth curve based on empirical age-height relationship and define the parameters to predict height at age and growth rate in the model.
Outputs saved in "VB_CHparams.RData"

3.	2_Initial_Population.R uses the parameters estimated above ("VB_CHparams.RData") to create the initial age, density and size-frequency conditions for each simulated population.
Outputs saved in “Initialization.RData" 

5.	Set_parameters.R  adds partial and total mortality probabilities from field data and gathers all parameters (i.e. from steps 1 & 2).
Outputs saved in "Set_parameters.RData" which is called for implementing the model

7.	Create_metapopulations.R  describe the functions to create the neighbouring populations and estimate the number of recruits arriving to the focal population following density-dependent recruitment as a function of colony size of reproductive adults. Functions saved in "Metapopulations.Rdata"

8.	Recruitment.R describes the function to add the total number of recruits that successfully settled (local + external) to the population

9.	Run_Main.R Sets the default conditions to run the model, such as number of simulations, years, scenarios and harvest conditions. 

10.	IBModel.R The IBM model function (AEmod).  The function and all parameters needed to run the default model are saved in "model_data.Rdata"

11.	Mod_Paralell.R Runs the model recruitment scenarios in paralell

12.	Plot_Outputs.R summarises model outputs and creates figures for paper

13.	Set_sensitivity.R Selects the parameter of interest to run sensitivity

14.	Run_Sensitivity.R  Runs model for sensitivity

15.	Plot_Sensitivity.R  Gathers all outputs and creates figure
