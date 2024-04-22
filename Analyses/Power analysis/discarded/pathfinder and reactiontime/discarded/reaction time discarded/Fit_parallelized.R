pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)

Run_poweranalysis = function(subjects, trials, effectsize_alpha, effectsize_beta){
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)
  message("Number of CPU cores in R: ", parallelly::availableCores())
  
  
  source(here::here("realshit","Power analysis","pathfinder","reaction time", "pathfinder_rt_datasets_scripts.R"))
  source(here::here("realshit","Power analysis","pathfinder","reaction time", "pathfinder_rt_fit_scripts.R"))
  
  
  
  # 
  subjects = c(5,6,7,8,9,10,12,15,20)
  trials = c(20,30,40,50,75,100,150)
  
#  subjects = c(5,6)
 # trials = c(20,30)
  
  effectsize_alpha = seq(0,1,by = 0.2)
  effectsize_beta = 0
  
  subjects = subjects
  trials = trials
  effect_size_alpha = effectsize_alpha
  effect_size_beta = effectsize_beta
  
  parameter = "both"
  replicate = 1:20
  
  parameters = expand.grid(subjects = subjects,
                           trials = trials,
                           parameter = parameter,
                           psi_only = F,
                           effect_size_alpha = effect_size_alpha,
                           effect_size_beta = effect_size_beta,
                           replicate = replicate) %>% 
    mutate(id = 1:nrow(.))
  
  data_list <- split(parameters, parameters$id)
  
  
  #cores = 16*2
  
  plan(multisession, workers = cores)
  
  possfit_model = possibly(.f = power_analysis_without_psi, otherwise = "Error")
  
  
  results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE),.progress = TRUE)
  
  
  saveRDS(results,here::here("realshit","Power analysis","pathfinder","reaction time", "results","effectsize_alpha = 0.4, effectsize_beta = 0.2.rds"))
  
  
  
  directory = here::here("realshit","Power analysis","pathfinder","results",paste0(parameter),paste0("results_","subjects=",subjects,"_trials=",trials))
  
  if(!dir.exists(directory)){
    dir.create(directory)
  }
  saveRDS(results,here::here(directory,paste0("results_","subjects=",subjects,"_trials=",trials,"_effectsizealpha=",effectsize_alpha,"_effectsizebeta=",effectsize_beta,".rds")))
  
  print("saved it! ")
}



if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  # If arguments are provided, use them
  subjects <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
  trials <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
  effectsize <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
} else {
  # If no arguments are provided, use default values
  subjects <- 0  # Set your default value for subjects
  trials <- 0   # Set your default value for trials
  effectsize <- 0  # Set your default value for effect size
}

print("runnning")
print("alpha")

print(subjects)
print(trials)
print(effectsize)

# Call your function with the extracted arguments
Run_poweranalysis(subjects, trials, effectsize_alpha,effectsize_beta)


