pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)

Run_poweranalysis = function(subjects, trials, effectsize_alpha, effectsize_beta){
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)
  message("Number of CPU cores in R: ", parallelly::availableCores())
  
  
  source(here::here("realshit","Power analysis","pathfinder", "pathfinder_datasets_scripts.R"))
  source(here::here("realshit","Power analysis","pathfinder", "pathfinder_fit_scripts.R"))
  
  
  # subjects = 10
  # trials = 100
  # effectsize_alpha = 0.5
  # effectsize_beta = 1
  
  subjects = subjects
  trials = trials
  effect_size_alpha = effectsize_alpha
  effect_size_beta = effectsize_beta
  
  parameter = "both"
  replicate = 1:100
  
  parameters = expand.grid(subjects = subjects,
                           trials = trials,
                           parameter = parameter,
                           psi_only = F,
                           effect_size_alpha = effect_size_alpha,
                           effect_size_beta = effect_size_beta,
                           replicate = replicate) %>% 
    mutate(id = 1:nrow(.))
  
  data_list <- split(parameters, parameters$id)
  
  
  
  results = list()
  i = 1
  for(i in 1:length(data_list)){
    print(i)
    print(paste0("PSI_ONLY = ",parameters$psi_only[1]))
    print(paste0("FITTING = ", nrow(parameters)))
    
    tryCatch({
      testing <- power_analysis_without_psi(data_list[[i]])
      results[[i]] <- testing
    }, error = function(e) {
      cat("Error in iteration", i, ":", conditionMessage(e), "\n")
      results[[i]] <- "Error"
    })
    
  }
  
  
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


