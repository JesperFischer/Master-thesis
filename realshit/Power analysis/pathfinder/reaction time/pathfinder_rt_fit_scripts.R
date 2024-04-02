
power_analysis_without_psi = function(parameters){
  
  
  source(here::here("realshit","Power analysis","pathfinder","reaction time", "pathfinder_rt_datasets_scripts.R"))
  
  source(here::here("realshit","Power analysis","pathfinder","reaction time", "pathfinder_rt_fit_scripts.R"))
  
  
  pattern = paste0("effect_size_alpha = ",
                   as.character(parameters$effect_size_alpha[1]),
                   " effect_size_beta = ",
                   as.character(parameters$effect_size_beta[1]),
                   " trials = ",as.character(parameters$trials[1]),
                   " subjects = ",as.character(parameters$subjects[1]))
  
  
  # Get a list of all files in the directory
  files <- sort(list.files(path = here::here("realshit","Power analysis","pathfinder","reaction time","datasets_single",pattern), full.names = TRUE))
  
  data = read.csv(files[parameters$id])
  
  
  #data = data %>% filter(participant_id < 20)
  
  #wrangling the data
  
  fit_model = function(data){
    
    mod_noncent = cmdstanr::cmdstan_model(here::here("realshit","Power analysis","pathfinder","reaction time","parameterized_NORT.stan"),stanc_options = list("O1"))
    mod_noncent = cmdstanr::cmdstan_model(here::here("realshit","Power analysis","pathfinder","reaction time","parameterized_NORT_halfparam.stan"),stanc_options = list("O1"))
    
    data_trans = data %>% 
      mutate(session = ifelse(sessions == 1, 0 ,ifelse(sessions == 2, 1, 2)))
    
    data_trans = transform_data_to_stan(data_trans)
    
    #data = data %>% filter(sessions == 1 |sessions == 3)
    
    data_trans = data_trans %>% arrange(sessions, participant_id)
    
    data_stan = list(T = nrow(data_trans),
                     S = length(unique(data_trans$participant_id)),
                     S_id = as.numeric(data_trans$participant_id ),
                     X = data_trans %>% .$X,
                     X_lapse = as.matrix(data.frame(int = rep(1,nrow(data_trans)))),
                     X_alpha = as.matrix(data.frame(int = rep(1,nrow(data_trans)),
                                                    session = data_trans %>% .$sessions)),
                     X_beta = as.matrix(data.frame(int = rep(1,nrow(data_trans)),
                                                   session = data_trans %>% .$sessions)),
                     N_alpha = 2,
                     N_beta = 2,
                     N_lapse = 1,
                     Y = data_trans %>% .$resp,
                     npx = data_trans %>% .$npx
    )
    
    
    fit_norm <- mod_noncent$sample(
      data = data_stan,
      iter_sampling = 1000,
      iter_warmup = 1000,
      chains = 4,
      init = 0,
      parallel_chains = 4,
      refresh = 500,
      adapt_delta = 0.8,
      max_treedepth = 10
    )
    
    diags_nocentered = data.frame(fit_norm$diagnostic_summary())
    rhat_nocentered = data.frame(fit_norm$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                                    "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% 
      .$rhat
    
    
    
    data_rows = transform_data_to_stan(data)
    data_rows = data_rows %>% arrange(sessions, participant_id)
    
    diags = data.frame(fit_norm$diagnostic_summary())
    
    
    alpha_draws = as_draws_df(fit_norm$draws("gm[2]")) %>% .$`gm[2]`
    beta_draws = as_draws_df(fit_norm$draws("gm[4]")) %>% .$`gm[4]`
    
    p_zero_alpha = sum(alpha_draws<0)/length(alpha_draws)
    
    
    
    if(data$parameter[1] == "both"){
      if(data$real_effectsize_beta[1] > 0){
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }else if(data$real_effectsize_beta[1] < 0){
        p_zero_beta = sum(beta_draws>0)/length(beta_draws)
      }else{
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }
    }
    if(data$parameter[1] == "beta"){
      if(data$real_effectsize_beta[1] > 0){
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }else if(data$real_effectsize_beta[1] < 0){
        p_zero_beta = sum(beta_draws>0)/length(beta_draws)
      }else{
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }
    }
    if(data$parameter[1] == "alpha"){
      p_zero_beta = sum(beta_draws<0)/length(beta_draws)
    }
    
    
    rows = data_rows %>% dplyr::group_by(sessions,participant_id) %>% 
      dplyr::summarize(n = n()) %>% dplyr::ungroup() %>% dplyr::mutate(rows = cumsum(n)) %>% .$rows
    
    
    alphas = data.frame(fit_norm$summary("alpha"))[rows,]
    
    alphas = alphas %>% dplyr::mutate(real_values = data_rows$alpha[rows],
                                      parameter = "alpha",
                                      participant_id = data_rows$participant_id[rows]
    )
    
    betas = data.frame(fit_norm$summary("beta"))[rows,]
    
    betas = betas %>% mutate(real_values = data_rows$beta[rows],
                             parameter = "beta",
                             participant_id = data_rows$participant_id[rows]
    )
    
    
    
    
    #saving
    difs = data.frame(fit_norm$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                         "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% mutate(iter = data$iter[1])%>% 
      mutate(obs_effectsize_alpha = data$observed_effectsize_alpha[1],
             obs_effectsize_beta = data$observed_effectsize_beta[1],
             subjects = data$subjects[1],
             trials = data$trials[1],
             sim_alphacor = data$sim_alphacor[1],
             sim_betacor = data$sim_betacor[1],
             divergences = mean(diags$num_divergent),
             treedepths = mean(diags$num_max_treedepth),
             p_alpha = p_zero_alpha,
             p_beta = p_zero_beta,
             real_effectsize_alpha = data$real_effectsize_alpha[1],
             real_effectsize_beta = data$real_effectsize_beta[1],
             centered_model = data.frame(as_draws_df(fit_norm$draws("centered")))[1,1])
    
    
    
    
    
    
    return(list(difs))
    
  }
  
  norts = fit_model(data)
  
  
  return(norts)
  
  
  fit_model_rt = function(data){
    
    mod_noncent_rt = cmdstanr::cmdstan_model(here::here("realshit","Power analysis","pathfinder","reaction time","parameterized_RT.stan"),stanc_options = list("O1"))
    
    data_trans = data %>% 
      mutate(session = ifelse(sessions == 1, 0 ,ifelse(sessions == 2, 1, 2)))
    
    data_trans = transform_data_to_stan(data_trans)
    
    #data = data %>% filter(sessions == 1 |sessions == 3)
    
    data_trans = data_trans %>% arrange(sessions, participant_id)
    
    data_stan = list(T = nrow(data_trans),
                     S = length(unique(data_trans$participant_id)),
                     S_id = as.numeric(data_trans$participant_id),
                     X = data_trans %>% .$X,
                     X_alpha = as.matrix(data.frame(int = rep(1,nrow(data_trans)),
                                                    session = data %>% .$sessions)),
                     X_beta = as.matrix(data.frame(int = rep(1,nrow(data_trans)),
                                                   session = data_trans %>% .$sessions)),
                     
                     Y = data_trans %>% .$resp,
                     RT = data_trans %>% .$rts,
                     min_RT = data_trans %>% group_by(participant_id) %>% dplyr::summarize(min = min(rts)) %>% .$min
    )
    
    
    fit_norm <- mod_noncent_rt$sample(
      data = data_stan,
      iter_sampling = 1000,
      iter_warmup = 1000,
      chains = 4,
      parallel_chains = 4,
      refresh = 500,
      init = 0,
      adapt_delta = 0.8,
      max_treedepth = 10
    )
    
    diags_nocentered = data.frame(fit_norm$diagnostic_summary())
    rhat_nocentered = data.frame(fit_norm$summary(c(paste0("gm[",1:9,"]"),paste0("tau_u[",1:9,"]")))) %>% 
      .$rhat
    
    
    data_rows = transform_data_to_stan(data)
    data_rows = data_rows %>% arrange(sessions, participant_id)
    
    diags = data.frame(fit_norm$diagnostic_summary())
    
    
    alpha_draws = as_draws_df(fit_norm$draws("gm[2]")) %>% .$`gm[2]`
    beta_draws = as_draws_df(fit_norm$draws("gm[4]")) %>% .$`gm[4]`
    
    p_zero_alpha = sum(alpha_draws<0)/length(alpha_draws)
    
    
    
    if(data$parameter[1] == "both"){
      if(data$real_effectsize_beta[1] > 0){
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }else if(data$real_effectsize_beta[1] < 0){
        p_zero_beta = sum(beta_draws>0)/length(beta_draws)
      }else{
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }
    }
    if(data$parameter[1] == "beta"){
      if(data$real_effectsize_beta[1] > 0){
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }else if(data$real_effectsize_beta[1] < 0){
        p_zero_beta = sum(beta_draws>0)/length(beta_draws)
      }else{
        p_zero_beta = sum(beta_draws<0)/length(beta_draws)
      }
    }
    if(data$parameter[1] == "alpha"){
      p_zero_beta = sum(beta_draws<0)/length(beta_draws)
    }
    
    
    rows = data_rows %>% dplyr::group_by(sessions,participant_id) %>% 
      dplyr::summarize(n = n()) %>% dplyr::ungroup() %>% dplyr::mutate(rows = cumsum(n)) %>% .$rows
    
    
    
    #saving
    difs = data.frame(fit_norm$summary(c(paste0("gm[",1:9,"]"),paste0("tau_u[",1:9,"]")))) %>% 
      mutate(iter = data$iter[1])%>% 
      mutate(obs_effectsize_alpha = data$observed_effectsize_alpha[1],
             obs_effectsize_beta = data$observed_effectsize_beta[1],
             subjects = data$subjects[1],
             trials = data$trials[1],
             sim_alphacor = data$sim_alphacor[1],
             sim_betacor = data$sim_betacor[1],
             divergences = mean(diags$num_divergent),
             treedepths = mean(diags$num_max_treedepth),
             p_alpha = p_zero_alpha,
             p_beta = p_zero_beta,
             real_effectsize_alpha = data$real_effectsize_alpha[1],
             real_effectsize_beta = data$real_effectsize_beta[1],
             centered_model = data.frame(as_draws_df(fit_norm$draws("centered")))[1,1])
    
    
    
    
    
    
    return(list(difs))
  }
  
  rts = fit_model_rt(data)
  
  return(list(norts,rts))
  
  
}


transform_data_to_stan = function(data){
  
  names(data)
  
  data = data %>% group_by(X,participant_id, sessions,alpha,beta)%>%
    summarise(X=mean(X),
              npx=n(),
              resp=resp,
              rts = rts)%>%
    ungroup()
  
  return(data)
}
