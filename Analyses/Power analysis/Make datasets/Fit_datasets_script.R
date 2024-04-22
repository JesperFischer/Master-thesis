
power_analysis_without_psi = function(parameters){
  
  
  source(here::here("realshit","Power analysis", "Make_datasets_scripts.R"))
  
  source(here::here("realshit","Power analysis", "Fit_datasets_script.R"))
  
  
  pattern = paste0("effect_size_alpha = ",
                   as.character(parameters$effect_size_alpha[1]),
                   " effect_size_beta = ",
                   as.character(parameters$effect_size_beta[1]),
                   " trials = ",as.character(parameters$trials[1]),
                   " subjects = ",as.character(parameters$subjects[1]))
  
  
  # Get a list of all files in the directory
  files <- sort(list.files(path = here::here("realshit","Power analysis","datasets",pattern), full.names = TRUE))
  
  data = read.csv(files[parameters$id])
  
  #wrangling the data
  
  fit_model = function(data){
    
    
    mod_noncent = cmdstanr::cmdstan_model(here::here("Stanmodels","Parameterized cummulative normal.stan"),stanc_options = list("O1"))
    
    data = data %>% 
      mutate(session = ifelse(sessions == 1, 0 ,ifelse(sessions == 2, 1, 2)))
    
    data = transform_data_to_stan(data)
    
    #data = data %>% filter(sessions == 1 |sessions == 3)
    
    data = data %>% arrange(sessions, participant_id)
    
    data_stan = list(T = nrow(data),
                     S = length(unique(data$participant_id)),
                     S_id = as.numeric(data$participant_id ),
                     X = data %>% .$X,
                     X_lapse = as.matrix(data.frame(int = rep(1,nrow(data)))),
                     X_alpha = as.matrix(data.frame(int = rep(1,nrow(data)),
                                                    session = data %>% .$sessions)),
                     X_beta = as.matrix(data.frame(int = rep(1,nrow(data)),
                                                   session = data %>% .$sessions)),
                     N_alpha = 2,
                     N_beta = 2,
                     N_lapse = 1,
                     Y = data %>% .$resp,
                     npx = data %>% .$npx
    )
    
    
    fit_nocentered <- mod_noncent$sample(
      data = data_stan,
      iter_sampling = 1000,
      iter_warmup = 1000,
      chains = 4,
      init = 0,
      parallel_chains = 4,
      refresh = 500,
      adapt_delta = 0.9,
      max_treedepth = 12
    )
    
    diags_nocentered = data.frame(fit_nocentered$diagnostic_summary())
    rhat_nocentered = data.frame(fit_nocentered$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                                          "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% 
      .$rhat
    
    return(fit_nocentered)
  }
  
  
  fit_norm = fit_model(data)
  
  
  data_rows = transform_data_to_stan(data)
  data_rows = data_rows %>% arrange(sessions, participant_id)
  
  diags = data.frame(fit_norm$diagnostic_summary())
  
  
  alpha_draws = as_draws_df(fit_norm$draws("gm[5]")) %>% .$`gm[5]`
  beta_draws = as_draws_df(fit_norm$draws("gm[2]")) %>% .$`gm[2]`
  
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
  
  
  rows = data_rows %>% group_by(sessions,participant_id) %>% 
    summarize(n = n()) %>% ungroup() %>% mutate(rows = cumsum(n)) %>% .$rows
  
  
  alphas = data.frame(fit_norm$summary("alpha"))[rows,]
  
  alphas = alphas %>% mutate(real_values = data_rows$alpha[rows],
                             parameter = "alpha",
                             participant_id = data_rows$participant_id[rows],
                             psiestimate = data %>% group_by(participant_id,sessions) %>% summarize(psi_estimate_threshold = last(Estimatedthreshold)) %>% .$psi_estimate_threshold
  )
  
  betas = data.frame(fit_norm$summary("beta"))[rows,]
  betas = betas %>% mutate(real_values = data_rows$beta[rows],
                           parameter = "beta",
                           participant_id = data_rows$participant_id[rows],
                           psiestimate = data %>% group_by(participant_id,sessions) %>% summarize(psi_estimate_slope = last(Estimatedslope)) %>% .$psi_estimate_slope)
  
  
  indi_estimates = rbind(alphas,betas)%>% mutate(iter = data$iter[1],
                                                 real_effectsize_alpha = data$real_effectsize_alpha[1],
                                                 real_effectsize_beta = data$real_effectsize_beta[1],
                                                 obs_effectsize_alpha = data$observed_effectsize_alpha[1],
                                                 obs_effectsize_beta = data$observed_effectsize_beta[1],
                                                 subjects = data$subjects[1],
                                                 trials = data$trials[1],
                                                 sim_alphacor = data$sim_alphacor[1],
                                                 sim_betacor = data$sim_betacor[1],
                                                 divergences = mean(diags$num_divergent),
                                                 treedepths = mean(diags$num_max_treedepth)
  )
  
  indi_estimates$session = session = rep(c(rep(1,length(unique(indi_estimates$participant_id))),rep(2,length(unique(indi_estimates$participant_id)))),2)
  

  #Plot!
  # indi_estimates %>% 
  #   mutate(session = rep(c(rep(1,length(unique(indi_estimates$participant_id))),rep(2,length(unique(indi_estimates$participant_id)))),2)) %>% 
  #   ggplot(aes(x = mean, y = real_values, xmin = q5, xmax = q95, col = session))+
  #   geom_pointrange()+facet_wrap(~parameter, scales = "free")
  
  
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
  
  
  
  return(list(difs,indi_estimates))
  
}


transform_data_to_stan = function(data){
  
  data = data %>% group_by(X,participant_id, sessions,alpha,beta)%>%
    summarise(X=mean(X),
              npx=n(),
              resp=sum(resp))%>%
    ungroup()
  
  return(data)
}
