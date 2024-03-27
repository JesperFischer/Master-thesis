




fit_pathfinder_static_logisticdata= function(parameters){
  
  N = parameters$trials
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  parm_ev = data.frame()
  
  
  x[1] = 0
  p[1] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (1/(1+exp(-(1.7/(parameters$real_beta_con))*(x[1]-parameters$real_alpha_con))))
  r[1] = rbinom(1,1,p[1])
  
  
  mod = cmdstanr::cmdstan_model(here::here("Stanmodels","pathfinder_logistic.stan"))
  
  for(i in 1:N){
    
    
    data_stan = list(trials = nrow((data.frame(r) %>% drop_na())),
                     x = data.frame(x) %>% drop_na() %>% .$x,
                     r = data.frame(r) %>% drop_na() %>% .$r)
    
    fit_psy <- mod$pathfinder(data = data_stan,refresh=0, init = 0)
    
    if(i < 5){
      x[i] = data.frame(fit_psy$summary("psy_alpha")) %>% .$mean
    }else{
      sign = sample(c(-1,1),1)
      x[i] = fit_psy$draws("psy_alpha")[sample(1:1000,1)]+sign*fit_psy$draws("psy_beta")[sample(1:1000,1)]
    }
    
    if(abs(x[i]) > 50){
      print("to high")
      x[i] = data.frame(fit_psy$summary("psy_alpha")) %>% .$mean
    }
    
    
    c_trial = data.frame(fit_psy$summary(c("psy_alpha","psy_beta","psy_lambda"))) %>% mutate(trials = i)
    
    parm_ev = rbind(c_trial,parm_ev)
    p[i] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (1/(1+exp(-(1.7/(parameters$real_beta_con))*(x[i]-parameters$real_alpha_con))))
    r[i] = rbinom(1,1,p[i])
    
  }
  
  df_trial = data.frame(X = x,
                        prob = p,
                        resp = r,
                        trials = 1:N) %>% mutate(trials = parameters$trials,
                                                 real_lambda_con = parameters$real_lambda_con,
                                                 real_alpha_con = parameters$real_alpha_con,
                                                 real_beta_con = parameters$real_beta_con,
                                                 real_lambda_uncon = parameters$real_lambda_uncon,
                                                 real_alpha_uncon = parameters$real_alpha_uncon,
                                                 real_beta_uncon = parameters$real_beta_uncon,
                                                 id =  rnorm(1,0,1000),
                                                 subs = parameters$subs,
                                                 sessions = parameters$sessions,
                                                 model = "normal")
  
  #df_trial = full_join(parm_ev,df_trial)
  
  
  return(list(df_trial))
}



fit_pathfinder_static_normaldata = function(parameters){
  
  N = parameters$trials
  
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  parm_ev = data.frame()
  
  x[1] = 0
  p[1] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (0.5+0.5*erf((x[1]-parameters$real_alpha_con)/(parameters$real_beta_con*sqrt(2))))
  r[1] = rbinom(1,1,p[1])
  
  
  mod = cmdstanr::cmdstan_model(here::here("Stanmodels","pathfinder.stan"))
  
  for(i in 1:N){
    
    
    data_stan = list(trials = nrow((data.frame(r) %>% drop_na())),
                     x = data.frame(x) %>% drop_na() %>% .$x,
                     r = data.frame(r) %>% drop_na() %>% .$r)
    
    fit_psy <- mod$pathfinder(data = data_stan,refresh=0, init = 0)
    
    if(i < 5){
      x[i] = data.frame(fit_psy$summary("psy_alpha")) %>% .$mean
    }else{
      sign = sample(c(-1,1),1)
      x[i] = fit_psy$draws("psy_alpha")[sample(1:1000,1)]+sign*fit_psy$draws("psy_beta")[sample(1:1000,1)]
    }
    
    if(abs(x[i]) > 50){
      print("to high")
      x[i] = data.frame(fit_psy$summary("psy_alpha")) %>% .$mean
    }
    
    
    c_trial = data.frame(fit_psy$summary(c("psy_alpha","psy_beta","psy_lambda"))) %>% mutate(trials = i)
    
    parm_ev = rbind(c_trial,parm_ev)
    p[i] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (0.5+0.5*erf((x[i]-parameters$real_alpha_con)/(parameters$real_beta_con*sqrt(2))))
    r[i] = rbinom(1,1,p[i])
    
  }
  
  df_trial = data.frame(X = x,
                        prob = p,
                        resp = r,
                        trials = 1:N) %>% mutate(trials = parameters$trials,
                                                 real_lambda_con = parameters$real_lambda_con,
                                                 real_alpha_con = parameters$real_alpha_con,
                                                 real_beta_con = parameters$real_beta_con,
                                                 real_lambda_uncon = parameters$real_lambda_uncon,
                                                 real_alpha_uncon = parameters$real_alpha_uncon,
                                                 real_beta_uncon = parameters$real_beta_uncon,
                                                 id =  rnorm(1,0,1000),
                                                 subs = parameters$subs,
                                                 sessions = parameters$sessions,
                                                 model = "normal")
  
  #df_trial = full_join(parm_ev,df_trial)
  
  
  return(list(df_trial))
}


get_params = function(subs, trials, beta, lambda, alpha){
  
  alpha_uncon = rnorm(subs,alpha,10)
  beta_uncon = rnorm(subs,beta,0.6)
  lambda_uncon = rnorm(subs,lambda,2)    

  
  
  parameters = data.frame(real_alpha_uncon = alpha_uncon,
                          real_lambda_uncon = lambda_uncon,
                          real_beta_uncon = beta_uncon,
                          trials = trials,
                          subs = 1:subs) %>% 
    mutate(id = 1:n())
  

  parameters$real_beta_con = exp(parameters$real_beta_uncon)
  parameters$real_lambda_con = brms::inv_logit_scaled(parameters$real_lambda_uncon) / 2
  parameters$real_alpha_con = parameters$real_alpha_uncon
  
  data_list2 <- split(parameters, parameters$id)
  
  
  plan(multisession, workers = 4)
  
  #adding safety for if something goes wrong then it just outputs "Error" instead of crashing
  
  possfit_model = possibly(.f = fit_pathfinder_static_normaldata, otherwise = "Error")
  
  #test that it works
  #test  = possfit_model(data_list2[[1]])
  
  #run simulations!
  results_normal <- future_map(data_list2, ~possfit_model(.x), .progress = TRUE, .options = furrr_options(seed = TRUE))
  
  trialwise_datanormal = map_dfr(results_normal, bind_rows)
  
  
  possfit_model = possibly(.f = fit_pathfinder_static_logisticdata, otherwise = "Error")
  
  #test that it works
  #test  = possfit_model(data_list2[[1]])
  
  #run simulations!
  #results_log <- future_map(data_list2, ~possfit_model(.x), .progress = TRUE, .options = furrr_options(seed = TRUE))
  
  #trialwise_datalog = map_dfr(results_log, bind_rows)
  
  
  #trialwise_datalog = trialwise_datalog %>% mutate(model = "logistic")

  
  trialwise_datanormal = trialwise_datanormal %>% mutate(model = "normal")
  
  #trialwise_data = rbind(trialwise_datalog, trialwise_datanormal)

  trialwise_data = trialwise_datanormal
  
  trialwise_data$simulated_alpha = alpha  
  trialwise_data$simulated_beta = beta
  trialwise_data$simulated_lapse = lambda
  
  
  return(trialwise_data)
}




get_ICC_psychometric = function(trialwise_data){
  
  sim_n_id = rnorm(1,0,10)
  
  
  write.csv(trialwise_data, here::here("data","Model recovery",paste0("subs = ", max(trialwise_data$subs),
                                                                " trials = ", trialwise_data$trials[1],
                                                                " id = ", rnorm(1,0,10),".csv")))

  
  
  mod_logs = cmdstanr::cmdstan_model(here::here("Stanmodels","Model_comparison_logistic.stan"))
  
  mod_norm = cmdstanr::cmdstan_model(here::here("Stanmodels","Model_comparison_normal.stan"))
  
  mod_gumbel = cmdstanr::cmdstan_model(here::here("Stanmodels","Model_comparison_gumbel.stan"))
  
  datanormal = trialwise_data %>% filter(model == "normal")
  
  
  df = datanormal %>% group_by(subs,X) %>% summarize(yn = sum(resp), n = n())
  

  datastan = list(Y = df$yn,
                  N = nrow(df),
                  npx = df$n,
                  S = length(unique(df$subs)),
                  S_id = df$subs,
                  X = df$X)
    
  
  fit_logs_data_norm <- mod_logs$sample(
    data = datastan,
    chains = 4,
    refresh = 500,
    init = 0,
    parallel_chains = 4,
    adapt_delta = 0.9,
    max_treedepth = 12)
  
  
  
  fit_gumb_data_norm <- mod_gumbel$sample(
    data = datastan,
    chains = 4,
    refresh = 500,
    init = 0,
    parallel_chains = 4,
    adapt_delta = 0.9,
    max_treedepth = 12)
  
  
  
  fit_norm_data_norm <- mod_norm$sample(
    data = datastan,
    chains = 4,
    refresh = 500,
    init = 0,
    parallel_chains = 4,
    adapt_delta = 0.9,
    max_treedepth = 12)
  
  
  loo_datanorm = data.frame(loo::loo_compare(list(logs = fit_logs_data_norm$loo(),norm = fit_norm_data_norm$loo(), gumbel = fit_gumb_data_norm$loo()))) %>% 
    mutate(trials = unique(trialwise_data$trials),
           subs = length(unique(trialwise_data$subs)),
           sim_id = sim_n_id,
           mean_div_normal = mean(fit_norm_data_norm$diagnostic_summary("divergences")$num_divergent),
           mean_div_logistic = mean(fit_logs_data_norm$diagnostic_summary("divergences")$num_divergent),
           mean_tree_normal = mean(fit_norm_data_norm$diagnostic_summary("treedepth")$num_max_treedepth),
           mean_tree_logistic = mean(fit_logs_data_norm$diagnostic_summary("treedepth")$num_max_treedepth),
           mean_div_gumbel = mean(fit_gumb_data_norm$diagnostic_summary("divergences")$num_divergent),
           mean_tree_gumbel = mean(fit_gumb_data_norm$diagnostic_summary("treedepth")$num_max_treedepth))
  
  individual_esti = fit_norm_data_norm$summary(c("gm","tau_u","alpha","beta")) %>% mutate(trials = unique(trialwise_data$trials),
                                                                                             subs = length(unique(trialwise_data$subs)),
                                                                                             sim_id = sim_n_id)
  
  simulated_esti = datanormal %>% dplyr::select(-c(X,prob,resp)) %>% distinct() %>% mutate(trials = unique(trialwise_data$trials),
                                                                                           subs = length(unique(trialwise_data$subs)),
                                                                                           sim_id = sim_n_id)
  
  
  # datalogs = trialwise_data %>% filter(model == "logistic")
  # 
  # df = datalogs %>% group_by(subs,X) %>% summarize(yn = sum(resp), n = n())
  # 
  # 
  # datastan = list(Y = df$yn,
  #                 N = nrow(df),
  #                 npx = df$n,
  #                 S = length(unique(df$subs)),
  #                 S_id = df$subs,
  #                 X = df$X)
  # 
  # fit_logs_data_logs <- mod_logs$sample(
  #   data = datastan,
  #   chains = 4,
  #   refresh = 500,
  #   init = 0,
  #   parallel_chains = 4,
  #   adapt_delta = 0.9,
  #   max_treedepth = 12)
  # 
  # 
  # fit_norm_data_logs <- mod_norm$sample(
  #   data = datastan,
  #   chains = 4,
  #   refresh = 500,
  #   init = 0,
  #   parallel_chains = 4,
  #   adapt_delta = 0.9,
  #   max_treedepth = 12)
  # 
  # loo_datalog = data.frame(loo::loo_compare(list(logs = fit_logs_data_logs$loo(),norm = fit_norm_data_logs$loo()))) %>%
  #   mutate(trials = unique(trialwise_data$trials),
  #          subs = length(unique(trialwise_data$subs)),
  #          sim_id = sim_n_id)
  # 
  # 
  loos = list(normaldata = loo_datanorm)
  
  return(list(loos, individual_esti, simulated_esti))
  
}


together = function(parameters){
  iccs = get_ICC_psychometric(get_params(parameters$subs,parameters$trials, parameters$beta, parameters$lamda, parameters$alpha))
  return(list(iccs))
  
}
