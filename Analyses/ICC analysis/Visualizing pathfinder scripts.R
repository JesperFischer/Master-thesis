

get_params_single = function(trials){
  
  alpha_uncon = rnorm(1,0,10)
  beta_uncon = rnorm(1,2,0.6)
  lambda_uncon = rnorm(1,-4,2)    
  
  
  
  parameters = data.frame(real_alpha_uncon = alpha_uncon,
                          real_lambda_uncon = lambda_uncon,
                          real_beta_uncon = beta_uncon,
                          trials = trials)
  
  
  parameters$real_beta_con = exp(parameters$real_beta_uncon)
  parameters$real_lambda_con = brms::inv_logit_scaled(parameters$real_lambda_uncon) / 2
  parameters$real_alpha_con = parameters$real_alpha_uncon
  
  
  
  pathfinder = get_pathfinder(parameters)
  
  uniform = get_uniform(parameters)
  
  psi = get_psi(parameters)
  
  return(list(pathfinder,uniform,psi))
}

# algorithm and function to obtain simulus and values as well as the fitted posterior distribution for the uniform approach
get_uniform = function(parameters){
  
  N = parameters$trials
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  mod = cmdstanr::cmdstan_model(here::here("Stanmodels","pathfinder.stan"))
  
  x = seq(-50,50,length.out = N)
  
  parm_ev = data.frame()
  for(i in 1:N){
    p[i] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (0.5+0.5*erf((x[i]-parameters$real_alpha_con)/(parameters$real_beta_con*sqrt(2))))
    r[i] = rbinom(1,1,p[i])
    
      
      data_stan = list(trials = nrow((data.frame(r[1:i]) %>% drop_na())),
                       x = data.frame(x[1:i]) %>% drop_na() %>% .$x,
                       r = data.frame(r[1:i]) %>% drop_na() %>% .$r)
      
      fit_psy <- mod$pathfinder(data = data_stan,refresh=0)
      
      
      c_trial = data.frame(fit_psy$summary(c("psy_alpha","psy_beta","psy_lambda"))) %>% mutate(trials = i)
      
      parm_ev = rbind(c_trial,parm_ev)
      
  }
  
  
  df_runif = data.frame(X = x,
                             prob = p,
                             resp = r,
                             real_lambda_con = parameters$real_lambda_con,
                             real_alpha_con = parameters$real_alpha_con,
                             real_beta_con = parameters$real_beta_con,
                             real_lambda_uncon = parameters$real_lambda_uncon,
                             real_alpha_uncon = parameters$real_alpha_uncon,
                             real_beta_uncon = parameters$real_beta_uncon,
                             id =  rnorm(1,0,1000),
                             model = "normal",
                             trials = 1:N)
  
  
  df_runif = inner_join(parm_ev,df_runif, by = "trials")
  
  reals = df_runif %>% dplyr::select(real_lambda_con,real_alpha_con,real_beta_con) %>% 
    rename(psy_lambda = real_lambda_con, psy_beta = real_beta_con, psy_alpha = real_alpha_con)%>% 
    pivot_longer(cols = c("psy_alpha","psy_beta","psy_lambda"), names_to = "variable")
  
  # df_runif %>% ggplot(aes(x = trials, y = mean, ymin = q5, ymax = q95, col = variable))+
  #   geom_pointrange() +
  #   facet_wrap(~variable, scales = "free")+
  #   geom_hline(data = reals, aes(yintercept = value))+
  #   theme_classic()
  # 
  
  df_runif$alpha_reals = unique(reals %>% filter(variable == "psy_alpha") %>% .$value)
  df_runif$beta_reals = unique(reals %>% filter(variable == "psy_beta") %>% .$value)
  df_runif$lambda_reals = unique(reals %>% filter(variable == "psy_lambda") %>% .$value)
  
return(df_runif)
  
}

# algorithm and function to obtain simulus and values as well as the fitted posterior distribution for the pathfinder approach
get_pathfinder = function(parameters){
  
  
  
  N = parameters$trials
  
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  parm_ev = data.frame()
  
  x[1] = 0
  p[1] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (0.5+0.5*erf((x[1]-parameters$real_alpha_con)/(parameters$real_beta_con*sqrt(2))))
  r[1] = rbinom(1,1,p[1])
  
  
  mod = cmdstanr::cmdstan_model(here::here("Stanmodels","pathfinder.stan"))

  tictoc::tic()
  for(i in 1:N){
    
    
    data_stan = list(trials = nrow((data.frame(r) %>% drop_na())),
                     x = data.frame(x) %>% drop_na() %>% .$x,
                     r = data.frame(r) %>% drop_na() %>% .$r)
    
    fit_psy <- mod$pathfinder(data = data_stan,refresh=0)
    
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
  t2 = tictoc::toc()
  
  df_pathfinder = data.frame(X = x,
                        prob = p,
                        resp = r,
                        real_lambda_con = parameters$real_lambda_con,
                        real_alpha_con = parameters$real_alpha_con,
                        real_beta_con = parameters$real_beta_con,
                        real_lambda_uncon = parameters$real_lambda_uncon,
                        real_alpha_uncon = parameters$real_alpha_uncon,
                        real_beta_uncon = parameters$real_beta_uncon,
                        id =  rnorm(1,0,1000),
                        model = "normal",
                        trials = 1:N)
  
  
  df_pathfinder = inner_join(parm_ev,df_pathfinder, by = "trials")
  
  reals = df_pathfinder %>% dplyr::select(real_lambda_con,real_alpha_con,real_beta_con) %>% 
    rename(psy_lambda = real_lambda_con, psy_beta = real_beta_con, psy_alpha = real_alpha_con)%>% 
    pivot_longer(cols = c("psy_alpha","psy_beta","psy_lambda"), names_to = "variable")
  
  # df_pathfinder %>% ggplot(aes(x = trials, y = mean, ymin = q5, ymax = q95, col = variable))+
  #   geom_pointrange() +
  #   facet_wrap(~variable, scales = "free")+
  #   geom_hline(data = reals, aes(yintercept = value))+
  #   theme_classic()
  
  
  df_pathfinder$alpha_reals = unique(reals %>% filter(variable == "psy_alpha") %>% .$value)
  df_pathfinder$beta_reals = unique(reals %>% filter(variable == "psy_beta") %>% .$value)
  df_pathfinder$lambda_reals = unique(reals %>% filter(variable == "psy_lambda") %>% .$value)
  
  return(list(df_pathfinder,t2))
  
}



# algorithm and function to obtain simulus and values as well as the fitted posterior distribution for the PSI approach
get_psi = function(parameters){
  
  python_script <- here::here("Python","Psi_script.py")
  
  alpha = parameters$real_alpha_con
  
  beta = parameters$real_beta_con
  
  lapse = parameters$real_lambda_con
  
  trials = as.integer(parameters$trials)
  

  library(reticulate)
  
  # Use reticulate to run the Python script with arguments
  
  source_python(python_script, convert = FALSE)
  tictoc::tic()
  d = get_stim(lapse, alpha, beta, trials)
  t2 = tictoc::toc()
  #Produces a warning about will be removed in future versions (but seems to be a thing with reticulate and not my code)
  dd = reticulate::py_to_r(d)
  dd = reticulate::py_to_r(d)%>% unnest()
  
  return(list(dd,t2))
    

  
}


# function that combines the above functions in order to compare the 3 different ways of computing the trial-by-trial stimulus values.
# This function also fits these trial-by-trial stimulus with the same model. Making the comparison between the methods fair.
get_params_single_multiple = function(params){
  
  sim_n_id = rnorm(1,0,10)
  
  alpha_uncon = rnorm(1,0,10)
  beta_uncon = rnorm(1,2,0.6)
  lambda_uncon = rnorm(1,-4,2)    
  
  parameters = data.frame(real_alpha_uncon = alpha_uncon,
                          real_lambda_uncon = lambda_uncon,
                          real_beta_uncon = beta_uncon,
                          trials = params$trials)
  
  
  parameters$real_beta_con = exp(parameters$real_beta_uncon)
  parameters$real_lambda_con = brms::inv_logit_scaled(parameters$real_lambda_uncon) / 2
  parameters$real_alpha_con = parameters$real_alpha_uncon
  
  
  
  pathfinder = get_pathfinder(parameters)
  
  timer = pathfinder[[2]]$toc - pathfinder[[2]]$tic
  pathfinder = pathfinder[[1]]
  
  mod = cmdstanr::cmdstan_model(here::here("Stanmodels","single_sub_psycho.stan"))
  
  
  datastan = list(N = nrow(pathfinder),
                  x = pathfinder$X,
                  y = pathfinder$resp)
  
  
  fit = mod$sample(data = datastan,
                   chains = 4,
                   refresh = 500,
                   iter_warmup = 1000,
                   iter_sampling = 1000,
                   parallel_chains = 4,
                   adapt_delta = 0.90,
                   max_treedepth = 12)
  
  
  
  df_pathfinder = data.frame(fit$summary(c("alpha","beta","lambda","alpha_unconstrained","beta_unconstrained","lambda_unconstrained"))) %>% 
    mutate(trials = unique(parameters$trials),
           sim_id = sim_n_id,
           real_lambda_con = unique(pathfinder$real_lambda_con),
           real_lambda_uncon = unique(pathfinder$real_lambda_uncon),
           real_beta_con = unique(pathfinder$real_beta_con),
           real_beta_uncon = unique(pathfinder$real_beta_uncon),
           real_alpha_con = unique(pathfinder$real_alpha_con),
           real_alpha_uncon = unique(pathfinder$real_alpha_uncon),
           id = unique(pathfinder$id),
           timer = timer,
           mean_div = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent)) %>% .$mean_div,
           tree_depth = data.frame(fit$diagnostic_summary()) %>% summarize(mean_treedepth = mean(num_max_treedepth)) %>% .$mean_treedepth)
  
  
  uniform = get_uniform(parameters)
  
  
  datastan = list(N = nrow(uniform),
                  x = uniform$X,
                  y = uniform$resp)
  
  
  fit = mod$sample(data = datastan,
                   chains = 4,
                   refresh = 500,
                   iter_warmup = 1000,
                   iter_sampling = 1000,
                   parallel_chains = 4,
                   adapt_delta = 0.90,
                   max_treedepth = 12)
  
  
  
  df_uniform = data.frame(fit$summary(c("alpha","beta","lambda","alpha_unconstrained","beta_unconstrained","lambda_unconstrained"))) %>% 
    mutate(trials = unique(parameters$trials),
           sim_id = sim_n_id,
           real_lambda_con = unique(uniform$real_lambda_con),
           real_lambda_uncon = unique(uniform$real_lambda_uncon),
           real_beta_con = unique(uniform$real_beta_con),
           real_beta_uncon = unique(uniform$real_beta_uncon),
           real_alpha_con = unique(uniform$real_alpha_con),
           real_alpha_uncon = unique(uniform$real_alpha_uncon),
           id = unique(uniform$id),
           timer = NA,
           mean_div = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent)) %>% .$mean_div,
           tree_depth = data.frame(fit$diagnostic_summary()) %>% summarize(mean_treedepth = mean(num_max_treedepth)) %>% .$mean_treedepth)
  
  
  
  
  psi = get_psi(parameters)
  
  
  timer = psi[[2]]$toc - psi[[2]]$tic
  psi = psi[[1]]
  
  datastan = list(N = nrow(psi),
                  x = psi$X,
                  y = psi$resp)
  
  
  fit = mod$sample(data = datastan,
                   chains = 4,
                   refresh = 500,
                   iter_warmup = 1000,
                   iter_sampling = 1000,
                   parallel_chains = 4,
                   adapt_delta = 0.90,
                   max_treedepth = 12)
  
  
  
  df_psi = data.frame(fit$summary(c("alpha","beta","lambda","alpha_unconstrained","beta_unconstrained","lambda_unconstrained"))) %>% 
    mutate(trials = unique(parameters$trials),
           sim_id = sim_n_id,
           real_lambda_con = unique(psi$lapse),
           real_lambda_uncon = NA,
           real_beta_con = unique(psi$beta),
           real_beta_uncon = NA,
           real_alpha_con = unique(psi$alpha),
           real_alpha_uncon = NA,
           id = NA,
           timer = timer,
           mean_div = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent)) %>% .$mean_div,
           tree_depth = data.frame(fit$diagnostic_summary()) %>% summarize(mean_treedepth = mean(num_max_treedepth)) %>% .$mean_treedepth)
  
  
  
  return(list(df_pathfinder, df_uniform, df_psi))
  
}
