


# getting stimulus values of agent based on the uniform approach (without fitting)
fit_static = function(parameters){
  
  N = parameters$trials
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  
  x = seq(-50,50,length.out = N)
  
  p = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (0.5+0.5*erf((x-parameters$real_alpha_con)/(parameters$real_beta_con*sqrt(2))))
  
  r = rbinom(N,1,p)
  
  
  df_trial = data.frame(X = x,
                        prob = p,
                        resp = r,
                        trials = parameters$trials,
                        real_lambda_con = parameters$real_lambda_con,
                        real_alpha_con = parameters$real_alpha_con,
                        real_beta_con = parameters$real_beta_con,
                        real_lambda_uncon = parameters$real_lambda_uncon,
                        real_alpha_uncon = parameters$real_alpha_uncon,
                        real_beta_uncon = parameters$real_beta_uncon,
                        id =  rnorm(1,0,1000),
                        model = "normal",
                        trials = 1:N)
  
  
  
  return(df_trial)
}

#given trials gives stimulus values of agent.
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
  
  
  trialwise_data = fit_static(parameters)
  

  return(trialwise_data)
}



# fits the single subject psychometric function and returns a dataframe with posterior parameters
get_ICC_psychometric_singe_fit = function(trialwise_data){
  
  sim_n_id = rnorm(1,0,10)
  
  #mod = cmdstanr::cmdstan_model(here::here("report","DDM","Stan Models","Psychometric_2_ICC.stan"))
  mod = cmdstanr::cmdstan_model(here::here("Stanmodels","single_sub_psycho.stan"))
  

  datastan = list(N = nrow(trialwise_data),
                  x = trialwise_data$X,
                  y = trialwise_data$resp)
  

  fit = mod$sample(data = datastan,
                   chains = 4,
                   refresh = 500,
                   iter_warmup = 1000,
                   iter_sampling = 1000,
                   parallel_chains = 4,
                   adapt_delta = 0.90,
                   max_treedepth = 12)
  
  
  
  df = data.frame(fit$summary(c("alpha","beta","lambda","alpha_unconstrained","beta_unconstrained","lambda_unconstrained"))) %>% 
            mutate(trials = unique(trialwise_data$trials),
                   sim_id = sim_n_id,
                   real_lambda_con = unique(trialwise_data$real_lambda_con),
                   real_lambda_uncon = unique(trialwise_data$real_lambda_uncon),
                   real_beta_con = unique(trialwise_data$real_beta_con),
                   real_beta_uncon = unique(trialwise_data$real_beta_uncon),
                   real_alpha_con = unique(trialwise_data$real_alpha_con),
                   real_alpha_uncon = unique(trialwise_data$real_alpha_uncon),
                   id = unique(trialwise_data$id),
                   mean_div = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent)) %>% .$mean_div,
                   tree_depth = data.frame(fit$diagnostic_summary()) %>% summarize(mean_treedepth = mean(num_max_treedepth)) %>% .$mean_treedepth)
  
  return(list(df))
  
}



#function to combine the above functions. Takes the number of trials in a dataframe (parameters) 
# and returns the fitted posterior distribution of the parameters.
together = function(parameters){
  iccs = get_ICC_psychometric_singe_fit(get_params_single(parameters$trials))
  return(list(iccs))
  
}



