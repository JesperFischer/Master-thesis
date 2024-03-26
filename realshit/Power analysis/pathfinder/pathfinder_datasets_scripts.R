power_analysis_v2 = function(parameters){
  
  source(here::here("realshit","Power analysis", "Make_datasets_scripts.R"))
  source(here::here("realshit","Power analysis","pathfinder", "pathfinder_datasets_scripts.R"))
  
  a = get_sub_paramers(parameters = data.frame(subjects = parameters$subjects,
                                               effect_size_alpha = parameters$effect_size_alpha,
                                               effect_size_beta = parameters$effect_size_beta
  ))
  
  
  effectsizedata_beta = a %>% group_by(sessions) %>% summarize(mean = mean(beta), sd = sd(beta))
  
  change_in_beta = (effectsizedata_beta[1,2]-effectsizedata_beta[2,2])[[1]]
  change_in_sd_cohens_beta = ((effectsizedata_beta[1,2]-effectsizedata_beta[2,2])/(sqrt((effectsizedata_beta[1,3]^2 + effectsizedata_beta[2,3]^2)/2)))[[1]]
  change_in_sd_cohens_beta
  
  effectsizedata_alpha = a %>% group_by(sessions) %>% summarize(mean = mean(alpha), sd = sd(alpha))
  
  change_in_alpha = (effectsizedata_alpha[1,2]-effectsizedata_alpha[2,2])[[1]]
  change_in_sd_cohens_alpha = -((effectsizedata_alpha[1,2]-effectsizedata_alpha[2,2])/(sqrt((effectsizedata_alpha[1,3]^2 + effectsizedata_alpha[2,3]^2)/2)))[[1]]
  change_in_sd_cohens_alpha
  
  #getting Psi stimuli
  df = a %>% arrange(sessions,participant_id) %>% mutate(trials = parameters$trials)
  
  df$alpha_obs_power = change_in_sd_cohens_alpha
  df$beta_obs_power = change_in_sd_cohens_beta
  
  print("Timer for Pathfinder")
  tictoc::tic()
  
  ########################################################################################################### GETTING PSI!
  df = df %>% mutate(n = 1:n())
  
  data_list2 <-  split(df, df$n)
  
  plan(multisession, workers = 3)
  
  #adding safety for if something goes wrong then it just outputs "Error" instead of crashing
  
  possfit_model = possibly(.f = get_pathfinder_stim, otherwise = "Error")
  
  results_normal <- future_map(data_list2, ~possfit_model(.x), .progress = TRUE, .options = furrr_options(seed = TRUE))
  
  data = map_dfr(results_normal, bind_rows)
  

  ###########################################################################################################
  tictoc::toc()
  
  #wrangling the data
  
  data = inner_join(df,
                    data %>% 
                      dplyr::select(resp,X,participant_id,sessions),
                    by = c("participant_id","sessions")
  )
  
  
  
  data = data %>% mutate(iter = parameters$id)%>% 
    mutate(real_effectsize_alpha = parameters$effect_size_alpha,
           real_effectsize_beta = parameters$effect_size_beta,
           subjects = parameters$subjects,
           trials = parameters$trials,
           parameter = parameters$parameter,
           observed_effectsize_alpha = change_in_sd_cohens_alpha,
           observed_effectsize_beta = change_in_sd_cohens_beta)
  
  directory = paste0("effect_size_alpha = " , parameters$effect_size_alpha,
                     " effect_size_beta = ", parameters$effect_size_beta,
                     " trials = ", df$trials[1],
                     " subjects = ", length(unique(df$participant_id))[1],
                     " iter = ", parameters$id[1],
                     " random id = ", round(rnorm(1,0,1000),2),".csv")
  
  
  
  pattern = paste0("effect_size_alpha = ",
                   as.character(parameters$effect_size_alpha[1]),
                   " effect_size_beta = ",
                   as.character(parameters$effect_size_beta[1]),
                   " trials = ",as.character(parameters$trials[1]),
                   " subjects = ",as.character(parameters$subjects[1]))
  
  if(!dir.exists(here::here("realshit","Power analysis","pathfinder","datasets",pattern))){
    dir.create(here::here("realshit","Power analysis","pathfinder","datasets",pattern))
  }
  
  
  write.csv(data,here::here("realshit","Power analysis","pathfinder","datasets",pattern,directory))
  
  if(parameters$psi_only == T){
    print("done with")
    print(parameters$id)
    return("done")
  }else{
    
    
    fit_model = function(data){
      
      source(here::here("Fitting functions","scripts","utilities.R"))
      
      mod_noncent = cmdstanr::cmdstan_model(here::here("Power analysis","Stan models","Simulations","Parameterized cummulative normal.stan"),stanc_options = list("O1"))
      
      mod_cent = cmdstanr::cmdstan_model(here::here("Power analysis","Stan models","Simulations","Parameterized cummulative normal_centered.stan"),stanc_options = list("O1"))
      
      data = data %>% 
        mutate(session = ifelse(sessions == 1, 0 ,1))
      
      data = transform_data_to_stan(data)
      
      
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
      
      
      
      
      #start with centered parameterization
      if(nrow(data) > 1500){
        
        #fitting
        fit_centered <- mod_cent$sample(
          data = data_stan,
          iter_sampling = 1000,
          iter_warmup = 1000,
          chains = 4,
          parallel_chains = 4,
          refresh = 500,
          adapt_delta = 0.8,
          max_treedepth = 10
        )
        
        diags_centered = data.frame(fit_centered$diagnostic_summary())
        rhat_centered = data.frame(fit_centered$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                                          "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% .$rhat
        
        
        
        if(max(diags_centered$num_divergent) >= 1 | max(rhat_centered) > 1.03){
          #fitting
          fit_nocentered <- mod_noncent$sample(
            data = data_stan,
            iter_sampling = 1000,
            iter_warmup = 1000,
            chains = 4,
            parallel_chains = 4,
            refresh = 500,
            adapt_delta = 0.8,
            max_treedepth = 10
          )
          
          diags_nocentered = data.frame(fit_nocentered$diagnostic_summary())
          rhat_nocentered = data.frame(fit_nocentered$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                                                "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% .$rhat
          
          if (max(diags_centered) > max(diags_nocentered)) {
            fit_norm <- fit_nocentered
          } else if (max(diags_centered) < max(diags_nocentered)) {
            fit_norm <- fit_centered
          } else {
            if (max(rhat_centered) > max(rhat_nocentered)) {
              fit_norm <- fit_nocentered
            } else if (max(rhat_centered) < max(rhat_nocentered)) {
              fit_norm <- fit_centered
            } else {
              fit_norm <- fit_nocentered
            }
          }
          
        }else{
          fit_norm <- fit_centered
        }
        
      }
      #start with noncentered parameterization
      if(nrow(data) <= 1500){
        
        fit_nocentered <- mod_noncent$sample(
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
        
        diags_nocentered = data.frame(fit_nocentered$diagnostic_summary())
        rhat_nocentered = data.frame(fit_nocentered$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                                              "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% .$rhat
        
        if(max(diags_nocentered$num_divergent) >= 1 | max(rhat_nocentered) > 1.03){
          #fitting
          
          fit_centered <- mod_cent$sample(
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
          
          diags_centered = data.frame(fit_centered$diagnostic_summary())
          rhat_centered = data.frame(fit_centered$summary(c("gm[1]","gm[2]","gm[3]","gm[4]","gm[5]",
                                                            "tau_u[1]","tau_u[2]","tau_u[3]","tau_u[4]","tau_u[5]"))) %>% .$rhat
          
          if (max(diags_centered) > max(diags_nocentered)) {
            fit_norm <- fit_nocentered
          } else if (max(diags_centered) < max(diags_nocentered)) {
            fit_norm <- fit_centered
          } else {
            if (max(rhat_centered) > max(rhat_nocentered)) {
              fit_norm <- fit_nocentered
            } else if (max(rhat_centered) < max(rhat_nocentered)) {
              fit_norm <- fit_centered
            } else {
              fit_norm <- fit_nocentered
            }
          }
          
        }else{
          fit_norm <- fit_nocentered
        }
      }
      
      
      return(fit_norm)
    }
    
    
    data_rows = transform_data_to_stan(data)
    
    
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
    
    
    rows = data_rows %>% group_by(participant_id, sessions) %>% summarize(n = n()) %>% ungroup() %>% mutate(rows = cumsum(n)) %>% .$rows
    
    
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
  
}



erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


get_moving_estimates = function(effectsize, sd, mean, correlation, correlation_upper, correlation_lower, subs){
  
  
  var_alpha1 = sd^2
  var_alpha2 = 1.5*var_alpha1
  
  mu_difference = effectsize * sqrt((var_alpha2 + var_alpha1) / 2)
  
  sd_difference = sqrt(var_alpha2-var_alpha1)
  
  library(faux)
  d = data.frame(correlation = -5)
  
  while(!(d$correlation > correlation_lower & d$correlation < correlation_upper)){
    dd = rnorm_multi(n = subs, 
                     mu = c(mean,mean),
                     sd = c(sd, sd),
                     r = correlation, 
                     varnames = c("param_1", "param_2"),
                     empirical = FALSE) %>% 
      mutate(param_2 = param_2 + rnorm(subs, mu_difference, sd_difference))
    
    d$correlation = cor.test(dd$param_1, dd$param_2)$estimate
  }
  return(dd)
  
}

get_stationary_estimates = function(sd, mean, correlation, correlation_upper, correlation_lower, subs){
  
  library(faux)
  d = data.frame(correlation = -5)
  
  while(!(d$correlation > correlation_lower & d$correlation < correlation_upper)){
    dd = rnorm_multi(n = subs, 
                     mu = c(mean,mean),
                     sd = c(sd, sd),
                     r = correlation, 
                     varnames = c("param_1", "param_2"),
                     empirical = FALSE)
    
    d$correlation = cor.test(dd$param_1, dd$param_2)$estimate
  }
  return(dd)
  
}


get_sub_paramers = function(parameters){
  
  alphas = get_moving_estimates(effectsize = parameters$effect_size_alpha,
                                sd = 10.8,
                                mean = -8.3,
                                correlation = 0.52,
                                correlation_upper = 0.6,
                                correlation_lower = 0.4,
                                subs = parameters$subjects)
  
  
  betas = exp(get_moving_estimates(effectsize = parameters$effect_size_beta,
                                   sd = 0.31,
                                   mean = 2.3,
                                   correlation = 0.4,
                                   correlation_upper = 0.6,
                                   correlation_lower = 0.17,
                                   subs = parameters$subjects))
  
  
  
  lambdas = brms::inv_logit_scaled(get_stationary_estimates(
    sd = 1.3,
    mean = -4.2,
    correlation = 0,
    correlation_upper = 0.8,
    correlation_lower = -0.8,
    subs = parameters$subjects)) / 2
  
  
  lapse = data.frame(lambdas) %>% mutate(participant_id = 1:parameters$subjects) %>% 
    pivot_longer(c("param_1","param_2"), names_to = "sessions",values_to = "lapse")
  
  alpha = data.frame(alphas) %>% mutate(participant_id = 1:parameters$subjects) %>% 
    pivot_longer(c("param_1","param_2"), names_to = "sessions",values_to = "alpha")
  
  beta = data.frame(betas) %>% mutate(participant_id = 1:parameters$subjects) %>% 
    pivot_longer(c("param_1","param_2"), names_to = "sessions",values_to = "beta") 
  
  parameters2 = inner_join(inner_join(lapse,alpha,by = join_by(sessions, participant_id)),
                           beta,by = join_by(sessions, participant_id))%>% 
    mutate(id = 1:nrow(.),sessions = ifelse(sessions == "param_1", 1,ifelse(sessions == "param_2",2,NA)))
  
  return(parameters2)
}

get_pathfinder_stim = function(parameters){
  
  N = parameters$trials
  
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  parm_ev = data.frame()
  
  x[1] = 0
  p[1] = parameters$lapse + (1 - 2 * parameters$lapse) * (0.5+0.5*erf((x[1]-parameters$alpha)/(parameters$beta*sqrt(2))))
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
    p[i] = parameters$lapse + (1 - 2 * parameters$lapse) * (0.5+0.5*erf((x[i]-parameters$alpha)/(parameters$beta*sqrt(2))))
    r[i] = rbinom(1,1,p[i])
    
  }
  
  df_trial = data.frame(X = x,
                        prob = p,
                        resp = r,
                        trials = 1:N) %>% mutate(trials = parameters$trials,
                                                 lapse = parameters$lapse,
                                                 alpha = parameters$alpha,
                                                 beta = parameters$beta,
                                                 id =  rnorm(1,0,1000),
                                                 participant_id = parameters$participant_id,
                                                 sessions = parameters$sessions,
                                                 model = "normal")
  
  #df_trial = full_join(parm_ev,df_trial)
  
  
  return(list(df_trial))
  
  
  
  return(d)
  
}



#function to get confidence intervals used in plot_intervals.
ci = function(x){
  list = list(which(cumsum(x)/sum(x) > 0.025)[1],
              last(which(cumsum(x)/sum(x) < 0.975)))
  return(list)
  
}

#function to get the line intervals from the exteropost and interopost dataframe
get_line_intervals = function(data, parameter){
  
  
  if(parameter == "alpha"){
    confidence = ci(rowMeans(data))
  }else if(parameter == "beta"){
    confidence = ci(colMeans(data))
  }
  rg = seq(-50.5,50.5, by = 1)
  rg = seq(0.1, 25, by = 0.1)
  
  
  upper = rg[confidence[[1]]]
  lower = rg[confidence[[2]]]
  
  
  data = data.frame(upper = upper,lower = lower)
  
  return(data)
}

