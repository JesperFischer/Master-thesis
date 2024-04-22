



fit_pathfinder_static_v2 = function(parameters){
  
  N = parameters$trials
  
  
  r = array(NA,N)
  p = array(NA,N)
  x = array(NA,N)
  
  parm_ev = data.frame()
  
  x[1] = 0
  p[1] = parameters$real_lambda_con + (1 - 2 * parameters$real_lambda_con) * (0.5+0.5*erf((x[1]-parameters$real_alpha_con)/(parameters$real_beta_con*sqrt(2))))
  r[1] = rbinom(1,1,p[1])
  
  
  mod = cmdstanr::cmdstan_model(here::here("report","DDM","Stan Models","pathfinder.stan"))
  
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
                        subs = parameters$subs,
                        sessions = parameters$sessions,
                        model = "normal",
                        trials = 1:N)
  
  #df_trial = full_join(parm_ev,df_trial)
  
  
  return(list(df_trial))
}


get_params = function(subs, trials){
  
  alpha_uncon = rnorm(subs,0,10)
  beta_uncon = rnorm(subs,3,0.6)
  lambda_uncon = rnorm(subs,-4,2)    
  
  beta_uncon = rnorm(subs,3,0.6)
  beta_uncon = rnorm(subs,3,0.6)
  
  
  
  replicate = 1:2
  
  parameters1 = NULL
  
  
  
  parameters = data.frame(real_alpha_uncon = alpha_uncon,
                          real_lambda_uncon = lambda_uncon,
                          real_beta_uncon = beta_uncon,
                          trials = trials)
  
  
  for(i in 1:(length(replicate))){
    parameters1 = rbind(parameters1,parameters %>% mutate(sessions = i))
  }
  
  parameters = parameters1 %>% 
    mutate(id = 1:nrow(.)) %>% 
    mutate(subs = rep(1:subs,2))
  
  parameters$real_beta_con = exp(parameters$real_beta_uncon)
  parameters$real_lambda_con = brms::inv_logit_scaled(parameters$real_lambda_uncon) / 2
  parameters$real_alpha_con = parameters$real_alpha_uncon
  
  data_list2 <- split(parameters, parameters$id)
  
  
  plan(multisession, workers = 4)
  
  #adding safety for if something goes wrong then it just outputs "Error" instead of crashing
  
  possfit_model = possibly(.f = fit_pathfinder_static_v2, otherwise = "Error")
  
  #test that it works
  #test  = possfit_model(data_list2[[1]])
  
  #run simulations!
  results <- future_map(data_list2, ~possfit_model(.x), .progress = TRUE, .options = furrr_options(seed = TRUE))
  
  
  trialwise_data = map_dfr(results, bind_rows)
  
  # write.csv(trialwise_data, here::here("datasets",paste0("subs = ", max(trialwise_data$subs),
  #                                                      " trials = ", trialwise_data$trials[1],
  #                                                      " id = ", rnorm(1,0,10),".csv")))
  # 
  # 
  return(trialwise_data)
}




get_ICC_psychometric = function(trialwise_data){
  
  sim_n_id = rnorm(1,0,10)
  
  #  mod = cmdstanr::cmdstan_model(here::here("report","DDM","Stan Models","Psychometric_2_ICC.stan"))
  mod = cmdstanr::cmdstan_model(here::here("report","DDM","Stan Models","Psychometric_2_ICC_MSE.stan"))
  
  
  # 
  # write.csv(trialwise_data, here::here("beta=3 datasets",paste0("subs = ", max(trialwise_data$subs),
  #                                                               " trials = ", trialwise_data$trials[1],
  #                                                               " id = ", rnorm(1,0,10),".csv")))
  # 
  
  
  params = trialwise_data %>% arrange(sessions, subs)
  
  sesTsubs = length(unique(params$subs))*length(unique(params$sessions))
  
  data_stan = list(T = nrow(params %>% filter(subs == 1 & sessions == 1)),
                   P = length(unique(params$subs)),
                   S = length(unique(params$sessions)),
                   resp = matrix(params %>% arrange(sessions, subs) %>% .$resp, nrow = sesTsubs, byrow = T),
                   X = matrix(params %>% arrange(sessions, subs) %>% .$X, nrow = sesTsubs, byrow = T),
                   
                   real_alpha = params %>% group_by(sessions, subs) %>% summarize(alpha = mean(real_alpha_uncon)) %>% .$alpha,
                   real_lambda = params %>% group_by(sessions, subs) %>% summarize(lambda = mean(real_lambda_uncon)) %>% .$lambda,
                   real_beta = params %>% group_by(sessions, subs) %>% summarize(beta = mean(real_beta_uncon)) %>% .$beta
                   
                   
  )
  
  #trialwise_data <- read_csv("datasets/subs = 14 trials = 90 id = 1.78080543202527.csv")
  # trialwise_data %>% pivot_longer(cols = c("real_lambda","real_alpha","real_beta")) %>% group_by(name) %>% summarize(mean = mean(value), sd = sd(value))
  # 
  # vector = fit$summary("beta") %>% .$mean
  # pairs <- lapply(1:(length(vector)/2), function(i) vector[c(i, i+14)])
  # pairwise_sd <- function(pair) {
  #   sd(pair)
  # }
  # 
  # # Calculate standard deviations of pairs
  # pairwise_sds <- unlist(lapply(pairs, pairwise_sd))
  # mean(pairwise_sds)
  
  fit <- mod$sample(
    data = data_stan,
    chains = 4,
    refresh = 500,
    init = 0,
    parallel_chains = 4,
    adapt_delta = 0.9,
    max_treedepth = 12)
  
  
  variances = rbind(data.frame(fit$summary("tau_b")) %>% 
                      mutate(variable = c("between_sub_var_beta","between_sub_var_lambda","between_sub_var_alpha")),
                    data.frame(fit$summary("tau_w")) %>% 
                      mutate(variable = c("within_sub_var_beta","within_sub_var_lambda","within_sub_var_alpha"))
  )
  
  
  residual = data.frame(fit$summary(c("resid_alpha_var","resid_beta_var","resid_lambda_var")))
  
  extras = rbind(residual,variances)
  
  
  df = rbind(rbind(data.frame(fit$summary("ICC")) %>% mutate(variable = c("ICC_beta","ICC_lambda","ICC_alpha")) %>% 
                     mutate(residual_variance = F),
                   rbind(fit$summary("ICC_alpha"),fit$summary("ICC_lambda"),fit$summary("ICC_beta"))%>% 
                     mutate(variable = c("ICC_alpha","ICC_lambda","ICC_beta"),residual_variance = T)
  ),extras %>% mutate(residual_variance = NA)) %>% 
    mutate(diag = data.frame(fit$diagnostic_summary()) %>% summarize(mean_div = mean(num_divergent)),
           tree_depth = data.frame(fit$diagnostic_summary()) %>% summarize(mean_treedepth = mean(num_max_treedepth)),
           trials = unique(params$trials),
           subs = length(unique(params$subs)),
           sessions = length(unique(params$sessions)),
           sim_id = sim_n_id)
  
  
  indi_estimates = rbind(data.frame(fit$summary(c("alpha","beta","lambda")))) %>% 
    mutate(id_model = gsub(".*?\\[(\\d+)\\].*?", "\\1", variable),
           parameter = gsub("\\[.*", "", variable),
           subs = rep(rep(1:length(unique(params$subs)), 3),2),
           sessions = rep(rep(1:2, each = length(unique(params$subs))),3),
           n_subs = length(unique(params$subs)),
           sim_id = sim_n_id)
  
  trial_estimates = inner_join(indi_estimates, params) %>% distinct()
  
  #correlations:
  correlation_no_unc = trial_estimates %>% dplyr::select(-c("X","resp","prob","trials.1")) %>% distinct() %>% 
    pivot_longer(cols = c("real_lambda_uncon","real_alpha_uncon","real_beta_uncon"), names_to = "real_parameternames",values_to = "real_values") %>% 
    filter(paste0("real_",parameter,"_uncon") == real_parameternames) %>% 
    group_by(parameter) %>% 
    summarize(correlation_mean = cor.test(mean, real_values)$estimate[1],
              correlation_q2 = cor.test(mean, real_values)$conf.int[[1]],
              correlation_q97 = cor.test(mean, real_values)$conf.int[[2]]
    )%>% mutate(ucnertainty_added = F)
  
  samples = 2000
  
  correlation_unc = trial_estimates %>% dplyr::select(-c("X","resp","prob","trials.1")) %>% distinct() %>% 
    pivot_longer(cols = c("real_lambda_uncon","real_alpha_uncon","real_beta_uncon"), names_to = "real_parameternames",values_to = "real_values") %>% 
    filter(paste0("real_",parameter,"_uncon") == real_parameternames) %>% 
    rowwise() %>% 
    mutate(new_estimates = list(rnorm(samples, mean,sd)),
           id = list(1:samples)) %>%  unnest() %>% group_by(id, parameter) %>%
    summarize(correlation = cor.test(mean, new_estimates)$estimate[1]) %>% 
    group_by(parameter) %>% 
    summarize(correlation_mean = mean(correlation),
              correlation_q2 = HDInterval::hdi(correlation)[[1]],
              correlation_q97 = HDInterval::hdi(correlation)[[2]]
    ) %>% mutate(ucnertainty_added = T)
  
  
  
  correlations = rbind(correlation_no_unc,correlation_unc) %>% 
    mutate(variable = paste0("ICC_",parameter),
           trials = unique(params$trials),
           subs = length(unique(params$subs)),
           sessions = length(unique(params$sessions)),
           sim_id = sim_n_id)
  
  #plot:
  # trial_estimates %>% dplyr::select(-c("X","resp","prob","trials.1")) %>% distinct() %>% pivot_longer(cols = c("real_lambda","real_alpha","real_beta"), names_to = "real_parameternames",values_to = "real_values") %>% filter(paste0("real_",parameter) == real_parameternames) %>% ggplot(aes(x = mean, y = real_values))+geom_point()+facet_wrap(~parameter, scales = "free")
  
  
  
  
  return(list(df, trial_estimates,correlations))
  
}


together = function(parameters){
  iccs = get_ICC_psychometric(get_params(parameters$subs,parameters$trials))
  return(list(iccs))
  
}
