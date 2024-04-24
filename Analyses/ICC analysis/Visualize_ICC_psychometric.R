
Correlations = function(results, string){
  
  results2 <- unlist(results, recursive = FALSE)
  
  params_correlation <- map_dfr(results2, 3) %>% mutate(string = string)
  
  return(params_correlation)
  
}



trialwise = function(results, parametera){
  
  results2 <- unlist(results, recursive = FALSE)
  
  trial <- map_dfr(results2, 2)
  
  
  q = trial %>% filter(parameter == parametera) %>% dplyr::select(-c("resp","prob","X","trials.1")) %>% distinct() %>% 
    pivot_longer(cols = c("real_beta_uncon","real_lambda_uncon","real_alpha_uncon")) %>% 
    filter(paste0("real_",parameter, "_uncon") == name)
  
  return(q)
}


ICC_raw = function(results, string){
  
  results2 <- unlist(results, recursive = FALSE)
  
  ICC_raw <- map_dfr(results2, 1) %>% mutate(string = string) 
  
  return(ICC_raw)
  
}


ICC = function(results,string){
  
  results2 <- unlist(results, recursive = FALSE)
  
  # Step 2: Use map_dfr to row-bind all dataframes
  params_icc <- map_dfr(results2, 1)
  
  
  
  i = params_icc %>% filter(!grepl("ICC",variable)) %>% 
    unnest() %>% dplyr::select(variable,mean,sd,q5,q95, residual_variance, trials,subs, sim_id,mean_div ) %>% 
    rename(parameter = variable, uncertainty = residual_variance) %>% mutate(correlation = F)
  
  
  ii_real = i %>% mutate(color = ifelse(grepl("alpha",parameter),"alpha",ifelse(grepl("beta",parameter),"beta","lambda")),
                         facet = ifelse(grepl("resid",parameter),"resid",ifelse(grepl("between",parameter),"between","within")),
  )
  
  
  parameters = c("alpha","beta","lambda")
  i = 0
  
  ICCS = list()
  for(para in parameters){
    i = i+1
    iia = ii_real %>% dplyr::filter(color == para)
    
    ii = iia %>% arrange(subs,trials) %>% dplyr::select(facet, subs, trials, mean, sd)
    
    ##ii = ii %>% group_by(subs,trials,facet) %>% summarize(mean = mean(mean), sd = mean(sd))
    
    ii2 = ii %>% rename(mu = mean, sigma = sd) %>% pivot_wider(names_from = "facet",values_from = c(mu,sigma)) %>% unnest() %>% 
      mutate(mu_between = mu_between^2, mu_within = mu_within^2) 
      # group_by(subs, trials) %>% 
      # summarize(mu_resid = mean(mu_resid),
      #           mu_between  = mean(mu_between),
      #           mu_within   = mean(mu_within),
      #           sigma_resid  = mean(sigma_resid ),
      #           sigma_between = mean(sigma_between ),
      #           sigma_within = mean(sigma_within)
      # )
    
    # ii2 %>% mutate(ICC = mu_between / (mu_between+ mu_within + mu_resid)) %>% ggplot(aes(x = trials, y = ICC, col = as.factor(subs)))+geom_point()
    # 
    # params_icc %>% filter(variable == "ICC_alpha", residual_variance == T) %>% arrange(subs,trials) %>% ggplot(aes(x = trials, y = mean, col = as.factor(subs)))+geom_point()
    
    
    mod = cmdstanr::cmdstan_model(here::here("Stanmodels","Exponential decrease measurement.stan"))
    
    datastan = list(N = nrow(ii2), 
                    trials = ii2$trials,
                    means_residual_obs = ii2$mu_resid,
                    means_between_obs = ii2$mu_between,
                    means_within_obs = ii2$mu_within,
                    sds_residual = ii2$sigma_resid,
                    sds_between = ii2$sigma_between,
                    sds_within = ii2$sigma_within,
                    sub_id = as.numeric(as.factor(ii2$subs)),
                    N_subs = length(unique(ii2$subs)))
    
    fit = mod$sample(data = datastan,
                     chains = 4,
                     refresh = 500,
                     iter_warmup = 2000,
                     iter_sampling = 2000,
                     parallel_chains = 4,
                     adapt_delta = 0.90,
                     max_treedepth = 12)
    
    
    
    ICCS[i] = list(fit$summary("ICC") %>% mutate(subs = as.factor(ii2$subs), trials = ii2$trials) %>% mutate(parameter = para))
  }
  
  ICC_real = map_dfr(ICCS, bind_rows) %>% mutate(string = string) 
  
  return(ICC_real)
  
  }