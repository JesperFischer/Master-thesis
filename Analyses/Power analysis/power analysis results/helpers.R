get_df = function(subjects,trials,parameter){
  files = list.files(here::here("Analyses","Power analysis","power analysis results",parameter), recursive = T)
  
  realfiles = files[grepl(paste0("trials=",trials),files)]
  
  realfiles = realfiles[grepl(paste0("subjects=",subjects),realfiles)]
  
  if(length(realfiles) == 0){
    print("error")
    return(NA)
  }
  df = data.frame()
  for(i in 1:length(realfiles)){
    effectsize = str_sub(realfiles[i], nchar(realfiles[i])-6,nchar(realfiles[i])-4)
    result = readRDS(here::here("Analyses","Power analysis","power analysis results", parameter,realfiles[i]))
    tester = map_dfr(result, 1) %>% mutate(effectsize = as.numeric(effectsize), parameters = parameter)
    df = rbind(df,tester)
  }
  
  return(df)
}


get_full_df = function(parameter){
  dirs = list.files(here::here("Analyses","Power analysis","power analysis results",parameter))
  #dirs = dirs[!grepl("subjects=50",dirs)]
  #dirs = dirs[!grepl("subjects=10_trials=25",dirs)]
  
  matches <- str_match(dirs, "results_subjects=(\\d+)_trials=(\\d+)")
  
  # Convert matched values to numeric
  subjects <- as.numeric(matches[, 2])
  trials <- as.numeric(matches[, 3])
  
  
  df_list <- lapply(seq(length(dirs)), function(i) {
    get_df(subjects[i], trials[i], parameter)
  })
  
  return(do.call(rbind, df_list))
  
}


get_modeling_data = function(df, parameter, signi_level = 0.05){
  
  #getting the number of significant results in each:
  alpha = signi_level
  
  p_zero = paste0("p_",parameter)
  
  group_change = paste0("effectsize_",parameter)
  
  if(parameter == "alpha"){
    
    ff = df %>% filter(variable == "mu_g[4]") %>% 
      mutate(significant = ifelse(get(p_zero)<alpha, 1, 0),
             subjects = as.factor(subjects),
             trials = as.factor(trials)) %>% 
      select(group_change, significant, subjects, trials) %>% 
      rename(group_change = group_change)
  }
  
  
  
  if(parameter == "beta"){
    
    ff = df %>% filter(variable == "mu_g[5]") %>% 
      mutate(significant = ifelse(get(p_zero)<alpha, 1, 0),
             subjects = as.factor(subjects),
             trials = as.factor(trials)) %>% 
      select(group_change, significant, subjects, trials) %>% 
      rename(group_change = group_change)
    
  }
  
  ff = ff %>% arrange(trials, subjects)
  
  ff$trials = droplevels(ff$trials)
  ff$subjects = droplevels(ff$subjects)
  
  return(ff)
}



fit_results_independent_trials_and_subjects = function(df,parameter, signi_level = 0.05){
  
  
  ff = get_modeling_data(df,parameter, signi_level = signi_level)
  ff = ff %>% filter(group_change > 0)
  
  
  #making the contrast matrix:
  
  wide_df <- ff %>% arrange(trials, subjects)%>% dplyr::select(trials,subjects) %>% 
    ungroup() %>% mutate(index = 1:nrow(.), value = 1) %>% 
    pivot_wider(names_from = c(subjects,trials), values_from = value, values_fill = 0)
  
  #Removing idex for the pivotwider to work
  
  TnS = as.matrix(wide_df[,-1])
  

  # fitting in stan
  standata = list(x = (ff$group_change),
                  N = nrow(ff),
                  nT_T_nS = ncol(TnS),
                  TnS = as.matrix(TnS),
                  y = ff$significant)
  
  
  mod = cmdstanr::cmdstan_model(here::here("Analyses","Power analysis","power analysis results","stanmodels","indpendent trials and subject effects.stan"),stanc_options = list("O1"))
  
  
  fit_cond_cmd <- mod$sample(
    data = standata, 
    chains = 4, 
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  
  #looking at the parameter values of this psychometric to gain infomation on these
  
  ddalpha = data.frame(fit_cond_cmd$summary("alphas")) %>% 
    mutate(subject_trials = colnames(TnS)) %>% 
    separate(subject_trials, c("subjects","trials"), sep = "_", convert = TRUE) %>% 
    mutate(trials = as.factor(trials), subjects = as.factor(subjects), parameter = "alpha")
  
  
  ddbeta = data.frame(fit_cond_cmd$summary("betas")) %>% 
    mutate(subject_trials = colnames(TnS)) %>% 
    separate(subject_trials, c("subjects","trials"), sep = "_", convert = TRUE) %>% 
    mutate(trials = as.factor(trials), subjects = as.factor(subjects), parameter = "beta",
           mean = exp(mean), q5 = exp(q5), q95 = exp(q95))
  
  df = rbind(ddalpha,ddbeta)
  return(list(df,fit_cond_cmd))
}

fit_results_independent_trials_and_subjects_v2 = function(df,parameter, signi_level = 0.05){
  
  
  ff = get_modeling_data(df,parameter, signi_level = signi_level)
  
  ff = ff %>% filter(group_change > 0)
  
  
  #making the contrast matrix:
  
  wide_df <- ff %>% arrange(trials, subjects)%>% dplyr::select(trials,subjects) %>% ungroup() %>% mutate(index = 1:nrow(.), value = 1) %>% 
    pivot_wider(names_from = c(subjects,trials), values_from = value, values_fill = 0)
  
  #Removing idex for the pivotwider to work
  
  TnS = as.matrix(wide_df[,-1])
  
  
  # fitting in stan
  standata = list(x = (ff$group_change),
                  N = nrow(ff),
                  nT_T_nS = ncol(TnS),
                  TnS = as.matrix(TnS),
                  y = ff$significant)
  
  
  mod = cmdstanr::cmdstan_model(here::here("Analyses","Power analysis","power analysis results","stanmodels","indpendent trials and subject effects weibull.stan"),stanc_options = list("O1"))
  
  
  fit_cond_cmd <- mod$sample(
    data = standata, 
    chains = 4, 
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  
  #looking at the parameter values of this psychometric to gain infomation on these
  
  ddalpha = data.frame(fit_cond_cmd$summary("alphas")) %>% 
    mutate(subject_trials = colnames(TnS)) %>% 
    separate(subject_trials, c("subjects","trials"), sep = "_", convert = TRUE) %>% 
    mutate(trials = as.factor(trials), subjects = as.factor(subjects), parameter = "alpha",
           mean = exp(mean), q5 = exp(q5), q95 = exp(q95))
  
  
  ddbeta = data.frame(fit_cond_cmd$summary("betas")) %>% 
    mutate(subject_trials = colnames(TnS)) %>% 
    separate(subject_trials, c("subjects","trials"), sep = "_", convert = TRUE) %>% 
    mutate(trials = as.factor(trials), subjects = as.factor(subjects), parameter = "beta",
           mean = exp(mean), q5 = exp(q5), q95 = exp(q95))
  
  df = rbind(ddalpha,ddbeta)
  return(list(df,fit_cond_cmd))
}

fit_results_independent_trials_and_subjects_v3 = function(df,parameter, signi_level = 0.05){
  
  
  ff = get_modeling_data(df,parameter, signi_level = signi_level)
  
  ff = ff %>% filter(group_change > 0)
  
  
  #making the contrast matrix:
  
  wide_df <- ff %>% arrange(trials, subjects)%>% dplyr::select(trials,subjects) %>% ungroup() %>% mutate(index = 1:nrow(.), value = 1) %>% 
    pivot_wider(names_from = c(subjects,trials), values_from = value, values_fill = 0)
  
  #Removing idex for the pivotwider to work
  
  TnS = as.matrix(wide_df[,-1])
  
  
  # fitting in stan
  standata = list(x = (ff$group_change),
                  N = nrow(ff),
                  nT_T_nS = ncol(TnS),
                  TnS = as.matrix(TnS),
                  y = ff$significant)
  
  
  mod = cmdstanr::cmdstan_model(here::here("Analyses","Power analysis","power analysis results","stanmodels","independent trials and subjects logs.stan"),stanc_options = list("O1"))
  
  
  fit_cond_cmd <- mod$sample(
    data = standata, 
    chains = 4, 
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  
  #looking at the parameter values of this psychometric to gain infomation on these
  
  ddalpha = data.frame(fit_cond_cmd$summary("alphas")) %>% 
    mutate(subject_trials = colnames(TnS)) %>% 
    separate(subject_trials, c("subjects","trials"), sep = "_", convert = TRUE) %>% 
    mutate(trials = as.factor(trials), subjects = as.factor(subjects), parameter = "alpha")
  
  
  ddbeta = data.frame(fit_cond_cmd$summary("betas")) %>% 
    mutate(subject_trials = colnames(TnS)) %>% 
    separate(subject_trials, c("subjects","trials"), sep = "_", convert = TRUE) %>% 
    mutate(trials = as.factor(trials), subjects = as.factor(subjects), parameter = "beta",
           mean = exp(mean), q5 = exp(q5), q95 = exp(q95))
  
  df = rbind(ddalpha,ddbeta)
  return(list(df,fit_cond_cmd))
}



fit_exponential_decays = function(df,parameter, ntrials, indi_parameters, signi_level = 0.05, desired_power = 0.8){
  
  
  ff = get_modeling_data(df,parameter, signi_level = signi_level, desired_power = desired_power)
  
  ff = ff %>% filter(trials == ntrials)
  
  ff = droplevels(ff)
  
  standata = list(x = ff$group_change,
                  N = nrow(ff),
                  Subs = unique(as.numeric(as.character(ff$subjects))),
                  n_Subs = length(unique(as.numeric(as.character(ff$subjects)))),
                  Subs_id = as.numeric(as.factor(ff$subjects)),
                  Trials = ff$trials,
                  y = ff$significant)
  
  
  mod = cmdstanr::cmdstan_model(here::here("Analyses","Power analysis","power analysis results","stanmodels","exponential decay with asymptote.stan"),stanc_options = list("O1"))
  #mod = cmdstanr::cmdstan_model(here::here("realshit","Power analysis","pathfinder","reaction time","power_analysis_psi","stanmodels","exponential decay with asymptote_both.stan"),stanc_options = list("O1"))
  #mod = cmdstanr::cmdstan_model(here::here("realshit","Power analysis","pathfinder","reaction time","power_analysis_psi","stanmodels","exponential decay without asymptote.stan"),stanc_options = list("O1"))
  
  #mod = cmdstanr::cmdstan_model(here::here("Power analysis","Stan models","Visualization","exponential decay without asymptote.stan"),stanc_options = list("O1"))
  
  
  fit <- mod$sample(
    data = standata, 
    chains = 4, 
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  
  exp_decay = function(x,int,expo,asym){
    
    return(int * exp(-expo * x) + asym)
    
  }
  
  get_decay_line = function(fit,parameter, nsubs){
    
    params = fit$summary(c(paste0("int_",parameter),paste0("expo_",parameter),paste0("asym_",parameter)))
    
    params = params %>% select(variable,mean) %>% pivot_wider(names_from = variable, values_from = mean)
    
    params = params %>% 
      mutate(x = list(5:nsubs), y = list(exp_decay(5:nsubs,get(paste0("int_",parameter)), get(paste0("expo_",parameter)), get(paste0("asym_",parameter))))) %>% 
      unnest() %>% 
      mutate(parameter = parameter) %>% 
      dplyr::select(x,y,parameter)
    
  }
  
  
  dcalpha = get_decay_line(fit,"alpha",nsubs = max(df_alpha$subjects))
  dcbeta = get_decay_line(fit,"beta",nsubs = max(df_alpha$subjects))
  
  dfparameters = rbind(dcalpha,dcbeta)
  
  parameters = indi_parameters %>% filter(trials == ntrials) %>% 
    mutate(subjects = as.numeric(as.character(subjects))) %>% 
    ggplot()+
    geom_pointrange(aes(x = subjects, y = mean, ymin = q5, ymax = q95, col = trials), width = 0.2, position=position_dodge(width=0.3))+
    facet_wrap(~parameter, scales = "free")+
    geom_line(data = dfparameters, aes(x = x, y = y))
  
  
  # next we might that the model fitted individually to the significance on differing trials and subjects and 
  # get the desired power for each subject / trial combination and  plot this
  
  
  return(list(parameters, dfparameters, fit))
  
}




fit_exponential_decays = function(df,parameter, ntrials, indi_parameters, signi_level = 0.05, desired_power = 0.8){
  
  
  ff = get_modeling_data(df,parameter, signi_level = signi_level, desired_power = desired_power)
  
  ff = ff %>% filter(trials == ntrials)
  
  #ff = ff %>% filter(trials != 10) %>% filter(subjects != 100 & group_change > 0)
  
  #ff = ff %>% filter(trials %in% c(100,150)) %>% filter(subjects %in% c(10,20,40)& group_change > 0)
  
  
  ff = ff %>% filter(group_change > 0) %>% filter(!(trials == 100 & subjects == 40))
  
  #ff = ff %>% filter(group_change > 0)
  subs = scale(as.numeric(as.character(ff$subjects)))[,1]
  trials = scale(as.numeric(as.character(ff$trials)))[,1]
  int = scale(as.numeric(as.character(ff$subjects)))[,1] * scale(as.numeric(as.character(ff$trials)))[,1]
  
  ff %>% ggplot(aes(x = group_change, y = significant, col = subjects))+geom_smooth()+facet_wrap(~trials)
  
  
  design_matrix = data.frame(subs = subs,
                             trials = trials,
                             interaction = int)
  
  design_matrix1 = data.frame(subs = subs,
                             trials = trials,
                             interaction = int)
  
  
  
  
  design_matrix = data.frame(subs = as.numeric(as.character(ff$subjects)),
                             trials = as.numeric(as.character(ff$trials)),
                             int = as.numeric(as.character(ff$trials)) * as.numeric(as.character(ff$subjects)))
  
  
  
  
  design_matrix1 = data.frame(subs = as.numeric(as.character(ff$subjects)),
                             trials = as.numeric(as.character(ff$trials)),
                             int = as.numeric(as.character(ff$trials)) * as.numeric(as.character(ff$subjects)))
  
  
  
  
  standata = list(x = ff$group_change,
                  N = nrow(ff),
                  design_matrix = design_matrix,
                  design_matrix1 = design_matrix1,
                  param = ncol(design_matrix),
                  param1 = ncol(design_matrix1),
                  y = ff$significant)
  
  
  mod = cmdstanr::cmdstan_model(here::here("Analyses","Power analysis","power analysis results","stanmodels","exponential decay with asymptote_both.stan"),stanc_options = list("O1"))

  mod = cmdstanr::cmdstan_model(here::here("Analyses","Power analysis","power analysis results","stanmodels","weibull2.stan"),stanc_options = list("O1"))
  
  
  fit <- mod$sample(
    data = standata, 
    chains = 4,
    refresh = 50,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 8
  )
  
  mcmc_trace(fit$draws(c("expo_alpha","expo_beta","intercept_beta","intercept_alpha")))
  
  exp_decay = function(x,int,expo,asym){
    
    return(int * exp(-expo * x) + asym)
    
  }
  
  get_decay_line = function(fit,parameter, nsubs){
    
    params = fit$summary(c(paste0("int_",parameter),paste0("expo_",parameter),paste0("asym_",parameter)))
    
    params = params %>% select(variable,mean) %>% pivot_wider(names_from = variable, values_from = mean)
    
    params = params %>% 
      mutate(x = list(5:nsubs), y = list(exp_decay(5:nsubs,get(paste0("int_",parameter)), get(paste0("expo_",parameter)), get(paste0("asym_",parameter))))) %>% 
      unnest() %>% 
      mutate(parameter = parameter) %>% 
      dplyr::select(x,y,parameter)
    
  }
  
  
  dcalpha = get_decay_line(fit,"alpha",nsubs = max(df_alpha$subjects))
  dcbeta = get_decay_line(fit,"beta",nsubs = max(df_alpha$subjects))
  
  dfparameters = rbind(dcalpha,dcbeta)
  
  parameters = indi_parameters %>% filter(trials == ntrials) %>% 
    mutate(subjects = as.numeric(as.character(subjects))) %>% 
    ggplot()+
    geom_pointrange(aes(x = subjects, y = mean, ymin = q5, ymax = q95, col = trials), width = 0.2, position=position_dodge(width=0.3))+
    facet_wrap(~parameter, scales = "free")+
    geom_line(data = dfparameters, aes(x = x, y = y))
  
  
  # next we might that the model fitted individually to the significance on differing trials and subjects and 
  # get the desired power for each subject / trial combination and  plot this
  
  
  return(list(parameters, dfparameters, fit))
  
}




