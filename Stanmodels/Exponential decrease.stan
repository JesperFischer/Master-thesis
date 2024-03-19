
data{
  
  int N;
  int N_subs;
  
  array[N] int sub_id;
  
  vector[N] trials;
  
  vector[N] means_between;
  vector[N] means_within;
  vector[N] means_residual;
}

parameters{
  
  // real a;
  vector[N_subs] b_uncon_between;
  vector[N_subs] sigma_uncon_between;
  
  vector[N_subs] b_uncon_within;
  vector[N_subs] sigma_uncon_within;
  
  vector[N_subs] b_uncon_resid;
  vector[N_subs] sigma_uncon_resid;
  
}

transformed parameters{
  
  vector[N_subs] b_between = b_uncon_between;
  vector[N_subs] sigma_between = exp(sigma_uncon_between);

  vector[N_subs] b_within = exp(b_uncon_within);
  vector[N_subs] sigma_within = exp(sigma_uncon_within);

  vector[N_subs] b_resid = exp(b_uncon_resid);
  vector[N_subs] sigma_resid = exp(sigma_uncon_resid);
  
  
  vector[N] mu_within;
  vector[N] mu_resid;
  
  for(n in 1:(N)){
    mu_within[n] = exp(b_within[sub_id[n]]*(1/trials[n]))-1;
    mu_resid[n] = exp(b_resid[sub_id[n]]*(1/trials[n]))-1;
  }
}



model{
  
  b_uncon_between ~ normal(10,10);
  sigma_uncon_between ~ normal(0,3);

  b_uncon_within ~ normal(0,3);
  sigma_uncon_within ~ normal(0,3);


  b_uncon_resid ~ normal(0,3);
  sigma_uncon_resid ~ normal(0,3);
  
  
  for(n in 1:N){
    
    means_between[n] ~ normal(b_between[sub_id[n]],sigma_between[sub_id[n]]);
    
    means_within[n] ~ normal(mu_within[n],sigma_within[sub_id[n]]);
    
    means_residual[n] ~ normal(mu_resid[n],sigma_resid[sub_id[n]]);
  }
}



generated quantities{
  
  
  vector[N] ICC;
  
  for(n in 1:N){
    
    ICC[n] = b_between[sub_id[n]] / (b_between[sub_id[n]] + mu_within[n] + mu_resid[n]);
      
  }
  
  
  
}
