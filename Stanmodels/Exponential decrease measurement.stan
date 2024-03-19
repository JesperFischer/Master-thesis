
data{
  
  int N;
  int N_subs;
  
  array[N] int sub_id;
  
  vector[N] trials;
  
  vector[N] means_between_obs;
  vector[N] means_within_obs;
  vector[N] means_residual_obs;
  
  vector[N] sds_between;
  vector[N] sds_within;
  vector[N] sds_residual;
  
}

parameters{
  
  //vector[N] means_between;
  // vector[N] means_within;
  // vector[N] means_residual;


  // real a;
  vector[N_subs] b_uncon_between;
  vector[N_subs] sigma_uncon_between;
  
  vector[N_subs] b_uncon_within_mu;
  vector[N_subs] b_uncon_within_sigma;

  vector[N_subs] b_uncon_resid_mu;
  vector[N_subs] b_uncon_resid_sigma;
  
}

transformed parameters{
  
  vector[N_subs] b_within_mu = exp(b_uncon_within_mu);
  vector[N_subs] b_within_sigma = exp(b_uncon_within_sigma);


  vector[N_subs] b_resid_mu = exp(b_uncon_resid_mu);
  vector[N_subs] b_resid_sigma = exp(b_uncon_resid_sigma);

  vector[N] b_between;
  vector[N] sigma_between;


  vector<lower=0>[N] mu_within;
  vector<lower=0>[N] mu_resid;

  vector<lower=0>[N] sigma_resid;
  vector<lower=0>[N] sigma_within;



  for(n in 1:(N)){
    mu_within[n] = exp(b_within_mu[sub_id[n]]*(1/trials[n]))-1;
    mu_resid[n] = exp(b_resid_mu[sub_id[n]]*(1/trials[n]))-1;
    
    
    sigma_within[n] = exp(b_within_sigma[sub_id[n]]*(1/trials[n]))-1;
    sigma_resid[n] = exp(b_resid_sigma[sub_id[n]]*(1/trials[n]))-1;
    
    b_between[n] = b_uncon_between[sub_id[n]];
    sigma_between[n] = exp(sigma_uncon_between[sub_id[n]]);
  }
  

  
}



model{
  
  b_uncon_between ~ normal(100,100);
  sigma_uncon_between ~ normal(0,5);

  b_uncon_within_mu ~ normal(0,10);
  b_uncon_within_sigma ~ normal(0,10);

  b_uncon_resid_mu ~ normal(0,10);
  b_uncon_resid_sigma ~ normal(0,10);

  
   // means_between_obs ~ normal(means_between, sds_between);
   // means_within_obs ~ normal(means_within, sds_within);
   // means_residual_obs ~ normal(means_residual, sds_residual);

  
  for(n in 1:N){
    
    means_between_obs[n] ~ normal(b_between[n],sigma_between[n]);

    means_within_obs[n] ~ normal(mu_within[n],sigma_within[n]);

    means_residual_obs[n] ~ normal(mu_resid[n],sigma_resid[n]);
  }
}



generated quantities{


  vector[N] ICC;

  for(n in 1:N){

    ICC[n] = b_between[n] / (b_between[n] + mu_within[n] + mu_resid[n]);

  }



}
