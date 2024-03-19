
data{
  
  int N;
  int N_subs;
  
  array[N] int sub_id;
  vector[N] trials;
  vector[N] means;
  
}

parameters{
  
  // real a;
  vector[N_subs] b_uncon;
  vector[N_subs] sigma_uncon;
}

transformed parameters{
  
  vector[N_subs] b = exp(b_uncon);
  vector[N_subs] sigma = exp(sigma_uncon);
  
  
  vector[N] mu;
  
  for(n in 1:N){
    mu[n] = 1-exp(-b[sub_id[n]]*trials[n]);
  }
}



model{
  
  
  // a ~ normal(0,10);
  b_uncon ~ normal(0,3);
  sigma_uncon ~ normal(0,3);
  
  
  for(n in 1:N){
    means[n] ~ normal(mu[n],sigma[sub_id[n]]);
}
}
