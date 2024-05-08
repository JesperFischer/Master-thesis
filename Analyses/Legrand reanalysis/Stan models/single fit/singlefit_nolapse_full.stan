data {
  int<lower=0> N;
  array[N] int<lower=0> n;
  array[N] int <lower=0> y;
  vector[N] x;
  vector[N] RT;
  real minRT;
  vector[N] conf;
  
  
  
}
parameters {
  real alpha;
  real beta_unconstrained;
  real intercept;
  real betart;
  real sigma_unconstrained;
  real ndt_unconstrained;
  
  real intercept_conf;
  real beta_conf;
  real sigma_conf_unconstrained;

  
  
}

transformed parameters{
  
  real <lower=0> beta = exp(beta_unconstrained);
  real <lower=0> sigma = exp(sigma_unconstrained);
  real <lower=0> ndt = inv_logit(ndt_unconstrained) * minRT;
  
  real <lower=0> sigma_conf = exp(sigma_conf_unconstrained);
  
  
}

model {
  vector[N] p = 0.5+0.5*erf((x-alpha) / (beta*sqrt(2)));
  
  for (i in 1:N){
    y[i] ~ bernoulli(p[i]);
    (RT[i] - ndt) ~ lognormal(intercept+betart*p[i]*(1-p[i]), sigma);
    conf[i] ~ beta_proportion(inv_logit(intercept_conf+beta_conf*p[i]*(1-p[i])), sigma_conf);
  }
  
  alpha ~ normal(0, 20);
  beta_unconstrained ~ normal(0,3);
  intercept ~ normal(0, 3);
  betart ~ normal(0, 3);
  sigma_unconstrained ~ normal(0, 3);
  ndt_unconstrained ~ normal(-3, 3);
    
  intercept_conf ~ normal(0, 3);
  beta_conf ~ normal(0, 3);
  sigma_conf_unconstrained ~ normal(0, 3);
    
    
    
}
