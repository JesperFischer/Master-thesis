data {
  int<lower=0> N;
  array[N] int<lower=0> n;
  array[N] int <lower=0> y;
  vector[N] x;
  vector[N] RT;
  real minRT;
  
  
}
parameters {
  real alpha;
  real beta_unconstrained;
  real intercept;
  real betart;
  real sigma_unconstrained;
  real ndt_unconstrained;
  
}

transformed parameters{
  
  real <lower=0> beta = exp(beta_unconstrained);
  real <lower=0> sigma = exp(sigma_unconstrained);
  real <lower=0> ndt = inv_logit(ndt_unconstrained) * minRT;
  
}

model {
  vector[N] p = 0.5+0.5*erf((x-alpha) / (beta*sqrt(2)));
  
  for (i in 1:N){
    y[i] ~ bernoulli(p[i]);
    (RT[i] - ndt) ~ lognormal(intercept+betart*p[i]*(1-p[i]), sigma);
  }
  
  alpha ~ normal(0, 20);
  beta_unconstrained ~ normal(0,3);
  intercept ~ normal(0, 3);
  betart ~ normal(0, 3);
  sigma_unconstrained ~ normal(0, 3);
  ndt_unconstrained ~ normal(-3, 3);
    
}
