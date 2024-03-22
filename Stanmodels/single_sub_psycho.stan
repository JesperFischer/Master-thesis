data {
  int<lower=0> N;
  array[N] int <lower=0> y;
  vector[N] x;
}
parameters {
  real alpha_unconstrained;
  real beta_unconstrained;
  real lambda_unconstrained;
}


transformed parameters{
  
  real alpha = alpha_unconstrained;
  real<lower=0> beta = exp(beta_unconstrained);
  real<lower=0,upper=0.5> lambda = inv_logit(lambda_unconstrained) / 2;
  
  
}


model {
  
  for (i in 1:N){
        y[i] ~ bernoulli(lambda  + (1 - 2 * lambda) * (0.5+0.5*erf((x[i]-alpha)/(beta*sqrt(2)))));
  }
  
  alpha_unconstrained ~ normal(0,20);
  beta_unconstrained ~ normal(0,2);
  lambda_unconstrained ~ normal(-4,2);
  
    
}
