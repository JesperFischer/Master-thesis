

// The input data is a vector 'y' of length 'N'.
data {
  
  int param;

  int<lower=0> N;
  vector[N] x;
  array[N] int y;
  
  matrix[N,param] design_matrix;
  
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {


  vector [param] expo_alpha;
  vector [param] expo_beta;
  
  // real <lower = 0> int_alpha;
  // real <lower = 0> int_beta;

  real asym_alpha;
  real asym_beta;  
  
  real intercept_alpha;
  real intercept_beta;
  
}


transformed parameters {
  
  vector[N] X_norm;
  vector[N] alpha;

  vector[N] beta;
  
  real mu_intercept_alpha = exp(intercept_alpha);
  real mu_intercept_beta = exp(intercept_beta);
  
  real mu_asym_alpha = exp(asym_alpha);
  real mu_asym_beta = exp(asym_beta);
  
  
  vector[N] mu_expo_alpha = design_matrix * expo_alpha;
  vector[N] mu_expo_beta = design_matrix * expo_beta;

  // vector[N] mu_asym_alpha = exp(design_matrix1 * asym_alpha);
  
  // vector[N] mu_asym_beta = exp(design_matrix1 * asym_beta);


  alpha =  mu_intercept_alpha * exp(-mu_expo_alpha) + mu_asym_alpha;

  beta =  mu_intercept_beta * exp(-mu_expo_beta) + mu_asym_beta;
  
  X_norm = 1/(1+exp(-(1/beta) .* (x - alpha)));
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  target += normal_lpdf(expo_alpha | 1, 3);
  target += normal_lpdf(expo_beta | 1, 3);

  
  // target += normal_lpdf(int_alpha | 0, 3);
  // target += normal_lpdf(int_beta | 0, 3);
  // 
  target += normal_lpdf(intercept_alpha | 0, 3);
  target += normal_lpdf(intercept_beta | 0, 3);
  
  target += normal_lpdf(asym_alpha | -3, 3);
  target += normal_lpdf(asym_beta | -3, 3);
  

  
  
   y ~ bernoulli(X_norm);
}


generated quantities{

  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = bernoulli_lpmf(y[n] | X_norm[n]);
  }


}

