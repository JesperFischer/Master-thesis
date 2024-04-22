

// The input data is a vector 'y' of length 'N'.
data {
  
  int param;
  int param1;

  int<lower=0> N;
  vector[N] x;
  array[N] int y;
  
  matrix[N,param] design_matrix;
  matrix[N,param1] design_matrix1;
  
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {


  vector [param] expo_alpha;
  vector [param] expo_beta;
  
  // real <lower = 0> int_alpha;
  // real <lower = 0> int_beta;

  vector [param1] asym_alpha;
  vector [param1] asym_beta;
  
  real intercept;
}


transformed parameters {
  
  vector<lower=0>[N] X_norm;
  vector<lower=0>[N] alpha;

  vector<lower=0>[N] beta;
  
  vector[N] mu_expo_alpha = exp(design_matrix * expo_alpha);
  vector[N] mu_expo_beta = exp(design_matrix * expo_beta);

  vector[N] mu_asym_alpha = exp(design_matrix1 * asym_alpha);
  vector[N] mu_asym_beta = exp(design_matrix1 * asym_beta);
  
  real<lower=0> intercept_mu = exp(intercept);
  
  
  alpha =  intercept_mu * exp(-mu_expo_alpha) + mu_asym_alpha;

  beta =  intercept_mu * exp(-mu_expo_beta) + mu_asym_beta;
  
  
  X_norm = 1-exp(-(x ./ alpha)^beta);
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  target += normal_lpdf(expo_alpha | 0, 5);
  target += normal_lpdf(expo_beta | 0, 5);
  
  target += normal_lpdf(intercept | 0, 5);
  
  
  // target += normal_lpdf(int_alpha | 0, 3);
  // target += normal_lpdf(int_beta | 0, 3);
  // 
  target += normal_lpdf(asym_alpha | 0, 3);
  target += normal_lpdf(asym_beta | 0, 3);
  

  
  
   y ~ bernoulli(X_norm);
}

