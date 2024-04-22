

// The input data is a vector 'y' of length 'N'.
data {
  
  int<lower=0> N;
  
  vector[N] y;
  vector[N] y_se;
  
  vector[N] trials;
  
  
  
  
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  real alpha1;
  real beta1;
  real sigma1;
  
}


transformed parameters {
  
  vector[N] X_norm;
  

  real alpha = exp(alpha1);
  real beta = beta1;
  real sigma = exp(sigma1);
  
  
  // vector[N] mu_expo_alpha = design_matrix * expo_alpha;
  // vector[N] mu_expo_beta = design_matrix * expo_beta;

  // vector[N] mu_asym_alpha = exp(design_matrix1 * asym_alpha);
  
  // vector[N] mu_asym_beta = exp(design_matrix1 * asym_beta);


  X_norm = 1-alpha*(1/(trials^beta));

  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  target += normal_lpdf(alpha1 | -1, 3);
  target += normal_lpdf(beta1 | -1, 3);
  target += normal_lpdf(sigma1 | -1, 3);
  
  
   y ~ normal(X_norm,sqrt(y_se+square(sigma)));
}

