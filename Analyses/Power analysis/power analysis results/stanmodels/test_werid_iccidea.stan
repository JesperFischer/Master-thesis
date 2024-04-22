

// The input data is a vector 'y' of length 'N'.
data {
  
  int N;
  vector[N] trials;
  vector[N] y;
  vector[N] y_se;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {


  real expo_alpha;
  
  real sigma_un;
  real intercept_alpha;
  
}


transformed parameters {
  
  vector[N] X_norm;
  vector[N] alpha;

  vector[N] beta;
  
  real mu_intercept_alpha = exp(intercept_alpha);
  
  real sigma = exp(sigma_un);
  
  // vector[N] mu_expo_alpha = design_matrix * expo_alpha;
  // vector[N] mu_expo_beta = design_matrix * expo_beta;

  // vector[N] mu_asym_alpha = exp(design_matrix1 * asym_alpha);
  
  // vector[N] mu_asym_beta = exp(design_matrix1 * asym_beta);

  alpha =  mu_intercept_alpha * (trials ^ expo_alpha);
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  target += normal_lpdf(expo_alpha | -1, 3);

  // target += normal_lpdf(int_alpha | 0, 3);
  // target += normal_lpdf(int_beta | 0, 3);
  // 
  target += normal_lpdf(intercept_alpha | 2, 2);

  target += normal_lpdf(sigma_un | 2, 2);

  // target += normal_lpdf(asym_alpha | 0, 3);
  // target += normal_lpdf(asym_beta | 0, 3);
  

  
  
   y ~ normal(alpha,sqrt(y_se+square(sigma)));
   // y ~ normal(alpha,sigma);
   
   
}
