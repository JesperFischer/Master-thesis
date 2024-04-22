

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] x;
  array[N] int y;
  // 
  int n_Subs;
  vector[n_Subs] Subs;
  array[N] int Subs_id;
  
  vector[N] Trials;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {


  real expo_alpha_un;
  real expo_beta_un;

  real int_alpha_un;
  real int_beta_un;

  real asym_alpha_un;
  real asym_beta_un;

}


transformed parameters {
  
  vector[N] X_norm;
  vector[n_Subs] alpha;
  
  vector<lower=0>[n_Subs] beta;
  
  real <lower = 0> expo_alpha = exp(expo_alpha_un);
  real <lower = 0> expo_beta= exp(expo_beta_un);

  real <lower = 0> int_alpha= exp(int_alpha_un);
  real <lower = 0> int_beta= exp(int_beta_un);

  real <lower = 0> asym_alpha= exp(asym_alpha_un);
  real <lower = 0> asym_beta= exp(asym_beta_un);
  
  
  

  alpha =  int_alpha * exp(-expo_alpha * Subs) + asym_alpha;

  beta =  int_beta * exp(-expo_beta * Subs) + asym_beta;
  
  for(i in 1:N){
    
    
    X_norm[i] = 0.5+0.5 * erf((x[i]-alpha[Subs_id[i]]) ./ (beta[Subs_id[i]] * sqrt(2)));
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  target += normal_lpdf(expo_alpha_un | 0, 10);
  target += normal_lpdf(expo_beta_un | 0, 10);
  
  target += normal_lpdf(int_alpha_un | 0, 10);
  target += normal_lpdf(int_beta_un | 0, 10);
  
  target += normal_lpdf(asym_alpha_un | 0, 10);
  target += normal_lpdf(asym_beta_un | 0, 10);
  

   y ~ bernoulli(X_norm);
}

