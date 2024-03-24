data { 
  
  int<lower=0> trials;
  
  array[trials] real x;
  
  array[trials] int r;
  } 


parameters { 
  real psy_alpha;
  real<lower=0> psy_beta_uncon; 
  real<lower=0, upper = 0.5> psy_lambda;
  } 





transformed parameters{ 
  
  array[trials] real phi;
  
  real psy_beta = exp(psy_beta_uncon);
  for(i in 1:trials){
    phi[i] = psy_lambda + (1 - 2 * psy_lambda) *  (1/(1+exp(-psy_beta*(x[i]-psy_alpha))));
    }
    } 






model { 
  
  target += normal_lpdf(psy_alpha | 0,20);
  
  target += normal_lpdf(psy_beta_uncon | 0, 3);
  
  target += beta_proportion_lpdf(psy_lambda | 0.05, 5) / 2;
  
  
  for(i in 1:trials){
    target += bernoulli_lpmf(r[i] | phi[i]); 
    }
    
}

