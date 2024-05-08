
data {
  int<lower=1> N_s1; // total number of observations (all subjects) (session1)
  int<lower=1> N_s2; // total number of observations (all subjects) (session2)
                      
  int<lower=1> C;  // total number of sessions (all subjects) (for now 2)
  int<lower=0> S;  //Total number of subjects

  array[max(N_s1,N_s2),C] int S_id;

  array[max(N_s1,N_s2),C] int Y;
  
  matrix[max(N_s1,N_s2), C] X;  // design matrix (first column being intercept i.e. 1)

}
transformed data{
  int N = 6;
  
}


parameters {
  // hierarchical group level means 
  vector [N] gm;

  // hierarchical group level deviations
  vector<lower = 0>[N]  tau_u;

  // Subject-level estimate matrix 
  matrix[N, S] z_expo;
  
  
  // for the cholesky decomposition
  cholesky_factor_corr[N] L_u;
  
}



transformed parameters{
  //getting subject level estimates for ease.
  matrix[max(N_s1,N_s2),C] X_normal = rep_matrix(0, max(N_s1,N_s2),C);
  matrix[max(N_s1,N_s2),C] phi = rep_matrix(0, max(N_s1,N_s2),C);
  matrix[max(N_s1,N_s2),C] mu = rep_matrix(0, max(N_s1,N_s2),C);
  
  
  
  matrix[S,2] psyalpha = rep_matrix(0, S,2);

  matrix<lower = 0>[S,2] psybeta = rep_matrix(0, S,2);

  matrix<lower = 0, upper = 1>[S,2] lapse = rep_matrix(0, S,2);

  
  matrix[S, N] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  matrix[S, N] param;
  
  // adding the participant deviation to the group means
  for(n in 1:N){
    param[,n]= gm[n] + indi_dif[,n];
  }
  
  
  psyalpha[,1] = param[,1];
  psyalpha[,2] = param[,4];
  
  psybeta[,1] = exp(param[,2]);
  psybeta[,2] = exp(param[,5]);
  
  lapse[,1] = inv_logit(param[,3]) / 2;
  lapse[,2] = inv_logit(param[,6]) / 2;
  


  
  
    for (n in 1:(max(N_s1,N_s2))){
     for (c in 1:C){
      X_normal[n,c] = 0.5+0.5*erf((X[n,c]-psyalpha[S_id[n,c],c])/(psybeta[S_id[n,c],c]*sqrt(2)));
      //probability of responding 1.
      phi[n,c] = lapse[S_id[n,c],c] + (1 - 2 * lapse[S_id[n,c],c]) * X_normal[n,c];
    }
  }
  
  
}


model{
  
  //priors
  target += normal_lpdf(gm[1] | 0,20);
  
  target += normal_lpdf(gm[2] | 0,5);
  
  target += normal_lpdf(gm[3] | -3,5);
  
  target += normal_lpdf(gm[4] | 0,20);

  target += normal_lpdf(gm[5] | 0,5);
  
  target += normal_lpdf(gm[6] | -3,5);


  target += std_normal_lpdf(to_vector(z_expo));
    target += lkj_corr_cholesky_lpdf(L_u | 2);

  
  target += normal_lpdf(tau_u[1] | 0, 10)-normal_lccdf(0 | 0, 10);

  target += normal_lpdf(tau_u[2] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[3] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[4] | 0, 10)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[5] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[6] | 0, 3)-normal_lccdf(0 | 0, 3);


  for (n in 1:(max(N_s1,N_s2))){
   for (c in 1:C){
     target += bernoulli_lpmf(Y[n,c] | phi[n,c]);

    }
  }
}


generated quantities{
  
  real prior_mu_cond;
  
  matrix[max(N_s1,N_s2),C] log_lik;
  
  matrix[N,N] correlation_matrix;
  
  
  correlation_matrix = L_u*L_u';
  
  

  
  for (n in 1:(max(N_s1,N_s2))){
   for (c in 1:C){
     log_lik[n,c] = bernoulli_lpmf(Y[n,c] | phi[n,c]);

    }
  }
  
  
  
}
