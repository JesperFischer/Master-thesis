
data {
  int<lower=1> N_s1; // total number of observations (all subjects) (session1)
  int<lower=1> N_s2; // total number of observations (all subjects) (session2)
                      
  int<lower=1> C;  // total number of sessions (all subjects) (for now 2)
  int<lower=0> S;  //Total number of subjects

  matrix[S,C] min_RT;

  array[max(N_s1,N_s2),C] int S_id;

  array[max(N_s1,N_s2),C] int Y;
  
  matrix <lower = 0> [max(N_s1,N_s2), C] RT;
  
  matrix <lower = 0> [max(N_s1,N_s2), C] Conf;
  
  
  matrix[max(N_s1,N_s2), C] X;  // design matrix (first column being intercept i.e. 1)

}
transformed data{
  int N = 20;
  
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
  matrix[max(N_s1,N_s2),C] mu_conf = rep_matrix(0, max(N_s1,N_s2),C);
  
  
  
  matrix[S,2] psyalpha = rep_matrix(0, S,2);

  matrix<lower = 0>[S,2] psybeta = rep_matrix(0, S,2);

  matrix<lower = 0, upper = 1>[S,2] lapse = rep_matrix(0, S,2);
  
  matrix[S,2] intercept = rep_matrix(0, S,2);
  
  matrix[S,2] beta = rep_matrix(0, S,2);
  
  matrix<lower = 0>[S,2] sigma = rep_matrix(0, S,2);
  
  matrix<lower = 0 >[S,2] ndt = rep_matrix(0, S,2);
  
  matrix[S,2] intercept_conf = rep_matrix(0, S,2);
  
  matrix[S,2] beta_conf = rep_matrix(0, S,2);
  
  matrix[S,2] sigma_conf = rep_matrix(0, S,2);
  
  
  matrix[S, N] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  matrix[S, N] param;
  
  // adding the participant deviation to the group means
  for(n in 1:N){
    param[,n]= gm[n] + indi_dif[,n];
  }
  
  
  
  
  psyalpha[,1] = param[,1];
  psyalpha[,2] = param[,11];
  
  psybeta[,1] = exp(param[,2]);
  psybeta[,2] = exp(param[,12]);
  
  lapse[,1] = inv_logit(param[,3]) / 2;
  lapse[,2] = inv_logit(param[,13]) / 2;
  
  intercept[,1] = param[,4];
  intercept[,2] = param[,14];
  
  beta[,1] = exp(param[,5]);
  beta[,2] = exp(param[,15]);
  
  sigma[,1] = exp(param[,6]);
  sigma[,2] = exp(param[,16]);
  
  ndt[,1] = inv_logit(param[,7]) .* min_RT[,1];
  ndt[,2] = inv_logit(param[,17])  .* min_RT[,2];
  
  intercept_conf[,1] = param[,8];
  intercept_conf[,2] = param[,18];

  beta_conf[,1] = param[,9];
  beta_conf[,2] = param[,19];

  sigma_conf[,1] = exp(param[,10]);
  sigma_conf[,2] = exp(param[,20]);
  
  
    for (n in 1:(max(N_s1,N_s2))){
     for (c in 1:C){
      X_normal[n,c] = 0.5+0.5*erf((X[n,c]-psyalpha[S_id[n,c],c])/(psybeta[S_id[n,c],c]*sqrt(2)));
      //probability of responding 1.
      phi[n,c] = lapse[S_id[n,c],c] + (1 - 2 * lapse[S_id[n,c],c]) * X_normal[n,c];
      //mean for shiftedlognormal RT likelihood
      mu[n,c] = intercept[S_id[n,c],c] + beta[S_id[n,c],c]*phi[n,c]*(1-phi[n,c]);
      
      mu_conf[n,c] = inv_logit(intercept_conf[S_id[n,c],c] + beta_conf[S_id[n,c],c]*phi[n,c]*(1-phi[n,c]));
      
    }
  }
  
  
}


model{
  
  //priors
  target += normal_lpdf(gm[1] | 0,20);
  target += normal_lpdf(gm[11] | 0,20);
  
  target += normal_lpdf(gm[2] | 0,5);
  target += normal_lpdf(gm[12] | 0,5);
  
  target += normal_lpdf(gm[3] | -3,5);
  target += normal_lpdf(gm[13] | -3,5);
  
  target += normal_lpdf(gm[4] | 0,5);
  target += normal_lpdf(gm[14] | 0,5);

  target += normal_lpdf(gm[5] | 0,5);
  target += normal_lpdf(gm[15] | 0,5);

  target += normal_lpdf(gm[6] | 0,5);
  target += normal_lpdf(gm[16] | 0,5);

  target += normal_lpdf(gm[7] | -2,5);
  target += normal_lpdf(gm[17] | -2,5);

  target += normal_lpdf(gm[8] | 0,5);
  target += normal_lpdf(gm[18] | 0,5);

  target += normal_lpdf(gm[9] | 0,5);
  target += normal_lpdf(gm[19] | 0,5);

  target += normal_lpdf(gm[10] | 0,5);
  target += normal_lpdf(gm[20] | 0,5);



  target += std_normal_lpdf(to_vector(z_expo));
    target += lkj_corr_cholesky_lpdf(L_u | 2);

  
  target += normal_lpdf(tau_u[1] | 0, 10)-normal_lccdf(0 | 0, 10);
  target += normal_lpdf(tau_u[11] | 0, 10)-normal_lccdf(0 | 0, 10);

  target += normal_lpdf(tau_u[2] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[12] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[3] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[13] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[4] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[14] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[5] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[15] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[6] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[16] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[7] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[17] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[8] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[18] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[9] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[19] | 0, 3)-normal_lccdf(0 | 0, 3);

  target += normal_lpdf(tau_u[10] | 0, 3)-normal_lccdf(0 | 0, 3);
  target += normal_lpdf(tau_u[20] | 0, 3)-normal_lccdf(0 | 0, 3);



  for (n in 1:(max(N_s1,N_s2))){
   for (c in 1:C){
     target += bernoulli_lpmf(Y[n,c] | phi[n,c]);
     target += lognormal_lpdf(RT[n,c] - ndt[S_id[n,c],c] | mu[n,c], sigma[S_id[n,c],c]);
     target += beta_proportion_lpdf(Conf[n,c] | mu_conf[n,c], sigma_conf[S_id[n,c],c]);
     
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
     log_lik[n,c] += lognormal_lpdf(RT[n,c] - ndt[S_id[n,c],c] | mu[n,c], sigma[S_id[n,c],c]);
     log_lik[n,c] += beta_proportion_lpdf(Conf[n,c] | mu_conf[n,c], sigma_conf[S_id[n,c],c]);
     

    }
  }
  
  
  
}
