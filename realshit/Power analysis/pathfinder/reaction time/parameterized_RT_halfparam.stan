data{

  //Constants
  int<lower=1> T; // Total number of trials in the data
  int<lower=1> S; // Total number of subjects in the data
  array[T] int S_id; //n vector of integeres that signify participant numbers
  
  array[T] int Y; // Vector of binary responses

  vector[T] X; // Vector of deltaBPM values that match the binary response  


  matrix[T,2] X_alpha;
  
  matrix[T,2] X_beta;
  

  vector[S] min_RT; // Vector of deltaBPM values that match the binary response  




  vector<lower=0>[T] RT; // Vector of deltaBPM values that match the binary response
  

}
transformed data{
  int<lower=1> N=9; //Number of free parameters

}
parameters{
  
  vector[N] gm;  // Group means 

  vector<lower = 0>[N]  tau_u;   // Between participant scales

  cholesky_factor_corr[N] L_u;    // Between participant cholesky decomposition

  matrix[N, S] z_expo;    // Participant deviation from the group means

}
transformed parameters{
  ///Recomposition of parameters (i.e. calculating the linear combination of parameters and design matrix.)


  // vector[S] psyalpha;
  // vector[S] psyalpha_dif
  
  matrix[S,2] alpha_p;
  
  // vector[S] psybeta;
  // vector[S] psybeta_dif;
  // 
  matrix[S,2] beta_p;

  
  vector[S] lapse;
  
  vector[S] betart;
  vector[S] intercept;
  vector[S] sigma;
  vector[S] ndt;


  vector[T] trial_psyalpha;
  vector[T] trial_psybeta;

profile("parameters") {

  // Extracting individual deviations for each subject for each parameter
  matrix[S, N] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  matrix[S, N] param;
  
  // adding the participant deviation to the group means
  for(n in 1:N){
    param[,n]= gm[n] + indi_dif[,n];
  }
  
    
  alpha_p = param[,1:2];
  
  beta_p = param[,3:4];

  lapse = inv_logit(param[,5]) / 2;
  
  intercept = param[,6];
  
  betart = exp(param[,7]);
  
  sigma = exp(param[,8]);
  
  ndt = inv_logit(param[,9]) .* min_RT;
  
}
  
  
}

model{
  
  // Defining priors.
  

  vector[T] prob;  
  
  vector[T] rt_mu;
  vector[T] rt_sd;

  vector[T] ndts;

  gm[1] ~ normal(0,20); //global mean of beta
  
  gm[2] ~ normal(0,20); //global mean of beta
  
  gm[3] ~ normal(0,3); //global mean of lapse

  gm[4] ~ normal(0,3); //global mean of alpha

  gm[5] ~ normal(-4, 2); //global mean of beta
  
  gm[6] ~ normal(0,3); //global mean of alpha
  
  gm[7] ~ normal(0,3); //global mean of alpha

  gm[8] ~ normal(0,3); //global mean of alpha

  gm[9] ~ normal(-4,5); //global mean of alpha
  

  to_vector(z_expo) ~ std_normal();
  
  tau_u[1] ~ normal(10 , 10);
  tau_u[2] ~ normal(10 , 10);
  
  tau_u[3] ~ normal(0 , 3);
  tau_u[4] ~ normal(0, 3);
  tau_u[5] ~ normal(0, 3);
  tau_u[6] ~ normal(0, 3);
  tau_u[7] ~ normal(0, 3);
  tau_u[8] ~ normal(0, 3);
  tau_u[9] ~ normal(0, 3);
    
    
  L_u ~ lkj_corr_cholesky(2);
  
  // Computing the likelihood. The cummulative normal is used here:
  

  profile("loop") {
  for(n in 1:T){

    prob[n] = lapse[S_id[n]] + (1 - 2 * lapse[S_id[n]]) * ((0.5+0.5 * erf((X[n]-dot_product(X_alpha[n,], alpha_p[S_id[n],])) / (exp(dot_product(X_beta[n,], beta_p[S_id[n],])) * sqrt(2)))));
    rt_mu[n] = intercept[S_id[n]] + betart[S_id[n]]*prob[n]*(1-prob[n]);
    rt_sd[n] = sigma[S_id[n]];
    
    ndts[n] = ndt[S_id[n]];
    }
    
  }
  
  profile("likelihood") {
  Y ~ bernoulli(prob);
  RT ~ lognormal(rt_mu - ndts , rt_sd);
  }
  
}

generated quantities{
  // Calculating the correlation matrix from the cholesky decomposition.
  real centered  = 0;
  matrix[N,N] correlation_matrix = L_u*L_u';

  
}
