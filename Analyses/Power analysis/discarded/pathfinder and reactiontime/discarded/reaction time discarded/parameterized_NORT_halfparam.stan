data{
  //Constants
  int<lower=1> T; //n Trials
  int<lower=1> S; //n participants
  array[T] int S_id; //n participants
  
  int<lower=1> N_alpha;
  int<lower=1> N_beta;
  int<lower=1> N_lapse;
  

  matrix[T,N_alpha] X_alpha;
  matrix[T,N_beta] X_beta;
  matrix[T,N_lapse] X_lapse;
  
  
  array[T] int Y;
  vector[T] X;
  
  array[T] int npx;  // total number of observations per X

  

}
transformed data{
  int<lower=1> N=N_beta+N_lapse+N_alpha; //n free parameters
  int<lower=1> N_centered=N_beta+N_alpha; //n free parameters

}
parameters{
  // Group means 
  vector[N] gm;
  // Between participant scales
  vector<lower = 0>[N]  tau_u;
  // Between participant cholesky decomposition
  cholesky_factor_corr[N_centered] L_u;
  // Participant deviation 
  matrix[N_lapse, S] z_expo;  
  
  array[S] vector[N_centered] param;        // Centered individual parameters


}
transformed parameters{


  ///Recomposition

  vector[T] alpha;
  vector[T] beta;
  vector[T] lapse;
  
  
  matrix[N_lapse,S] lapse_p =  gm[N] + z_expo * tau_u[N];
  
  
  
  for(n in 1:T){

    alpha[n] = X_alpha[n,1] * param[S_id[n],1] + X_alpha[n,2] * param[S_id[n],2];
    

    beta[n] = exp(-2*(X_beta[n,1] * param[S_id[n],3] + X_beta[n,2] * param[S_id[n],4]));
    
    lapse[n] = inv_logit(dot_product(X_lapse[,n], lapse_p[,S_id[n]])) / 2;
    
    }
    
}

model{



  gm[1] ~ normal(0,20); //global mean of alpha

  gm[2] ~ normal(0,20); //global mean of alpha

  gm[3] ~ normal(0,3); //global mean of beta
  
  gm[4] ~ normal(0, 3); //global mean of beta
  
  gm[5] ~ normal(-4,2); //global mean of lapse


  to_vector(z_expo) ~ std_normal();
  
  tau_u[1] ~ normal(10, 10);
  tau_u[2] ~ normal(10, 10);
  
  tau_u[3] ~ normal(0, 3);
  tau_u[4] ~ normal(0, 3);
  tau_u[5] ~ normal(0, 3);
    
    
  param ~ multi_normal_cholesky(gm[1:N_centered], diag_pre_multiply(tau_u[1:N_centered], L_u)); 

  L_u ~ lkj_corr_cholesky(2);





  Y ~ binomial(npx, lapse + (1 - 2 * lapse) .* (0.5+0.5*erf(sqrt(beta / 2) .* (X-alpha))));
 
   
}

generated quantities{
  
  matrix[N_centered,N_centered] correlation_matrix;
  correlation_matrix = L_u*L_u';

  
  
  vector[N] log_lik;
  

  for (n in 1:N){
    log_lik[n] = binomial_lpmf(Y[n] | npx[n], lapse[n] + (1 - 2 * lapse[n]) .* (0.5+0.5*erf(sqrt(beta[n] / 2) .* (X[n]-alpha[n]))));
  
  }
}
