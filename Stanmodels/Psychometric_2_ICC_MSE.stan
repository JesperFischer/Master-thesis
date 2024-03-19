data {
  int<lower=1> T;
  int<lower=1> S; //n sessions
  int<lower=1> P; //n participants

  array[S*P,T] int resp;
  matrix[S*P,T] X;
  
  vector[S*P] real_alpha;
  vector[S*P] real_beta;
  vector[S*P] real_lambda;
}
transformed data{
  int N = 3;
  
}

parameters {

    // Group means (of parameters)
  vector[N] mu_g;
  // Between participant scales
  vector<lower = 0>[N]  tau_b;
  // Between participant cholesky decomposition
  cholesky_factor_corr[N] L_b;
  
  // Participant deviation (from the group level) 
  matrix[N, P] z_p;  
  
  // Within participant scales
  vector<lower = 0>[N]  tau_w;
  // Within participant cholesky decomposition
  cholesky_factor_corr[N] L_w;
  // Session deviation (from own group level)
  matrix[N, S*P] z_s; 
  
}

transformed parameters{
  
  
  matrix[N, P] delta_mu_p;  
  matrix[N, P] mu_p;  
  matrix[N, S*P] mu_p_rep;  
  matrix[N, S*P] fp_s;


  vector[S*P] lambda;
  vector[S*P] alpha;
  vector[S*P] beta;  

  


  ///Recomposition (deviations from the group levels)
  
  delta_mu_p = diag_pre_multiply(tau_b, L_b) * z_p;

  
  for(idx in 1:P){
    mu_p[,idx] = mu_g + delta_mu_p[,idx];
  }
  
  
  // adding another as there are 2 sesssions
  
  mu_p_rep=append_col(mu_p,mu_p);
  
  // deviations from individual level


  fp_s = mu_p_rep + diag_pre_multiply(tau_w, L_w) * z_s;
  
  

  
  beta = to_vector(fp_s[1,]);
  lambda = to_vector(fp_s[2,]);
  alpha = to_vector(fp_s[3,]);

}


model {
  
  int p_s;
  mu_g[1] ~ normal(0,3);
  mu_g[2] ~ normal(-4,2);
  mu_g[3] ~ normal(0,10);
  
  //Between participant scales
  tau_b[1] ~ normal(0,3);
  tau_b[2] ~ normal(0,3);
  tau_b[3] ~ normal(0,10);
  
  //Between participant cholesky decomposition
  L_b ~ lkj_corr_cholesky(2);
  
  // Participant deviation 
  to_vector(z_p) ~ std_normal();  
  
  // Within participant scales
  tau_w[1] ~ normal(0,3);
  tau_w[2] ~ normal(0,3);
  tau_w[3] ~ normal(0,10);

  // Within participant cholesky decomposition
  L_w ~ lkj_corr_cholesky(2);
  
  // Session deviation 
  to_vector(z_s) ~ std_normal();  
  
  for(s in 1:S){
    for(p in 1:P){
          p_s = p;
          if(s == 2){
            p_s = p+P;
          }
            
            target += bernoulli_lpmf(resp[p_s,] | inv_logit(lambda[p_s]) / 2 + (1 - 2 * inv_logit(lambda[p_s]) / 2) * (0.5 + 0.5 * erf((X[p_s,] - alpha[p_s]) / (exp(beta[p_s]) * sqrt(2)))));
      
        }
    }
}

generated quantities{
  
  matrix[N,N] Corr_B;
  matrix[N,N] Corr_W;

  real resid_alpha_var = mean(square(real_alpha-alpha));
  real resid_lambda_var = mean(square(real_lambda-lambda));
  real resid_beta_var = mean(square(real_beta-beta));

  vector[N] ICC = square(tau_b) ./ (square(tau_b) + square(tau_w));
  
  real ICC_alpha = square(tau_b[3]) ./ (square(tau_b[3]) + square(tau_w[3]) + resid_alpha_var);
  
  real ICC_lambda = square(tau_b[2]) ./ (square(tau_b[2]) + square(tau_w[2]) + resid_lambda_var);
  
  real ICC_beta = square(tau_b[1]) ./ (square(tau_b[1]) + square(tau_w[1]) + resid_beta_var);  

  Corr_B = L_b * L_b';
  Corr_W = L_w * L_w';

  
}
