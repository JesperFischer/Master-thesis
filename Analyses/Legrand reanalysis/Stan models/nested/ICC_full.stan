
data {
  int<lower=1> N_s1; // total number of observations (all subjects) (session1)
  int<lower=1> N_s2; // total number of observations (all subjects) (session2)
                      
  int<lower=1> S;  // total number of sessions (all subjects) (for now 2)
  int<lower=0> P;  //Total number of subjects

  vector[P*S] min_RT;


  array[S*P] int t_p_s_p_ses;
//  array[max(N_s1,N_s2)*S] int npx;  // total number of observations per X

  array[P*S, max(t_p_s_p_ses)] int Y;
  

  matrix[P*S , max(t_p_s_p_ses)] RT;

  matrix[P*S , max(t_p_s_p_ses)] Conf;

  
  matrix[P*S , max(t_p_s_p_ses)] X;  // design matrix (first column being intercept i.e. 1)

}

transformed data{
  int <lower=1> N_param;
  N_param = 10; // number of free parameters
  int idx;
}

parameters {
  // hierarchical group level means 
  vector [N_param] mu_g;
  // between subject variance
  vector <lower=0> [N_param] tau_b;
  // cholesky
  cholesky_factor_corr[N_param] L_b;
  // Difference for each subject:
  matrix[N_param, P] z_p;
  
  // within subject variance
  vector <lower=0> [N_param] tau_w;
  // cholesky
  cholesky_factor_corr[N_param] L_w;
  // Difference for each sessions:
  matrix[N_param, P*S] z_s;
  
}



transformed parameters{
  
  matrix[N_param, P] delta_mu_p;  //difference between global and participant
  matrix[N_param, P] mu_p;        //mean per participant
  matrix[N_param, P*S] mu_p_rep;  //mean per participant (for both session) (just append the mu_p)
  
  matrix[N_param, P*S] delta_fp_s; //difference per participant per session
  
  matrix[N_param, P*S] fp_s; //parameter value for the subject at the session in unbounded scale:
  

  matrix[P*S,max(N_s1,N_s2)] X_normal;
  matrix[P*S,max(N_s1,N_s2)] mu;
  matrix[P*S,max(N_s1,N_s2)] mu_rt;
  matrix[P*S,max(N_s1,N_s2)] mu_conf;
  
  
  vector[P*S] psyalpha;
  vector[P*S] psybeta;
  vector[P*S] lapse;
  
  vector[P*S] intercept;
  vector[P*S] beta;
  vector[P*S] sigma;
  vector[P*S] ndt;
  
  vector[P*S] intercept_conf;
  vector[P*S] beta_conf;
  vector[P*S] sigma_conf;
  
  
  
  delta_mu_p = diag_pre_multiply(tau_b, L_b) * z_p; //difference between global and participant
  
  delta_fp_s =  diag_pre_multiply(tau_w, L_w) * z_s; //difference within participant per session
  
  for(parm in 1:N_param){
    mu_p[parm,] = mu_g[parm] + delta_mu_p[parm,]; //adding the subject deviation to the global to get subject wise parameter estimates (mean):
  }
  
  mu_p_rep=append_col(mu_p,mu_p); // now we append it to itself such that its the right length and we can add the sessionwise deviation: (this should be done as many times as you have sessions)
  
  for(param in 1:N_param){
    fp_s[param,] = mu_p_rep[param,] + delta_fp_s[param,]; //adding session wise deviation
  }

  //getting subject estimates:
  psyalpha = to_vector(fp_s[1,]);
  psybeta = to_vector(exp(fp_s[2,]));
  lapse = to_vector(inv_logit(fp_s[3,]));
  
  intercept = to_vector(fp_s[4,]);
  beta = to_vector(exp(fp_s[5,]));
  sigma = to_vector(exp(fp_s[6,]));
  ndt = to_vector(inv_logit(fp_s[7,])) .* min_RT;
  
  intercept_conf = to_vector(fp_s[8,]);
  beta_conf = to_vector(fp_s[9,]);
  sigma_conf = to_vector(exp(fp_s[10,]));


  //model

  for(ses in 1:S){
    for(p in 1:P){

      X_normal[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] = 0.5 + 0.5 * erf((X[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] - psyalpha[(ses-1)*P+p]) / (psybeta[(ses-1)*P+p] * sqrt(2)));
      
      mu[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]]  = lapse[(ses-1)*P+p] + (1 - 2 * lapse[(ses-1)*P+p]) * X_normal[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]];

      mu_rt[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] = intercept[(ses-1)*P+p] + beta[(ses-1)*P+p] * mu[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] .* (1 - (mu[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]]));
      
      mu_conf[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] = intercept_conf[(ses-1)*P+p] + beta_conf[(ses-1)*P+p] * mu[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] .* (1 - (mu[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]]));
      
    }
    
  }
}
    



model{
  
  
  target += normal_lpdf(mu_g[1] | 0,20);
  target += normal_lpdf(mu_g[2] | 0,5);
  target += normal_lpdf(mu_g[3] | -4, 2);
  target += normal_lpdf(mu_g[4] | 0,5);
  target += normal_lpdf(mu_g[5] | 0,5);
  target += normal_lpdf(mu_g[6] | 0,5);
  target += normal_lpdf(mu_g[7] | 0,5);
  target += normal_lpdf(mu_g[8] | 0,5);
  target += normal_lpdf(mu_g[9] | 0,5);
  target += normal_lpdf(mu_g[10] | 0,5);
  
  target += std_normal_lpdf(to_vector(z_p));
  target += std_normal_lpdf(to_vector(z_s));
  
  target += normal_lpdf(tau_b[1] | 0, 20)-normal_lccdf(0 | 0, 20);
  target += normal_lpdf(tau_b[2:10] | 0, 5)-normal_lccdf(0 | 0, 5);

  target += normal_lpdf(tau_w[1] | 0, 10)-normal_lccdf(0 | 0, 10);
  target += normal_lpdf(tau_w[2:10] | 0, 5)-normal_lccdf(0 | 0, 5);
  
  target += lkj_corr_cholesky_lpdf(L_b | 2);
  target += lkj_corr_cholesky_lpdf(L_w | 2);

  
  

  for(ses in 1:S){
    for(p in 1:P){
     target += bernoulli_lpmf(Y[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] | mu[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]]);
     
     target += lognormal_lpdf(RT[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] - ndt[(ses-1)*P+p] |  mu_rt[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] , sigma[(ses-1)*P+p]);
     target += beta_proportion_lpdf(Conf[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]] |  inv_logit(mu_conf[(ses-1)*P+p,1:t_p_s_p_ses[(ses-1)*P+p]]) , sigma_conf[(ses-1)*P+p]);

    }
    
  }

}



generated quantities{

  matrix[N_param,N_param] between_correlation;
  matrix[N_param,N_param] within_correlation;
  

  vector[N_param] ICC;
  
  ICC = square(tau_b) ./ (square(tau_b) + square(tau_w));
  
  between_correlation = L_b*L_b';
  within_correlation = L_w*L_w';
  
    
}
