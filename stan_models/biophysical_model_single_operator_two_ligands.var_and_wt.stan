// Fit dose-response curves to Phillips lab model for allosteric TFs
//     this version of the model calculates shifts of the biophysical parameters relative to the WT
//     It uses data for the variant and for the WT
//     It is meant to be used with GROQ-Seq data with 2 ligands

data {
  int<lower=1> N_lig;         // number of non-zero ligand concentrations for each ligand
  int<lower=1> N_zero;        // number of replicate measurements at zero ligand
  
  vector[N_lig] x1;           // non-zero ligand 1 concentrations
  vector[N_lig] x2;           // non-zero ligand 2 concentrations
  
  vector[N_zero] log_y0;      // log10 of gene expression at zero ligand
  vector[N_zero] log_y0err;   // estimated error
  
  vector[N_lig] log_y1;       // log10 of gene expression at each concentration for ligand 1
  vector[N_lig] log_y1err;    // estimated error of log10 gene expression at each concentration
  
  vector[N_lig] log_y2;       // log10 of gene expression at each concentration for ligand 1
  vector[N_lig] log_y2err;    // estimated error of log10 gene expression at each concentration
  
  // Matching data for wild-type
  int<lower=1> N_lig_wt;         // number of non-zero ligand concentrations for each ligand
  int<lower=1> N_zero_wt;        // number of replicate measurements at zero ligand
  
  vector[N_lig_wt] x1_wt;           // non-zero ligand 1 concentrations
  vector[N_lig_wt] x2_wt;           // non-zero ligand 2 concentrations
  
  vector[N_zero_wt] log_y0_wt;      // log10 of gene expression at zero ligand
  vector[N_zero_wt] log_y0err_wt;   // estimated error
  
  vector[N_lig_wt] log_y1_wt;       // log10 of gene expression at each concentration for ligand 1
  vector[N_lig_wt] log_y1err_wt;    // estimated error of log10 gene expression at each concentration
  
  vector[N_lig_wt] log_y2_wt;       // log10 of gene expression at each concentration for ligand 1
  vector[N_lig_wt] log_y2err_wt;    // estimated error of log10 gene expression at each concentration
  
  // y_max is the same for variant and wild-type
  real y_max;                 // geometric mean for prior on maximum gene expression value
  real y_max_prior_width;     // geometric std for prior on maximum gene expression value

  // Input data for priors on wild-type free energy parameters (single- and multi-operator models)
  // ligand 1 interactions
  real log_k_a_1_wt_prior_mean;
  real log_k_a_1_wt_prior_std;
  real log_k_i_1_wt_prior_mean;
  real log_k_i_1_wt_prior_std;
  // ligand 2 interactions
  real log_k_a_2_wt_prior_mean;
  real log_k_a_2_wt_prior_std;
  real log_k_i_2_wt_prior_mean;
  real log_k_i_2_wt_prior_std;
  
  real delta_eps_AI_wt_prior_mean;
  real delta_eps_AI_wt_prior_std;
  real delta_eps_RA_wt_prior_mean;
  real delta_eps_RA_wt_prior_std;
  
  // priors on mutational effects and epistasis  (single- and multi-operator models)
  real delta_prior_width;   // width of prior on "_mut" parameters (mutational effects)
  real eps_RA_prior_scale;  // scale factor for delta_eps_RA_mut - so that it can be made smaller for mutations far from the DNA
  
}

transformed data {

  // transformed data variable declarations shared by all models
  real ln_10;
  real hill_n;
  real N_NS; // number of non-sepcific binding sites (length of genome)
  real R;   // copy number of TF per cell
  
  hill_n = 2;
    
  ln_10 = log(10.0);
  
  R = 200;
  N_NS = 4600000;
  
}

parameters {
  real log_k_a_1_wt;         // log10 of IPTG binding affinity to active state
  real log_k_a_1_mut;        // shift relative to WT due to mutations
  
  real log_k_i_1_wt;         // log10 of IPTG binding affinity to inactive state
  real log_k_i_1_mut;        // shift relative to WT due to mutations
  
  real log_k_a_2_wt;         // log10 of IPTG binding affinity to active state
  real log_k_a_2_mut;        // shift relative to WT due to mutations
  
  real log_k_i_2_wt;         // log10 of IPTG binding affinity to inactive state
  real log_k_i_2_mut;        // shift relative to WT due to mutations
  
  real delta_eps_AI_wt;    // free energy difference between active and inactive states
  real delta_eps_AI_mut;   // shift relative to WT due to mutations
  
  real delta_eps_RA_wt;    // free energy for Active TF to operator
  real delta_eps_RA_mut;   // shift relative to WT due to mutations
  
  real log_y_max;          // log10 of maximum possible gene expression
  
  real<lower=0> sigma;      // scale factor for standard deviation of noise in log_y
  
}

transformed parameters {
  real K_A_1;
  real K_I_1;
  real log_k_a_1_var;
  real log_k_i_1_var;
  real K_A_2;
  real K_I_2;
  real log_k_a_2_var;
  real log_k_i_2_var;
  
  real K_A_1_wt;
  real K_I_1_wt;
  real K_A_2_wt;
  real K_I_2_wt;
  
  real delta_eps_AI_var;
  real delta_eps_RA_var;
  
  real log_mean_y0;
  vector[N_lig] log_mean_y1;
  vector[N_lig] log_mean_y2;
  
  real log_mean_y0_wt;
  vector[N_lig_wt] log_mean_y1_wt;
  vector[N_lig_wt] log_mean_y2_wt;
  
  log_k_a_1_var = log_k_a_1_wt + log_k_a_1_mut;
  log_k_i_1_var = log_k_i_1_wt + log_k_i_1_mut;
  log_k_a_2_var = log_k_a_2_wt + log_k_a_2_mut;
  log_k_i_2_var = log_k_i_2_wt + log_k_i_2_mut;
  
  delta_eps_AI_var = delta_eps_AI_wt + delta_eps_AI_mut;
  delta_eps_RA_var = delta_eps_RA_wt + delta_eps_RA_mut;
  
  K_A_1_wt = 10^log_k_a_1_wt;
  K_I_1_wt = 10^log_k_i_1_wt;
  K_A_2_wt = 10^log_k_a_2_wt;
  K_I_2_wt = 10^log_k_i_2_wt;
  
  K_A_1 = 10^log_k_a_1_var;
  K_I_1 = 10^log_k_i_1_var;
  K_A_2 = 10^log_k_a_2_var;
  K_I_2 = 10^log_k_i_2_var;
  
  {
	real fold_change;
    real c1;
    real c2;
    real c3;
    
    // variant
    c3 = R/N_NS * exp(-delta_eps_RA_var);
    c1 = 1;
    c2 = exp(-delta_eps_AI_var);
	
    fold_change = 1/(1 + (c1/(c1+c2))*c3);
	
    log_mean_y0 = log_y_max + log10(fold_change);
    
    // wild-type
    c3 = R/N_NS * exp(-delta_eps_RA_wt);
    c1 = 1;
    c2 = exp(-delta_eps_AI_wt);
	
    fold_change = 1/(1 + (c1/(c1+c2))*c3);
	
    log_mean_y0_wt = log_y_max + log10(fold_change);
  }
  
  //variant
  for (i in 1:N_lig) {
	real fold_change;
    real c1;
    real c2;
    real c3;
	
    c3 = R/N_NS * exp(-delta_eps_RA_var);
    c1 = (1 + x1[i]/K_A_1)^hill_n;
    c2 = ( (1 + x1[i]/K_I_1)^hill_n ) * exp(-delta_eps_AI_var);
	
    fold_change = 1/(1 + (c1/(c1+c2))*c3);
	
    log_mean_y1[i] = log_y_max + log10(fold_change);
	
    c1 = (1 + x2[i]/K_A_2)^hill_n;
    c2 = ( (1 + x2[i]/K_I_2)^hill_n ) * exp(-delta_eps_AI_var);
	
    fold_change = 1/(1 + (c1/(c1+c2))*c3);
	
    log_mean_y2[i] = log_y_max + log10(fold_change);
  }
	
  //wild-type
  for (i in 1:N_lig_wt) {
	real fold_change;
    real c1;
    real c2;
    real c3;
	
    c3 = R/N_NS * exp(-delta_eps_RA_wt);
    c1 = (1 + x1_wt[i]/K_A_1_wt)^hill_n;
    c2 = ( (1 + x1_wt[i]/K_I_1_wt)^hill_n ) * exp(-delta_eps_AI_wt);
	
    fold_change = 1/(1 + (c1/(c1+c2))*c3);
	
    log_mean_y1_wt[i] = log_y_max + log10(fold_change);
	
    c1 = (1 + x2_wt[i]/K_A_2_wt)^hill_n;
    c2 = ( (1 + x2_wt[i]/K_I_2_wt)^hill_n ) * exp(-delta_eps_AI_wt);
	
    fold_change = 1/(1 + (c1/(c1+c2))*c3);
	
    log_mean_y2_wt[i] = log_y_max + log10(fold_change);
  }
  
}

model {
  // priors on free energy params
  log_k_a_1_wt ~ normal(log_k_a_1_wt_prior_mean, log_k_a_1_wt_prior_std);
  log_k_i_1_wt ~ normal(log_k_i_1_wt_prior_mean, log_k_i_1_wt_prior_std);
  log_k_a_2_wt ~ normal(log_k_a_2_wt_prior_mean, log_k_a_2_wt_prior_std);
  log_k_i_2_wt ~ normal(log_k_i_2_wt_prior_mean, log_k_i_2_wt_prior_std);
  
  delta_eps_AI_wt ~ normal(delta_eps_AI_wt_prior_mean, delta_eps_AI_wt_prior_std);
  delta_eps_RA_wt ~ normal(delta_eps_RA_wt_prior_mean, delta_eps_RA_wt_prior_std);
  
  log_k_a_1_mut ~ normal(0, delta_prior_width/ln_10); // factor of 1/ln_10 is to compensate for use of log10 instead of ln
  log_k_i_1_mut ~ normal(0, delta_prior_width/ln_10); // factor of 1/ln_10 is to compensate for use of log10 instead of ln
  log_k_a_2_mut ~ normal(0, delta_prior_width/ln_10); // factor of 1/ln_10 is to compensate for use of log10 instead of ln
  log_k_i_2_mut ~ normal(0, delta_prior_width/ln_10); // factor of 1/ln_10 is to compensate for use of log10 instead of ln
  
  delta_eps_AI_mut ~ normal(0, delta_prior_width);
  delta_eps_RA_mut ~ normal(0, delta_prior_width*eps_RA_prior_scale);
  
  // prior on max output level
  log_y_max ~ normal(log10(y_max), y_max_prior_width);
  
  // prior on scale parameter for log-normal measurement error
  sigma ~ normal(0, 1);
  
  // model of the data (dose-response curve with noise)
  log_y0 ~ normal(log_mean_y0, sigma*log_y0err);
  log_y1 ~ normal(log_mean_y1, sigma*log_y1err);
  log_y2 ~ normal(log_mean_y2, sigma*log_y2err);
  log_y0_wt ~ normal(log_mean_y0_wt, sigma*log_y0err_wt);
  log_y1_wt ~ normal(log_mean_y1_wt, sigma*log_y1err_wt);
  log_y2_wt ~ normal(log_mean_y2_wt, sigma*log_y2err_wt);

}

generated quantities {
  real g_max;
  vector[19] x_out;
  
  vector[19] y1_out;
  vector[19] y2_out;
  vector[19] fc1_out;
  vector[19] fc2_out;
  vector[N_lig] dev_log_y1;       // log_mean_y - log_y
  vector[N_lig] dev_log_y2;       // log_mean_y - log_y
  
  real log_k_ratio_1_mut;        // shift relative to WT due to mutations for k_a : k_i ratio
  real log_k_ratio_2_mut;        // shift relative to WT due to mutations for k_a : k_i ratio
  
  log_k_ratio_1_mut = log_k_a_1_mut - log_k_i_1_mut;
  log_k_ratio_2_mut = log_k_a_2_mut - log_k_i_2_mut;
  
  dev_log_y1 = log_mean_y1 - log_y1;
  dev_log_y2 = log_mean_y2 - log_y2;
  
  g_max = 10^log_y_max;
  
  x_out[1] = 0;
  for (i in 2:19) {
    x_out[i] = 2^(i-2);
  }
  
  for (i in 1:19) {
    real c1;
    real c2;
    real c3;
	
    c3 = R/N_NS * exp(-delta_eps_RA_var);
    c1 = (1 + x_out[i]/K_A_1)^hill_n;
    c2 = ( (1 + x_out[i]/K_I_1)^hill_n ) * exp(-delta_eps_AI_var);
	
    y1_out[i] = g_max/(1 + (c1/(c1+c2))*c3);
    fc1_out[i] = 1/(1 + (c1/(c1+c2))*c3);
	
    c1 = (1 + x_out[i]/K_A_2)^hill_n;
    c2 = ( (1 + x_out[i]/K_I_2)^hill_n ) * exp(-delta_eps_AI_var);
	
    y2_out[i] = g_max/(1 + (c1/(c1+c2))*c3);
    fc2_out[i] = 1/(1 + (c1/(c1+c2))*c3);

  }
  
}
