// Log-Normal Survival Model for Bayesian Imputation
// Based on Moghaddam et al. (2022) methodology

data {
  // Data dimensions
  int<lower=0> n_obs;           // number of observed (uncensored) times
  int<lower=0> n_cens;          // number of censored times
  
  // Observed data
  vector<lower=0>[n_obs] t_obs;   // observed survival times
  vector<lower=0>[n_cens] t_cens; // censoring times
  
  // Prior parameters (passed as data for flexibility)
  // Location parameter (mean of log-times) - normal prior
  real mu_prior_mean;
  real<lower=0> mu_prior_sd;
  
  // Scale parameter (sd of log-times) - half-normal prior
  real<lower=0> sigma_prior_sd;
}

parameters {
  real mu;                      // Log-normal location parameter (mean of log-times)
  real<lower=1e-6> sigma;      // Log-normal scale parameter (sd of log-times)
}

model {
  // Priors
  mu ~ normal(mu_prior_mean, mu_prior_sd);
  sigma ~ normal(0, sigma_prior_sd);  // Half-normal (sigma constrained > 0)
  
  // Likelihood for observed times
  target += lognormal_lpdf(t_obs | mu, sigma);
  
  // Likelihood for censored times (survival function)
  target += lognormal_lccdf(t_cens | mu, sigma);
}

generated quantities {
  // One imputation per censored observation per MCMC iteration
  // This builds posterior distributions for each censored observation
  vector<lower=0>[n_cens] t_imputed;
  
  // Pointwise log-likelihoods for WAIC
  vector[n_obs] ll_obs;
  vector[n_cens] ll_cens;
  
  // Aggregate log-likelihood components
  real log_lik_obs = lognormal_lpdf(t_obs | mu, sigma);
  real log_lik_cens = lognormal_lccdf(t_cens | mu, sigma);
  real log_lik_total = log_lik_obs + log_lik_cens;
  
  // Fill pointwise log-likelihoods
  for (i in 1:n_obs) {
    ll_obs[i] = lognormal_lpdf(t_obs[i] | mu, sigma);
  }
  for (i in 1:n_cens) {
    ll_cens[i] = lognormal_lccdf(t_cens[i] | mu, sigma);
  }
  
  // Generate one imputation per censored observation
  for (i in 1:n_cens) {
    // Sample from conditional distribution f(t|T >= t_cens[i], mu, sigma)
    // For log-normal, we work in log-space then exponentiate
    real u = uniform_rng(0, 1);
    real log_t_cens = log(t_cens[i]);
    
    // CDF at censoring time, with numerical safeguards
    real Phi_c = normal_cdf(log_t_cens | mu, sigma);
    if (Phi_c < 1e-10) Phi_c = 1e-10;
    if (Phi_c > 0.999) Phi_c = 0.999;
    
    // Transform to truncated space and compute inverse CDF stably
    real p = Phi_c + u * (1 - Phi_c);
    if (p <= 0) p = 1e-10;
    if (p >= 1) p = 0.999999;
    
    real z = inv_Phi(p);
    real log_t_imputed = mu + sigma * z;
    t_imputed[i] = exp(log_t_imputed);
    
    // Numeric guard: ensure at least censoring time 
    if (t_imputed[i] < t_cens[i]) t_imputed[i] = t_cens[i] * (1.0 + 1e-12);
  }
} 
