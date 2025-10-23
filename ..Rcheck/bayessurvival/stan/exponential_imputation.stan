// Exponential Survival Model for Bayesian Imputation
// Based on Moghaddam et al. (2022) methodology

data {
  // Data dimensions
  int<lower=0> n_obs;           // number of observed (uncensored) times
  int<lower=0> n_cens;          // number of censored times
  
  // Observed data
  vector<lower=0>[n_obs] t_obs;   // observed survival times
  vector<lower=0>[n_cens] t_cens; // censoring times
  
  // Prior parameters (passed as data for flexibility)
  // Rate parameter priors (gamma)
  real<lower=0> rate_prior_shape;
  real<lower=0> rate_prior_rate;
}

parameters {
  real<lower=1e-8> rate;        // Exponential rate parameter (lambda), strictly positive
}

model {
  // Prior
  rate ~ gamma(rate_prior_shape, rate_prior_rate);
  
  // Likelihood for observed times
  target += exponential_lpdf(t_obs | rate);
  
  // Likelihood for censored times (survival function)
  target += exponential_lccdf(t_cens | rate);
}

generated quantities {
  // One imputation per censored observation per MCMC iteration
  // This builds posterior distributions for each censored observation
  vector<lower=0>[n_cens] t_imputed;
  
  // Pointwise log-likelihoods for WAIC
  vector[n_obs] ll_obs;
  vector[n_cens] ll_cens;
  
  // Aggregate log-likelihood components for completeness
  real log_lik_obs = exponential_lpdf(t_obs | rate);
  real log_lik_cens = exponential_lccdf(t_cens | rate);
  real log_lik_total = log_lik_obs + log_lik_cens;
  
  // Fill pointwise log-likelihoods
  for (i in 1:n_obs) {
    ll_obs[i] = exponential_lpdf(t_obs[i] | rate);
  }
  for (i in 1:n_cens) {
    ll_cens[i] = exponential_lccdf(t_cens[i] | rate);
  }
  
  // Generate one imputation per censored observation
  for (i in 1:n_cens) {
    // Conditional sampling: t = t_cens[i] - log(U) / rate, where U ~ Uniform(0,1)
    real u = uniform_rng(0, 1);
    t_imputed[i] = t_cens[i] - log(u) / rate;
    // Numeric guard: ensure at least censoring time 
    if (t_imputed[i] < t_cens[i]) t_imputed[i] = t_cens[i] * (1.0 + 1e-12);
  }
} 
