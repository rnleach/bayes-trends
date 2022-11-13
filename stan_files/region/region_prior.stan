functions {
  // The link function is not a traditional exponential because that would distort the meaning of 
  // the linear coefficients. However, the mean still MUST be positive for the Gamma distribution.
  // So the link function (practically) preserves the value of the linear model when it predicts a
  // mean value above zero. Otherwise it keeps the value very small, but still greater than zero.
  // This function effectively truncates predicted mean values at zero. So the linear parameters can
  // still be interpreted directly (with some care).
  real limit_positive(real x) {
    if (x >= 0) return x + 1.0 / 1000.0;
    return exp(x) / 1000.0;
  }
}
data {
  int<lower=0> N;                   // Length of the input data, number of points.
  array[N] real year;               // The year, normalized to mean=0 and std=1
}
generated quantities {
  // Hyper priors - make them broad enough to be learned from the data.
  real a_mu = normal_rng(1.0, 0.5);
  real<lower=0> a_sigma = exponential_rng(1.0);
  real by_mu = normal_rng(0.0, 0.6);
  real<lower=0> by_sigma = exponential_rng(1.0);
  real<lower=0> sigma_lam = exponential_rng(1.0);

  // Use the hyper-priors to determine the distribution
  real a = normal_rng(a_mu, a_sigma);
  real by = normal_rng(by_mu, by_sigma);   
  real<lower=0> sigma = exponential_rng(sigma_lam);

  // Prior predictions.
  vector[N] prior_pred;

  // Save the simulated means.
  vector[N] mu;

  for (n in 1:N) {
    mu[n] = limit_positive(a + by * year[n]);
    real betap = mu[n] / sigma^2;
    real alpha = mu[n] * betap;
    prior_pred[n] = gamma_rng(alpha, betap);
  }
}
