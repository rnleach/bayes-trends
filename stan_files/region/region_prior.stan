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
  real log_sigma_mu = normal_rng(0.0, 0.5);
  real<lower=0> log_sigma_sigma = exponential_rng(1.0);
  corr_matrix[3] rho = lkj_corr_rng(3, 2.0);

  // Use the hyper-priors to determine the distribution
  real a;
  real by;
  real log_sigma;
  {
    vector[3] MU = [a_mu, by_mu, log_sigma_mu]';
    vector[3] SIGMA = [a_sigma, by_sigma, log_sigma_sigma]';
    vector[3] params = multi_normal_rng(MU, quad_form_diag(rho, SIGMA));
    a = params[1];
    by = params[2];
    log_sigma = params[3];
  }

  // Prior predictions.
  vector[N] prior_pred;

  // Save the simulated means.
  vector[N] mu;

  for (n in 1:N) {
    mu[n] = limit_positive(a + by * year[n]);
    real betap = mu[n] / exp(log_sigma)^2;
    real alpha = mu[n] * betap;
    prior_pred[n] = gamma_rng(alpha, betap);
  }
}
