data {
  int<lower=0> N;                   // Length of the input data, number of points.
  array[N] real year;               // The year, normalized to offset 1990 and scale 30 years.
}
generated quantities {
  // Hyper priors - make them broad enough to be learned from the data.
  real a_mu = normal_rng(1.0, 0.5);
  real<lower=0> a_sigma = exponential_rng(2.0);
  real by_mu = normal_rng(0.0, 0.3);
  real<lower=0> by_sigma = exponential_rng(2.0);
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
    mu[n] = a + by * year[n];
    prior_pred[n] = normal_rng(mu[n], exp(log_sigma));
  }
}
