functions {
  // Our mean can never be less than zero because vapor pressure deficit cannot be less than zero.
  // But we need some small values so there is a gradient for solving the model, hence using the
  // exponential at negative values.
  real limit_positive(real x) {
    if (x >= 0) return x + 1.0 / 1000.0;
    return exp(x) / 1000.0;
  }
}
data {
  int<lower=0> N;        // Length of the input data.
  array[N] real year;    // The year, normalized to mean=0 and std=1
}
generated quantities {
  real a = gamma_rng(0.5, 1.0);     // The mean of the vapor pressure deficit data set should be 1
  real by = normal_rng(0.0, 0.1);   // The slope of the year term.
  real<lower=0> sigma = exponential_rng(1.0);

  vector[N] prior_pred;

  for (n in 1:N) {
    real mu = limit_positive(a + by * year[n]);
    real betap = mu / sigma^2;
    real alpha = mu * betap;
    prior_pred[n] = gamma_rng(alpha, betap);
  }
}
