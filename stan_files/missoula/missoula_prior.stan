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
  // The mean of the vapor pressure should be 1, so the mean the intercept should be near there too
  // absent a steep slope. A wide standard deviation means that the intercept can be adjusted to 
  // allow for a steep slope.
  real a = normal_rng(1.0, 1.0);

  // The slope of the year term. Mean should be zero to assume no preference for a moistening or
  // drying trend. Assume the trend is slow, so the slope should be pretty small.
  real by = normal_rng(0.0, 0.1);   

  // Prefer small standard deviations, but the tails of the exponential are "fat" enough to not
  // effectively rule out large standard deviations.
  real<lower=0> sigma = exponential_rng(1.0);

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
