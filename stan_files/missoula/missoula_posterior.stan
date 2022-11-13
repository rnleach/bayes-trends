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
  array[N] real<lower=0> vpd;       // The scaled (proportional to mean) vapor pressure deficit.
}
parameters {
  real a;
  real by;
  real<lower=0> sigma;
}
model {
  a ~ normal(1.0, 1.0);
  by ~ normal(0.0, 0.1);
  sigma ~ exponential(1.0);

  vector[N]  mu;
  vector[N]  betap;
  vector[N]  alpha;

  for (n in 1:N) {
    mu[n] = limit_positive(a + by * year[n]);
    betap[n] = mu[n] / sigma^2;
    alpha[n] = mu[n] * betap[n];
  }
  
  vpd ~ gamma(alpha, betap);

}
generated quantities {

  vector[N]  post_pred;
  vector[N]  mu;

  for (n in 1:N) {
    mu[n] = limit_positive(a + by * year[n]);
    real betap = mu[n] / sigma^2;
    real alpha = mu[n] * betap;
    post_pred[n] = gamma_rng(alpha, betap);
  }
}
