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
  array[N] real<lower=0> vpd;     // The scaled (proportional to mean) vapor pressure deficit.
}
parameters {
  real<lower=0> a;
  real<lower=0> by;
  real<lower=0> sigma;
}
model {
  a ~ gamma(0.5, 1.0);              // The mean of the vapor pressure deficit data set should be 1
  by ~ normal(0.0, 0.1);            // The slope of the year term.
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
