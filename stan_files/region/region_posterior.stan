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
  int<lower=1> N;                       // Length of the input data.
  int<lower=1> K;                       // The number of sites in the data set.
  array[N] real year;                   // The year, normalized to mean=0 and std=1
  array[N] int<lower=0, upper=K> site;  // Integer code to ID site.
  array[N] real<lower=0> vpd;           // The scaled (proportional to mean) vapor pressure deficit.
}
parameters {
  real a_mu;
  real<lower=0> a_sigma;
  real by_mu;
  real<lower=0> by_sigma;
  real<lower=0> sigma_lam;

  array[K] real a;                      // per site intercept
  array[K] real by;                     // per site slope
  array[K] real<lower=0> sigma;         // per site variability
}
model {
  a_mu ~ normal(0.0, 0.5);          // The mean of the vapor pressure deficit data set should be 1
  a_sigma ~ exponential(1.0);
  by_mu ~ normal(0.0, 0.6);         // The slope of the year term.
  by_sigma ~ exponential(1.0);
  sigma_lam ~ exponential(1.0);

  {
    for (k in 1:K) {
      a[k] ~ normal(a_mu, a_sigma);
      by[k] ~ normal(by_mu, by_sigma);
      sigma[k] ~ exponential(sigma_lam);
    }
  }

  {
    vector[N] mu;
    vector[N] betap;
    vector[N] alpha;
    for (n in 1:N) {
      mu[n] = limit_positive(a[site[n]] + by[site[n]] * year[n]);
      betap[n] = mu[n] / sigma[site[n]]^2;
      alpha[n] = mu[n] * betap[n];
    }
    vpd ~ gamma(alpha, betap);
  }
    
}
generated quantities {

  vector[N]  post_pred;
  vector[N]  mu;

  for (n in 1:N) {
    mu[n] = limit_positive(a[site[n]] + by[site[n]] * year[n]);
    real betap = mu[n] / sigma[site[n]]^2;
    real alpha = mu[n] * betap;
    post_pred[n] = gamma_rng(alpha, betap);
  }
}
