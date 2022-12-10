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
  int<lower=1> N;                       // Length of the input data.
  int<lower=1> K;                       // The number of sites in the data set.
  array[N] real year;                   // The year, normalized to mean=0 and std=1
  array[N] int<lower=0, upper=K> site;  // Integer code to ID site.
  array[N] real<lower=0> vpd;           // The scaled (proportional to mean) vapor pressure deficit.
}
parameters {
  array[K] real a;                      // per site intercept
  array[K] real by;                     // per site slope
  array[K] real log_sigma;              // per site variability

  real a_mu;                            // all sites intercept mean 
  real<lower=0> a_sigma;                // all sites intercept variability
  real by_mu;                           // all sites slope mean
  real<lower=0> by_sigma;               // all sites slope variability
  real log_sigma_mu;                    // all sites variability mean
  real <lower=0> log_sigma_sigma;       // all sites variability variability

  corr_matrix[3] rho;                   // Correlation matrix for intercept, slope, and variability
}
model {
  a_mu ~ normal(1.0, 0.5);
  a_sigma ~ exponential(1.0);
  by_mu ~ normal(0.0, 0.6);
  by_sigma ~ exponential(1.0);
  log_sigma_mu ~ normal(0.0, 0.5);
  log_sigma_sigma ~ exponential(1.0);
  rho ~ lkj_corr(2.0);

  {
    vector[3] MU;
    MU = [a_mu, by_mu, log_sigma_mu]';

    vector[3] SIGMA;
    SIGMA = [a_sigma, by_sigma, log_sigma_sigma]';

    array[K] vector[3] YY;
    for( k in 1:K) YY[k] = [a[k] , by[k], log_sigma[k]]';
    YY ~ multi_normal(MU , quad_form_diag(rho, SIGMA));
  }

  {
    vector[N] mu;
    vector[N] betap;
    vector[N] alpha;
    for (n in 1:N) {
      mu[n] = limit_positive(a[site[n]] + by[site[n]] * year[n]);
      betap[n] = mu[n] / exp(log_sigma[site[n]])^2;
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
    real betap = mu[n] / exp(log_sigma[site[n]])^2;
    real alpha = mu[n] * betap;
    post_pred[n] = gamma_rng(alpha, betap);
  }
}
