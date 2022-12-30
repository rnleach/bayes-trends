functions {

  vector merge_missing( array[] int miss_indexes, vector x_obs, vector x_miss) {
    int N = dims(x_obs)[1];
    int N_miss = dims(x_miss)[1];

    vector[N] merged;
    merged = x_obs;
    for ( i in 1:N_miss ) {
      // Requires a plus 1 assuming we come from a system with 0 based indexs instead of 1 based.
      merged[ miss_indexes[i] + 1 ] = x_miss[i];
    }

    return merged;
  }

  int year_to_index(int year, int min_year, int num_years) {
    int index = year - min_year + 1;
    if( index > num_years) {
      reject("Invalid year index greater than number of years!");
    }
    return index;
  }
  
}
data {
  int<lower=1> K;                       // The number of sites in the data set.
  int<lower=0> M;                       // The number of missing points.
  int<lower=1> N;                       // Length of the input data.
  array[N] int year;                    // The year.
  array[N] int<lower=0, upper=K> site;  // Integer code to ID site.
  array[N] real vpd;                    // The scaled (proportional to mean) vapor pressure deficit.

  array[M] int miss_indexes;            // The indexes of the missing values.
  array[M] real miss_means;             // The mean value to use for each missing value prior.
  array[M] real miss_stds;              // The standard deviation to use for each missing value prior.
}
transformed data {
  int<lower=2> NY = N %/% K;            // The number of years in the data. This requires that years
                                        // with missing data hav NaN values and are properly 
                                        // encoded in the 'miss_*' data values.

  int<lower=1900, upper=2023> min_year = min(year); 
  int<lower=1900, upper=2023> max_year = max(year);
  int num_years_int = max_year - min_year + 1;
  real num_years= num_years_int;
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

  vector[M] missing_vals;               // The missing values to be imputed.
}
model {
  // Set the hyper-priors for the linear relationship.
  a_mu ~ normal(1.0, 0.5);
  a_sigma ~ exponential(2.0);
  by_mu ~ normal(0.0, 0.3);
  by_sigma ~ exponential(2.0);
  log_sigma_mu ~ normal(0.0, 0.5);
  log_sigma_sigma ~ exponential(1.0);
  rho ~ lkj_corr(2.0);

  // Set the priors for each site and model the correlation between them.
  {
    vector[3] MU;
    MU = [a_mu, by_mu, log_sigma_mu]';

    vector[3] SIGMA;
    SIGMA = [a_sigma, by_sigma, log_sigma_sigma]';

    array[K] vector[3] YY;
    for( k in 1:K) YY[k] = [a[k] , by[k], log_sigma[k]]';
    YY ~ multi_normal(MU , quad_form_diag(rho, SIGMA));
  }

  // Set the priors for the missing data values.
  for( m in 1:M ) {
    missing_vals[m] ~ normal(miss_means[m], miss_stds[m]);
  }

  // Merge the missing data in and model the site level linear relationships.
  {
    vector[N] miss_merge = merge_missing(miss_indexes, to_vector(vpd), missing_vals);

    matrix[NY, K] vpd_by_year;
    matrix[NY, K] mu_by_year;
    matrix[NY, K] sigma_by_year;

    for( n in 1:N) {
      int y = year_to_index(year[n], min_year, num_years_int);
      int k = site[n];

      vpd_by_year[y, k] = miss_merge[n];
      mu_by_year[y, k] = a[k] + by[k] * (year[n] - min_year) / num_years;
      sigma_by_year[y, k] = exp(log_sigma[k]);
    }

    for( y in 1:NY) {
      vpd_by_year[y] ~ normal(mu_by_year[y], sigma_by_year[y]);
    }
  }
}
generated quantities {

  vector[N]  post_pred;
  vector[N]  mu;

  for (n in 1:N) {
    mu[n] = a[site[n]] + by[site[n]] * (year[n] - min_year) / num_years;
    real sigma = exp(log_sigma[site[n]]);
    post_pred[n] = normal_rng(mu[n], sigma);
  }
}
