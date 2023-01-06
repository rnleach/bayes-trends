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
  
  matrix cov_GPL2(matrix dist, matrix elev_diff, real sq_alpha, real sq_rho_d, real sq_rho_e, real delta) {
      int N = dims(dist)[1];
      matrix[N, N] K;
      for (i in 1:(N-1)) {
        K[i, i] = sq_alpha + delta;
        for (j in (i + 1):N) {
          K[i, j] = sq_alpha * exp(-(sq_rho_d * square(dist[i,j]) + sq_rho_e * square(elev_diff[i,j])));
          K[j, i] = K[i, j];
        }
      }
      K[N, N] = sq_alpha + delta;
      return K;
    }
}
data {
  int<lower=1> K;                       // The number of sites in the data set.
  int<lower=0> M;                       // The number of missing points.
  int<lower=1> N;                       // Length of the input data.
  array[N] int year;                    // The year.
  array[N] int<lower=0, upper=K> site;  // Integer code to ID site.
  array[N] real vpd;                    // The scaled (proportional to mean) vapor pressure deficit.

  matrix[K, K] dist_matrix;             // The distance between sites in km
  matrix[K, K] elev_matrix;             // The elevation difference between sites

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
  // real num_years= num_years_int;
}
parameters {
  array[K] real a;                      // per site intercept

  real a_mu;                            // all sites intercept mean 
  real<lower=0> a_sigma;                // all sites intercept variability

  real<lower=0> etasq;                  // Maximum covariance from Gaussian Mixture
  real<lower=0> rhosq_d;                // Length (distance) scale for Gaussian
  real<lower=0> rhosq_e;                // Length (elevation differnce) scale for Gaussian

  vector[M] missing_vals;               // The missing values to be imputed.

}
transformed parameters {
  matrix[K, K] L_SIGMA;
  matrix[K, K] SIGMA;

  SIGMA = cov_GPL2(dist_matrix, elev_matrix, etasq, rhosq_d, rhosq_e, 0.01);
  L_SIGMA = cholesky_decompose(SIGMA);
}
model {

  // Set the hyper-priors for the linear relationship.
  a_mu ~ normal(1.0, 0.5);
  a_sigma ~ exponential(2.0);

  // Set the hyper-priors for the covariance matrix.
  etasq ~ exponential(1);
  rhosq_d ~ exponential(100);
  rhosq_e ~ exponential(10);

  // Set the priors for each site and model the correlation between them.
  a ~ normal(a_mu, a_sigma);

  // Set the priors for the missing data values.
  for( m in 1:M ) {
    missing_vals[m] ~ normal(miss_means[m], miss_stds[m]);
  }

  // Merge the missing data in and rearrange the data into a years x site matrix.
  {
    vector[N] miss_merge = merge_missing(miss_indexes, to_vector(vpd), missing_vals);

    matrix[NY, K] vpd_by_year;
    matrix[NY, K] mu_by_year;

    for( n in 1:N) {
      int y = year_to_index(year[n], min_year, num_years_int);
      int s = site[n];

      vpd_by_year[y, s] = miss_merge[n];
      mu_by_year[y, s] = a[s];
    }

    for( y in 1:NY) {
      vpd_by_year[y] ~ multi_normal_cholesky(mu_by_year[y], L_SIGMA);
    }
  }
}
generated quantities {
  vector[N] post_pred;
  vector[N] mu;

  {
    matrix[NY, K] post_pred_matrix;
    matrix[NY, K] post_pred_mu_matrix;

    for( n in 1:N) {
      int k = site[n];
      int y = year_to_index(year[n], min_year, num_years_int);

      post_pred_mu_matrix[y, k] = a[k];
    }

    for( y in 1:NY) {
      post_pred_matrix[y] = multi_normal_cholesky_rng(
          post_pred_mu_matrix[y], L_SIGMA)';
    }

    for( n in 1:N) {
      int k = site[n];
      int y = year_to_index(year[n], min_year, num_years_int);
      
      mu[n] = post_pred_mu_matrix[y, k];
      post_pred[n] = post_pred_matrix[y, k];
    }

  }
}

/*
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data{
    int T[10];
    int society[10];
    int P[10];
    matrix[10,10] Dmat;
}
parameters{
    vector[10] z;
    real<lower=0> g;
    real<lower=0> b;
    real<lower=0> a;
    real<lower=0> etasq;
    real<lower=0> rhosq;
}
transformed parameters{
    vector[10] k;
    matrix[10,10] L_SIGMA;
    matrix[10,10] SIGMA;
    SIGMA = cov_GPL2(Dmat, etasq, rhosq, 0.01);
    L_SIGMA = cholesky_decompose(SIGMA);
    k = L_SIGMA * z;
}
model{
    vector[10] lambda;
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 2 );
    a ~ exponential( 1 );
    b ~ exponential( 1 );
    g ~ exponential( 1 );
    z ~ normal( 0 , 1 );
    for ( i in 1:10 ) {
        lambda[i] = (a * P[i]^b/g) * exp(k[society[i]]);
    }
    T ~ poisson( lambda );
}
*/
 
