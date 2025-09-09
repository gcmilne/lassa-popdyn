// ** Model of Mastomys natalensis body weight change (growth) ** //

data {
  // Number of observations
  int<lower=0> N;
  // Continuous response variable
  vector[N] y;
  // No. predictors for mu
  int<lower=0> K_mu;
  // No. predictors for log(sigma)
  int<lower=0> K_sigma;
  // Design matrix for mu (for weight & pop)
  matrix[N, K_mu-6] X;
  // Design matrix for log(sigma) (for weight)
  matrix[N, K_sigma] Z;
  // Total no. temperature-precipitation lag combinations
  int<lower=1> L;
  // Matrices for decomposed climate values at each lag combination
  matrix[N, L] t_seas_mat;    
  matrix[N, L] t_var_mat;
  matrix[N, L] t_trend_mat;
  matrix[N, L] p_seas_mat;
  matrix[N, L] p_var_mat;
  matrix[N, L] p_trend_mat;
  // Binary switches for predictors (1 = include, 0 = exclude)
  int<lower=0, upper=1> use_t_seas;
  int<lower=0, upper=1> use_t_var;
  int<lower=0, upper=1> use_t_trend;
  int<lower=0, upper=1> use_p_seas;
  int<lower=0, upper=1> use_p_var;
  int<lower=0, upper=1> use_p_trend;
}

parameters {
  // Intercept for mean growth
  real alpha_mu;
  // Intercept for growth variance
  real alpha_sigma;
  // Coefficients for mu
  vector[K_mu] beta;
   // Coefficients for log(sigma)
  vector[K_sigma] gamma;
  // Define lag weights as a simplex
  simplex[L] lag_weights;
}

transformed parameters {
  // Define demographic covariates as sub-vector of betas
  vector[2] beta_demog = beta[1:2];

  // Initialise storage for weighted sums of decomposed climate data
  vector[N] t_seas;
  vector[N] t_var;
  vector[N] t_trend;
  vector[N] p_seas;
  vector[N] p_var;
  vector[N] p_trend;

  // Conditionally compute weighted sums based on toggle parameters
  if (use_t_seas == 1)
    t_seas = t_seas_mat * lag_weights;
  else
    t_seas = rep_vector(0, N);

  if (use_t_var == 1)
    t_var = t_var_mat * lag_weights;
  else
    t_var = rep_vector(0, N);
    
  if (use_t_trend == 1)
    t_trend = t_trend_mat * lag_weights;
  else
    t_trend = rep_vector(0, N);

  if (use_p_seas == 1)
    p_seas = p_seas_mat * lag_weights;
  else
    p_seas = rep_vector(0, N);

  if (use_p_var == 1)
    p_var = p_var_mat * lag_weights;
  else
    p_var = rep_vector(0, N);
    
  if (use_p_trend == 1)
    p_trend = p_trend_mat * lag_weights;
  else
    p_trend = rep_vector(0, N);

  // Compute linear predictor for the mean
  real b_t_seas  = beta[3];
  real b_t_var   = beta[4];
  real b_t_trend = beta[5];
  real b_p_seas  = beta[6];
  real b_p_var   = beta[7];
  real b_p_trend = beta[8];

  vector[N] mu = alpha_mu + X * beta_demog +
                 b_t_seas  * t_seas  +
                 b_t_var   * t_var   +
                 b_t_trend * t_trend +
                 b_p_seas  * p_seas  +
                 b_p_var   * p_var   + 
                 b_p_trend * p_trend;
}

model {
  // Intercept for variance
  alpha_sigma ~ normal(0, 1);
  // Weight variance
  gamma ~ normal(0, 1);
  // Intercept for mean
  alpha_mu ~ normal(0, 1);
  // # Coefficients for mean
  beta ~ normal(0, 1);
  // Simplex prior for lag weights (dirichlet so weights sum to 1)
  lag_weights ~ dirichlet(rep_vector(1, L));
  // Define normal likelihood with non-equal variances
  y ~ normal(mu, exp(alpha_sigma + Z * gamma));
}
