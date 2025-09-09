// ** Model of Mastomys natalensis recruitment (pregnancy) ** //

data {
  int<lower=0> N;  // Number of observations
  int<lower=0, upper=1> y[N];  // Binary response variable
  int<lower=1> P;  // Number of covariates
  int<lower=1> L;  // Number of lag combinations
  matrix[N, L] t_seas_mat;
  matrix[N, L] t_var_mat;
  matrix[N, L] t_trend_mat;
  matrix[N, L] p_seas_mat;
  matrix[N, L] p_var_mat;
  matrix[N, L] p_trend_mat;
  matrix[N, 2] X;  // Design matrix for mu (for weight & pop)
  // Binary switches for predictors (1 = include, 0 = exclude)
  int<lower=0, upper=1> use_t_seas;
  int<lower=0, upper=1> use_t_var;
  int<lower=0, upper=1> use_t_trend;
  int<lower=0, upper=1> use_p_seas;
  int<lower=0, upper=1> use_p_var;
  int<lower=0, upper=1> use_p_trend;
}

parameters {
  real alpha;  // Intercept
  vector[P] beta;  // Coefficients
  simplex[L] lag_weights;  // Lag weight simplex
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
  
  // Construct the linear predictor dynamically
  real b_t_seas  = beta[3];
  real b_t_var   = beta[4];
  real b_t_trend = beta[5];
  real b_p_seas  = beta[6];
  real b_p_var   = beta[7];
  real b_p_trend = beta[8];
  
  vector[N] eta = alpha + X * beta_demog +
                 b_t_seas  * t_seas  +
                 b_t_var   * t_var   +
                 b_t_trend * t_trend +
                 b_p_seas  * p_seas  +
                 b_p_var   * p_var   +
                 b_p_trend * p_trend;
}

model {
  alpha ~ normal(0, 3);
  beta  ~ normal(0, 2);
  lag_weights ~ dirichlet(rep_vector(1, L));
  y ~ bernoulli_logit(eta);
}
