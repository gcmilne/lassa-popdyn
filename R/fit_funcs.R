#-----------------------------------------------------#
# Define functions for transmission parameter fitting #
#-----------------------------------------------------#

# 0. Script purpose ----
# - Defines functions for estimating transmission parameters of integral 
#   projection model using an approximate Bayesian procedure
# - Functions sourced in transmission fitting scripts


# Function to create a date vector given a start date & target end date
getDateVec <- function(start_date, end_sampling, dt) {
  
  # start_date: date to start vector from, in ymd format
  # end_sampling: target date to end vector, in ymd format
  # dt: no. days between each date element (numeric input)
  
  # Create date vector given inputs
  date_vec <- seq(start_date, end_sampling, by = dt)
  
  # Append extra dates to end of vector if sampling end date not within vector range
  while (max(date_vec) < end_sampling) {
    date_vec <- c(date_vec, max(date_vec) + dt)
  }
  
  # Return date vector
  return(date_vec)
  
}

# Specify transmission parameter priors
prior_trans_rate     <- function(x) dunif(x, 0.00001, 0.0015, log = TRUE)
prior_vertical_trans <- function(x) dbeta(x, shape1 = 20, shape2 = 2, log = TRUE) # Refs: 10.1093/pnasnexus/pgac114 & 10.1098/rsif.2024.0106

# set.seed(1)
# tmp <- rbeta(10000, 20, 2)
# quantile(tmp, c(0.025, 0.975))

# Function to run rodent model given some parameter values & calculate loglik given seroprevalence data
calcLogLik <- function(seroprev_data, stat_model, env_data, cmr_data, demog_pars, 
                       inf_pars, n_init, start_date, end_date, seed, sens = 1, spec = 1) {
  
  # Run model with given parameter values to calculate log-likelihood
  mod <- RunRodentModel(model = stat_model, env = env_data, cmr.data = cmr_data,
                        demog.pars = demog_pars, inf.pars = inf_pars,
                        n_init = n_init, start_date = start_date, end_date = end_date,
                        seed = seed, summarised = FALSE) %>%
    # Store time-specific modelled population size & seroprevalence (subadults & adults)
    group_by(time) %>%
    reframe(date, mod.pop, seroprev.sub, seroprev.adult) %>%
    distinct()
  
  # Get dates from data over which to interpolate model estimates (to match model dates to data dates)
  dates_df <- data.frame(date = as.Date(seq(min(seroprev_data$session.start),
                                            max(seroprev_data$session.start), by = 1)))
  
  # Perform right join to keep all dates in data
  result <- left_join(dates_df, mod, by = "date") %>%
    mutate(
      # Create flag for interpolation
      is.interpolated = case_when(
        is.na(mod.pop) ~ TRUE,
        .default = FALSE
      ),
      # Interpolate modelled seroprevalence
      mod.seroprev.sub.int = zoo::na.approx(seroprev.sub, na.rm = FALSE),
      mod.seroprev.ad.int  = zoo::na.approx(seroprev.adult, na.rm = FALSE)
    )
  
  # Put in long format for easier joining with data
  result_long <- result %>% 
    filter(!is.na(mod.seroprev.sub.int), !is.na(mod.seroprev.ad.int)) %>% 
    dplyr::select(date, contains(c("mod.seroprev"))) %>% 
    pivot_longer(!date, values_to = "mod.seroprev", names_to = "var") %>% 
    mutate(FER = if_else(str_detect(var, "sub"), FALSE, TRUE)) %>% 
    dplyr::select(!var)
  
  # Join to data
  mod.dat <- result_long %>%
    right_join(., seroprev_data, by = c("date" = "session.start", "FER")) %>%
    filter(!is.na(mod.seroprev))
  
  # Empty-case handling if NA in seroprevalence: empty tibble w/ same columns as normal case
  if (nrow(mod.dat) == 0) {
    empty_moddat <- tibble(
      date            = as.Date(character()),
      mod.seroprev    = numeric(),
      FER             = logical(),
      date_mid        = as.Date(character()),
      seropos         = numeric(),
      tested          = integer(),
      seroprev        = numeric(),
      mean.seroprev   = numeric(),
      lower.seroprev  = numeric(),
      upper.seroprev  = numeric(),
      loglik          = numeric()
    )
    return(empty_moddat)
  }

  # Calculate & return binomial log-likelihood
  sub_mod_seroprev <- mod.dat$mod.seroprev[mod.dat$FER == FALSE]
  ad_mod_seroprev <- mod.dat$mod.seroprev[mod.dat$FER == TRUE]
  
  loglik_sub <- sum(dbinom(  # sub-adults
    x = mod.dat$seropos[mod.dat$FER == FALSE],         # Data: no. seropositive
    size = mod.dat$tested[mod.dat$FER == FALSE],       # Data: no. tested
    prob = sens*sub_mod_seroprev + (1-spec)*(1-sub_mod_seroprev), # Modelled seroprevalence and sensitivity/specificity
    log = TRUE))
  
  loglik_ad <- sum(dbinom(  # adults
    x = mod.dat$seropos[mod.dat$FER == TRUE],         # Data: no. seropositive
    size = mod.dat$tested[mod.dat$FER == TRUE],       # Data: no. tested
    prob = sens*ad_mod_seroprev + (1-spec)*(1-ad_mod_seroprev), # Modelled seroprevalence and sensitivity/specificity
    log = TRUE))
  
  loglik <- loglik_sub + loglik_ad
  
  return(cbind(mod.dat, loglik))
  
}


# Define function to run grid search algorithm to calculate model fit for each pair of parameter values
GridSearch <- function(alpha_vals, beta_vals, seroprev.dat, demog.posts, morogoro.clim,
                       cmr.data, demog.pars, inf.pars, n_init, start_date, end_date, seed,
                       sens = 1, spec = 1,
                       loglik_func, alpha_prior_func, beta_prior_func, make_grid = TRUE) {

  # Initialize storage for parameter sets and log-likelihoods
  if (make_grid) {  # all combinations of parameter sets
    grid_results <- expand.grid(vertical_trans = alpha_vals, trans_rate = beta_vals)  # length = length(alpha_vals)^2
  } else {  # use pre-set parameter combinations
    grid_results <- data.frame(vertical_trans = alpha_vals, trans_rate = beta_vals)  # length = length(alpha_vals)
  }
  grid_results$loglik <- NA  # log-likelihood
  
  # Calculate posterior density for each parameter set
  grid_out <- purrr::map_dfr(1:nrow(grid_results), function(i) {
    
    # Set parameters for this grid point
    inf.pars$vertical_trans <- grid_results$vertical_trans[i]
    inf.pars$trans_rate     <- grid_results$trans_rate[i]
    
    # Calculate joint log-likelihood for current parameter set (& store modelled seroprevalence)
    fitted_mod <- loglik_func(
      seroprev_data = seroprev.dat,
      stat_model = demog.posts,
      env_data = morogoro.clim,
      cmr_data = cmr.data,
      demog_pars = demog.pars,
      inf_pars = inf.pars,
      n_init = n_init,
      start_date = start_date,
      end_date = end_date,
      seed = seed,
      sens = sens,
      spec = spec
    )
    
    # Extract log-likelihood, prior, posterior
    current_loglik <- fitted_mod$loglik[1]
    current_prior <- beta_prior_func(inf.pars$trans_rate) + alpha_prior_func(inf.pars$vertical_trans)
    posterior <- current_loglik + current_prior
    
    # Return everything
    tibble::tibble(
      vertical_trans  = inf.pars$vertical_trans,
      trans_rate      = inf.pars$trans_rate,
      joint_prior     = current_prior,
      joint_loglik    = current_loglik,
      joint_posterior = posterior,
      i               = i
    )
    
  })
  
  return(grid_out)
  
}
