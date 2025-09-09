#-----------------------------------------------------#
# Unscale fitted demographic process model posteriors #
#-----------------------------------------------------#

# 0. Script purpose ----
# - Reads in fitted demographic process model posteriors
# - Since these were fitted to unit-scaled covariates, this script transforms 
#   posteriors to be back on the original scale of data
# - These unscaled posteriors are subsequently used in the integral projection model


# 1. Load packages ----
pacman::p_load(dplyr, ggplot2, lubridate, stringr, tidyr, truncnorm)


# 2. Load data ----

# Scaled demographic data
survival    <- readRDS("data/scaled_survival_data.rds")
recruitment <- readRDS("data/scaled_recruitment_data.rds")
growth      <- readRDS("data/scaled_growth_data.rds")

# Scaled posterior estimates of climate-demographic parameters
post <- readRDS("Outputs/stan_fits/neat_posteriors_alldat_28d.rds")

# Scaled climate data
scaled_env <- readRDS("data/scaled_climate_data.rds")


# 3. Unscale posteriors ----

# Get 3 demographic processes (no male population size)
post_surv <- post %>% filter(process == "surv") %>% group_by(draw) %>% slice(1)
post_rec  <- post %>% filter(process == "rec")  %>% group_by(draw) %>% slice(1)
post_grow <- post %>% filter(process == "grow") %>% group_by(draw) %>% slice(1)

# Un-scale coefficients of survival
post_surv <- post_surv %>% 
  group_by(draw) %>% 
  mutate(
    beta_weight = beta_weight / survival$sd_weight[1],
    beta_pop    = beta_pop    / survival$sd_pop[1],
    beta_ttrend = beta_ttrend / scaled_env$sd_ttrend[1],
    beta_tseas  = beta_tseas  / scaled_env$sd_tseas[1],
    beta_tvar   = beta_tvar   / scaled_env$sd_tvar[1],
    beta_ptrend = beta_ptrend / scaled_env$sd_ptrend[1],
    beta_pseas  = beta_pseas  / scaled_env$sd_pseas[1],
    beta_pvar   = beta_pvar   / scaled_env$sd_pvar[1]
  )

# Un-scale intercept of survival
post_surv <- post_surv %>% 
  group_by(draw) %>% 
  mutate(
    alpha = alpha - sum(
      beta_weight * survival$mean_weight[1], 
        beta_pop * survival$mean_pop[1],
        beta_ttrend * scaled_env$mean_ttrend[1],
        beta_tseas  * scaled_env$mean_tseas[1],
        beta_tvar   * scaled_env$mean_tvar[1],
        beta_ptrend * scaled_env$mean_ptrend[1],
        beta_pseas  * scaled_env$mean_pseas[1],
        beta_pvar   * scaled_env$mean_pvar[1], 
      na.rm = TRUE
    )
  )

# Un-scale coefficients of recruitment
post_rec <- post_rec %>% 
  group_by(draw) %>% 
  mutate(
    beta_weight = beta_weight / recruitment$sd_weight[1],
    beta_pop    = beta_pop    / recruitment$sd_pop[1],
    beta_ttrend = beta_ttrend / scaled_env$sd_ttrend[1],
    beta_tseas  = beta_tseas  / scaled_env$sd_tseas[1],
    beta_tvar   = beta_tvar   / scaled_env$sd_tvar[1],
    beta_ptrend = beta_ptrend / scaled_env$sd_ptrend[1],
    beta_pseas  = beta_pseas  / scaled_env$sd_pseas[1],
    beta_pvar   = beta_pvar   / scaled_env$sd_pvar[1]
  )

# Un-scale intercept of recruitment
post_rec <- post_rec %>% 
  group_by(draw) %>% 
  mutate(
    alpha = alpha - sum(
      beta_weight   * recruitment$mean_weight[1],
        beta_pop    * recruitment$mean_pop[1],
        beta_ttrend * scaled_env$mean_ttrend[1],
        beta_tseas  * scaled_env$mean_tseas[1],
        beta_tvar   * scaled_env$mean_tvar[1],
        beta_ptrend * scaled_env$mean_ptrend[1],
        beta_pseas  * scaled_env$mean_pseas[1],
        beta_pvar   * scaled_env$mean_pvar[1], 
      na.rm = TRUE
    )
  )

# Un-scale coefficients of growth
post_grow <- post_grow %>% 
  group_by(draw) %>% 
  mutate(
    # mean
    beta_weight = growth$sd_W2[1] / growth$sd_W1[1]         * beta_weight,
    beta_pop    = growth$sd_W2[1] / growth$sd_pop[1]        * beta_pop,
    beta_ttrend = growth$sd_W2[1] / scaled_env$sd_ttrend[1] * beta_ttrend,
    beta_tseas  = growth$sd_W2[1] / scaled_env$sd_tseas[1]  * beta_tseas,
    beta_tvar   = growth$sd_W2[1] / scaled_env$sd_tvar[1]   * beta_tvar,
    beta_ptrend = growth$sd_W2[1] / scaled_env$sd_ptrend[1] * beta_ptrend,
    beta_pseas  = growth$sd_W2[1] / scaled_env$sd_pseas[1]  * beta_pseas,
    beta_pvar   = growth$sd_W2[1] / scaled_env$sd_pvar[1]   * beta_pvar,
    # variance
    beta_sigma  = beta_sigma / growth$sd_W1[1]
  )

# Un-scale intercepts of growth
post_grow <- post_grow %>% 
  group_by(draw) %>% 
  mutate(
    # mean
    alpha_mu = growth$mean_W2[1] + growth$sd_W2[1] * alpha_mu - sum(
      beta_weight   * growth$mean_W1[1],
        beta_pop    * growth$mean_pop[1],
        beta_ttrend * scaled_env$mean_ttrend[1],
        beta_tseas  * scaled_env$mean_tseas[1],
        beta_tvar   * scaled_env$mean_tvar[1],
        beta_ptrend * scaled_env$mean_ptrend[1],
        beta_pseas  * scaled_env$mean_pseas[1],
        beta_pvar   * scaled_env$mean_pvar[1],
      na.rm = TRUE
    ),
    # variance
    alpha_sigma = log(growth$sd_W2[1]) + alpha_sigma - sum(
      beta_sigma * growth$mean_W1[1],
      na.rm = TRUE
    )
  )

# Join unscaled posteriors together, re-join to climate lags & male population posterior
unscaled_posts <- bind_rows(post_surv, post_rec, post_grow) %>% 
  right_join(post %>% filter(process != "male_pop") %>% reframe(draw, process, lag_wgt_var, lag_wgt_val, t_lag, p_lag)) %>% 
  group_by(process, draw) %>% 
  fill(alpha, beta_weight, beta_pop, beta_ttrend, beta_tseas, beta_tvar, 
       beta_ptrend, beta_pseas, beta_pvar, beta_sigma, alpha_mu, alpha_sigma, .direction = "down")

male_post <- post %>% filter(process == "male_pop")

all_unscaled_posts <- bind_rows(male_post, unscaled_posts)
# dim(all_unscaled_posts)

# Save unscaled posteriors for subsequent analyses
saveRDS(all_unscaled_posts, "Outputs/unscaled_posteriors.rds")
