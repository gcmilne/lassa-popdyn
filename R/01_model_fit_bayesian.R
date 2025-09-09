#---------------------------------------------------------#
# Fit three process models given some demographic & climate data #
#---------------------------------------------------------#

# 0. Session details ----
# - Load rodent trapping (demographic) data & climate data
# - Find posterior distributions for model parameters for each demographic process 
# - Format & save posterior distributions for future workflow


# 1. Load packages ----
pacman::p_load(ggplot2, lubridate, dplyr, stringr, tidyr, rstan, bayesplot, tidybayes, 
               ggridges, patchwork, brms, loo, bridgesampling, boot, doParallel)

# Source custom ggplot theme
source("R/ggplot_theme.R")


# 2. Load data ----

# Function to apply unit scaling to covariates (mean of 0, SD of 1)
unitScale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Demographic data (from 00demo_formatting_Morogoro.R)
survival <- read.csv("Data/survival_data.csv", stringsAsFactors = FALSE) %>% 
  mutate(DATE = ymd(DATE),
         mean_pop = mean(POP),
         sd_pop = sd(POP),
         mean_weight = mean(WEIGHT),
         sd_weight = sd(WEIGHT),
         across(c(POP, WEIGHT), unitScale))

recruitment <- read.csv("Data/recruitment_data.csv", stringsAsFactors = FALSE) %>%
  mutate(DATE = ymd(DATE),
         mean_pop = mean(POP),
         sd_pop = sd(POP),
         mean_weight = mean(WEIGHT),
         sd_weight = sd(WEIGHT),
         across(c(POP, WEIGHT), unitScale))

growth <- read.csv("Data/growth_data.csv", stringsAsFactors = FALSE) %>%
  mutate(DATE1 = ymd(DATE1), DATE2 = ymd(DATE2),
         mean_pop = mean(POP),
         sd_pop = sd(POP),
         mean_W1 = mean(W1),
         sd_W1 = sd(W1),
         mean_W2 = mean(W2),
         sd_W2 = sd(W2),
         across(c(POP, W1, W2), unitScale))

# Minimum number alive (MNA) data (from demo_formatting_Morogoro.R)
MNA <- read.csv("Data/MNA.by.sex.csv") %>% 
  mutate(month = month(trap.date))

# Climatic data (from 00env_formatting.R)
env <- readRDS("Data/monthly.climate.1990-2024_separate_trend.rds") %>%
  filter(state == "Morogoro") %>% 
  ungroup() %>%
  dplyr::select(date, contains(c("seas", "var", "trend")))

# Apply unit scaling to climate data
env_data <- env %>%
  mutate(
    mean_tseas = mean(t_seas),
    sd_tseas   = sd(t_seas),
    mean_tvar  = mean(t_var),
    sd_tvar    = sd(t_var),
    mean_pseas = mean(p_seas),
    sd_pseas   = sd(p_seas),
    mean_pvar  = mean(p_var),
    sd_pvar    = sd(p_var),
    mean_ttrend = mean(t_trend),
    sd_ttrend = sd(t_trend),
    mean_ptrend = mean(p_trend),
    sd_ptrend = sd(p_trend),
    across(c(t_seas, t_var, p_seas, p_var, t_trend, p_trend), unitScale))

saveRDS(env_data,    "data/scaled_climate_data.rds")
saveRDS(survival,    "data/scaled_survival_data.rds")
saveRDS(recruitment, "data/scaled_recruitment_data.rds")
saveRDS(growth,      "data/scaled_growth_data.rds")

# Set seeds for reproducible data sampling & model fitting
sample_seed <- 1
model_seed  <- 1234


# 3. Define custom functions ----

# Function to generate lags in climate data in 2 dimensions
genBothLags <- function(demog_dat, clim_dat, lags, date_col) {
  
  # demog_dat: Demographic data 
  # clim_dat: Climate data
  # lags: climate lags (in days) in order: (1) temperature, (2) precipitation
  # date_col: column name for date in the data (used to join to lagged climate data)
  
  # Generate lagged temp data
  clim_lag_t <- clim_dat %>% 
    reframe(date) %>% 
    mutate(lag_date_t = ymd(date) - lags[1]) %>% 
    left_join(clim_dat, by = c("lag_date_t" = "date")) %>% 
    dplyr::select(date, t_seas, t_var, t_trend) %>% 
    na.omit() 
  
  # Generate lagged precipitation data
  clim_lag_p <- clim_dat %>% 
    reframe(date) %>% 
    mutate(lag_date_p = ymd(date) - lags[2]) %>% 
    left_join(clim_dat, by = c("lag_date_p" = "date")) %>% 
    dplyr::select(date, p_seas, p_var, p_trend) %>% 
    na.omit() 
  
  # Join both climate data together
  both_lag_clim <- full_join(clim_lag_t, clim_lag_p, by = "date") %>% na.omit()
  
  # Join to demographic data & select relevant columns
  demog_clim <- demog_dat %>% 
    left_join(both_lag_clim, by = setNames("date", date_col)) %>%
    reframe(t_seas, t_var, t_trend, p_seas, p_var, p_trend)
  
  # Return
  return(demog_clim)
  
}


# Function to generate lag array for all combination of lags for each climatic variable
genLagArray <- function(demog_dat, clim_dat, lag_grid, num_vars = 6) {
  
  # demog_dat: Demographic data 
  # clim_dat: Climate data
  # lag_grid: Grid of climate lags (in days) for temperature and precipitation
  # num_vars: number of climate variables
  
  # Set name for date column depending on input demographic data
  if (any(str_detect(colnames(demog_dat), "W2"))) {  # normal model (growth)
    date_col <- "DATE1"
  } else {  # binomial model (survival & pregnancy)
    date_col <- "DATE"
  }
  
  # For each lag combination, calculate each variable
  lags_array <- array(NA, dim = c(nrow(demog_dat), nrow(lag_grid), num_vars))  # Dimensions: N x (T*P lags) x 6
  
  # Populate the array for all variables & all lags
  for (i in 1:nrow(lag_grid)) {
    
    lagged_data <- genBothLags(demog_dat = demog_dat, clim_dat = clim_dat, lags = as.numeric(lag_grid[i, c("t_lag", "p_lag")]), date_col = date_col)
    
    lags_array[ , i, 1] <- lagged_data$t_seas  # Values for t_seas
    lags_array[ , i, 2] <- lagged_data$t_var   # Values for t_var
    lags_array[ , i, 3] <- lagged_data$t_trend # Values for t_trend
    lags_array[ , i, 4] <- lagged_data$p_seas  # Values for p_seas
    lags_array[ , i, 5] <- lagged_data$p_var   # Values for p_var
    lags_array[ , i, 6] <- lagged_data$p_trend # Values for p_trend
    
  }
  
  return(lags_array)
  
}


# 4. Fit demographic models ----

# Aim: fit each demographic process model with every possible combination of climate
#      covariates & then calculate marginal log-likelihood of each model

## (a) Compile Stan models ----
rec_mod  <- stan_model(file = "stan/recruitment_model_toggle.stan")  # for recruitment
grow_mod <- stan_model(file = "stan/grow_model_toggle.stan")  # for growth
surv_mod <- stan_model(file = "stan/surv_model_toggle.stan")  # for survival

# Define parallel fitting function
fitMods <- function(full_stan_data, clim_toggles, stan_mod, niter, nchains, ncores, 
                    seed, i_vec, nwarmup = floor(niter/2), model_name) {
  
  # Set no. parallel workers to not exceed computer's maximum no. cores
  num_workers <- parallel::detectCores()/ncores 
  cl <- makeCluster(num_workers)
  
  # Register parallel back-end
  registerDoParallel(cl)
  
  # Fit models in parallel
  foreach(i = i_vec, .packages = c("rstan", "bridgesampling", "tidyverse")) %dopar% {
    
    # Set climate toggles
    stan_dat <- full_stan_data
    stan_dat$use_t_seas  <- clim_toggles$use_t_seas[i]
    stan_dat$use_t_var   <- clim_toggles$use_t_var[i]
    stan_dat$use_t_trend <- clim_toggles$use_t_trend[i]
    stan_dat$use_p_seas  <- clim_toggles$use_p_seas[i]
    stan_dat$use_p_var   <- clim_toggles$use_p_var[i]
    stan_dat$use_p_trend <- clim_toggles$use_p_trend[i]
    
    # Fit model
    mod <- sampling(object = stan_mod, data = stan_dat, iter = niter, warmup = nwarmup,
                    chains = nchains, cores = ncores, seed = seed, verbose = FALSE,
                    include = FALSE, pars = c("t_seas", "t_var", "t_trend", "p_seas", "p_var", "p_trend",
                                              "beta_demog", "b_t_seas", "b_t_var", "b_t_trend",
                                              "b_p_seas", "b_p_var", "b_p_trend"))
    
    # Compute log marginal likelihood
    bridge_mod <- bridge_sampler(mod)

    # File paths
    model_path  <- paste0("Outputs/stan_fits/", model_name, "/fits/mod_", i, ".rds")
    bridge_path <- paste0("Outputs/stan_fits/", model_name, "/margLL/mod_", i, ".rds")
    
    # Save each to disk
    saveRDS(mod, model_path)
    saveRDS(bridge_mod, bridge_path)

  }
  
  # Stop the cluster after execution
  stopCluster(cl)
  
}

## (b) Create lagged climatic data ----

# Generate grid of climate lags
lag_int  <- 28
lag_seq  <- seq(0, 168, by = lag_int)
lag_grid <- setNames(expand.grid(lag_seq, lag_seq), c("t_lag", "p_lag")) %>% as.data.frame()

# Generate lagged climate data at dates of different demographic data
lag_arr_grow <- genLagArray(demog_dat = growth, clim_dat = env_data, lag_grid = lag_grid, num_vars = 6)
lag_arr_rec  <- genLagArray(demog_dat = recruitment, clim_dat = env_data, lag_grid = lag_grid, num_vars = 6)
lag_arr_surv <- genLagArray(demog_dat = survival, clim_dat = env_data, lag_grid = lag_grid, num_vars = 6)

# Get all possible combinations of climate covariates
clim_toggles <- crossing(use_t_seas = 0:1, use_t_var = 0:1, use_t_trend = 0:1, 
                         use_p_seas = 0:1, use_p_var = 0:1, use_p_trend = 0:1)

## (c) Initialise Stan data lists ----

# Growth
grow_stan <- list(N           = nrow(growth),
                  y           = growth$W2,
                  K_mu        = 8,
                  K_sigma     = 1,
                  X           = growth %>% dplyr::select(W1, POP),
                  Z           = growth %>% dplyr::select(W1),
                  L           = nrow(lag_grid),
                  t_seas_mat  = lag_arr_grow[, , 1] %>% bind_cols(),
                  t_var_mat   = lag_arr_grow[, , 2] %>% bind_cols(),
                  t_trend_mat = lag_arr_grow[, , 3] %>% bind_cols(),
                  p_seas_mat  = lag_arr_grow[, , 4] %>% bind_cols(),
                  p_var_mat   = lag_arr_grow[, , 5] %>% bind_cols(),
                  p_trend_mat = lag_arr_grow[, , 6] %>% bind_cols(),
                  use_t_seas  = 1,
                  use_t_var   = 1,
                  use_t_trend = 1,
                  use_p_seas  = 1,
                  use_p_var   = 1,
                  use_p_trend = 1)

# Survival
surv_stan <- list(N          = nrow(survival),
                  y          = survival$SURV,
                  P          = 8,
                  L          = nrow(lag_grid),
                  t_seas_mat  = lag_arr_surv[, , 1] %>% bind_cols(),
                  t_var_mat   = lag_arr_surv[, , 2] %>% bind_cols(),
                  t_trend_mat = lag_arr_surv[, , 3] %>% bind_cols(),
                  p_seas_mat  = lag_arr_surv[, , 4] %>% bind_cols(),
                  p_var_mat   = lag_arr_surv[, , 5] %>% bind_cols(),
                  p_trend_mat = lag_arr_surv[, , 6] %>% bind_cols(),
                  X          = survival %>% dplyr::select(WEIGHT, POP),
                  use_t_seas = 1,
                  use_t_var  = 1,
                  use_t_trend = 1,
                  use_p_seas = 1,
                  use_p_var  = 1,
                  use_p_trend = 1)

# Recruitment
rec_stan <- list(N           = nrow(recruitment),
                 y           = recruitment$REC,
                 P           = 8,
                 L           = nrow(lag_grid),
                 t_seas_mat  = lag_arr_rec[, , 1] %>% bind_cols(),
                 t_var_mat   = lag_arr_rec[, , 2] %>% bind_cols(),
                 t_trend_mat = lag_arr_rec[, , 3] %>% bind_cols(),
                 p_seas_mat  = lag_arr_rec[, , 4] %>% bind_cols(),
                 p_var_mat   = lag_arr_rec[, , 5] %>% bind_cols(),
                 p_trend_mat = lag_arr_rec[, , 6] %>% bind_cols(),
                 X           = recruitment %>% dplyr::select(WEIGHT, POP),
                 use_t_seas  = 1,
                 use_t_var   = 1,
                 use_t_trend = 1,
                 use_p_seas  = 1,
                 use_p_var   = 1,
                 use_p_trend = 1)

## (d) Fit models with different climate covariates ----

### Fit growth model ----
system.time( fitMods(full_stan_data = grow_stan, clim_toggles = clim_toggles, stan_mod = grow_mod, 
                     niter = 2000, nchains = 4, ncores = 4, seed = model_seed, model_name = "growth",
                     i_vec = 1:nrow(clim_toggles)) )

### Fit recruitment model ----
system.time( fitMods(full_stan_data = rec_stan, clim_toggles = clim_toggles, stan_mod = rec_mod, 
                     niter = 2000, nchains = 4, ncores = 4, seed = model_seed, model_name = "recruitment",
                     i_vec = 1:nrow(clim_toggles)) )

### Fit survival model ----
system.time( fitMods(full_stan_data = surv_stan, clim_toggles = clim_toggles, stan_mod = surv_mod, 
                     niter = 2000, nchains = 4, ncores = 4, seed = model_seed, model_name = "survival",
                     i_vec = 1:nrow(clim_toggles)) )  # 4.5 hours for 15 fits


## (e) Extract model convergence statistics ----

# Extract Rhat values to monitor model convergence
models <- c("growth", "recruitment", "survival")

# Create a named list-of-lists to store Rhat checks
rhats_pass <- vector("list", length(models))
names(rhats_pass) <- models

for (m in seq_along(models)) {
  
  print(paste("Extract rhat for", models[m], "model"))
  
  rhats <- logical(nrow(clim_toggles))
  
  for (i in seq_along(rhats)) {
    
    # Read one file at a time to manage memory
    model_path <- paste0("Outputs/stan_fits/",  models[m], "/fits/mod_", i, ".rds")
    mod_fit <- readRDS(model_path)
    
    summary_df <- summary(mod_fit)$summary
    rhats[i] <- all(summary_df[, "Rhat"] < 1.05)
    
    print(paste("iteration", i, "out of", nrow(clim_toggles), "complete"))
  }
  
  rhats_pass[[m]] <- rhats
  
  # Free memory
  rm(mod_fit, summary_df)
  gc()
}

saveRDS(rhats_pass, "Outputs/stan_fits/rhats_all_mods.rds")

# Turn into dataframe
tmp <- vector("list", length(models))
for (m in seq_along(models)){
  tmp[[m]] <- data.frame(
    model = models[m],
    run = seq_along(rhats_pass[[m]]),
    is_converged = rhats_pass[[m]]
  )
}
rhat_df <- bind_rows(tmp)

# Save
# saveRDS(rhat_df, "Outputs/stan_fits/rhats_all_mods.rds")
rhat_df <- readRDS("Outputs/stan_fits/rhats_all_mods.rds")


# 5. Plot marginal log-likelihoods of fitted demographic models ----
LL_dat <- NULL

# Load marginal log-likelihoods
for (m in seq_along(models)) {
  
  LL_tmp <- data.frame(model = rep(models[m], nrow(clim_toggles)), LL = rep(NA, nrow(clim_toggles)), run = 1:nrow(clim_toggles))
  
  for (i in 1:nrow(LL_tmp)) {
    # Read one file at a time to manage memory
    model_path   <- paste0("Outputs/stan_fits/",  models[m], "/margLL/mod_", i, ".rds")
    bridge_mod   <- readRDS(model_path)
    LL_tmp$LL[i] <- bridge_mod$logml
  }
  
  LL_dat <- bind_rows(LL_dat, LL_tmp)
}
rm(LL_tmp)

# Data to show if model has converged or not
rhat_df <- data.frame(is_converged = c(rhat_df$growth, rhat_df$recruitment, rhat_df$survival),
                      model = c(
                        rep("growth", length(rhat_df$growth)),
                        rep("recruitment", length(rhat_df$recruitment)),
                        rep("survival", length(rhat_df$survival))
                      )) %>% 
  group_by(model) %>% 
  mutate(run = row_number())

# Format data for plotting
LL_dat <- LL_dat %>% 
  left_join(clim_toggles %>% mutate(run = row_number())) %>% 
  left_join(rhat_df) %>% 
  group_by(model) %>% 
  mutate(LL_dif = LL - min(LL, na.rm = TRUE),
         bf_max = exp(max(LL, na.rm = TRUE) - LL)) %>% 
  mutate(model = case_when(
    model == "growth"      ~ "Growth",
    model == "recruitment" ~ "Pregnancy",
    model == "survival"    ~ "Survival"
  ),
  is_converged = if_else(is_converged == FALSE, "No", "Yes")) %>% 
  group_by(model) %>% 
  mutate(
    stars = case_when(
      !is.na(bf_max) & bf_max < 5 & bf_max != 1  ~ "**",
      !is.na(bf_max) & bf_max < 10 & bf_max >= 5 ~ "*",
      TRUE                                       ~ ""
    )
  )

# Make plots for each demographic process
for (i in unique(LL_dat$model)) {
  
  # Create a factor for run ordering based on a baseline model LL ranking
  baseline_order <- LL_dat %>%
    filter(model == i) %>%
    arrange(desc(LL)) %>%
    pull(run)
  
  # Set run as a factor with levels based on the growth order
  mod_LL_dat <- LL_dat %>%
    filter(model == i) %>% 
    arrange(desc(LL)) %>% 
    mutate(run_factor = factor(run, levels = baseline_order))
  
  # Reshape toggle parameters into long format for geom_tile()
  LL_dat_long <- mod_LL_dat %>%
    pivot_longer(cols = contains(c("seas", "var", "trend")), names_to = "covar", values_to = "include") %>% 
    mutate(run_factor = factor(run, levels = baseline_order),
           covar = case_when(
             str_detect(covar, "t_seas")  ~ "Temperature \n seasonality",
             str_detect(covar, "t_var")   ~ "Temperature \n variability",
             str_detect(covar, "t_trend") ~ "Temperature \n trend",
             str_detect(covar, "p_seas")  ~ "Precipitation \n seasonality",
             str_detect(covar, "p_var")   ~ "Precipitation \n variability",
             str_detect(covar, "p_trend") ~ "Precipitation \n trend"))
  
  # Model structure heatmap
  heatmap_plot <- LL_dat_long %>% 
    mutate(include = if_else(include==0, "No", "Yes")) %>% 
    ggplot(aes(x = run_factor, y = covar, fill = factor(include))) +
    geom_tile(color = "white") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "bottom") +
    scale_fill_viridis_d(option = "cividis", direction = -1) + 
    labs(x = "Model variant", y = "", fill = "Covariate included")
  
  # If all models converged, make all bars same colour
  if (all(mod_LL_dat$is_converged == "Yes")) bar_fill <- "black" else bar_fill <- c("lightgrey", "black")
  
  # Marginal Log-Likelihood bar plot
  dodge_width <- 0.9
  bar_plot <- mod_LL_dat %>%
    ggplot(aes(x = run_factor, y = LL_dif, fill = is_converged)) +
    geom_col(position = position_dodge(width = dodge_width)) +
    # Asterisks
    geom_text(aes(label = stars),
              position = position_dodge(width = dodge_width), vjust = -0.2, size = 3.5, color = "red") +
    labs(x = "", y = "Difference from minimum log marginal likelihood", fill = "Model converged") +
    scale_fill_manual(values = bar_fill) +
    facet_wrap(~model) + 
    coord_cartesian(clip = "off")
  
  # Combine & save plot
  print ( bar_plot / heatmap_plot + 
            plot_layout(heights = c(2, 1)) ) & # Adjust height proportions
    theme_Publication() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(paste0("Plots/model_fits/", i, "_climcovars.pdf"), device = cairo_pdf, height = 8, width = 12, units = "in")
  ggsave(paste0("Plots/model_fits/", i, "_climcovars.png"), bg = "white", height = 8, width = 12, units = "in", dpi = 1200)
}


# 6. Fit male population size as function of female pop size ----

# Model 1: fixed intercept & coefficient
brm1 <- brm(mna.male | trunc(lb = 0) ~ mna.female, data = MNA, family = gaussian(),
            chains = 4, cores = 4, iter = 2000, seed = model_seed)
pp_check(brm1)
loo1 <- loo(brm1)

# Model 2: random intercept for month, fixed coefficient
brm2 <- brm(mna.male | trunc(lb = 0) ~ mna.female + (1|month), data = MNA, family = gaussian(),
            chains = 4, cores = 4, iter = 2000, seed = model_seed)
pp_check(brm2)
loo2 <- loo(brm2)

# Model 3: fixed intercept, random coefficient for month
brm3 <- brm(mna.male | trunc(lb = 0) ~ mna.female + month, data = MNA, family = gaussian(),
            chains = 4, cores = 4, iter = 2000, seed = model_seed)
pp_check(brm3)
loo3 <- loo(brm3)

# Model 4: random intercept, random coefficient for month
brm4 <- brm(mna.male | trunc(lb = 0) ~ mna.female + month + (1|month), data = MNA, family = gaussian(),
            chains = 4, cores = 4, iter = 2000, seed = model_seed)
pp_check(brm4)
loo4 <- loo(brm4)

# Compare all models
loo_compare(loo1, loo2, loo3, loo4) # choose simplest model (model 2)

# Save fit of best-fitting model (random intercept)
saveRDS(brm2, file = "Outputs/stan_fits/mf_pop.rds")


# 7. Format & save posteriors of best-fitting models ----

# Function to rename posterior parameters with nicer names
renameVars <- function(data) {
  data %>%
    rename_with(
      ~ case_when(
        . == "beta[1]"  ~ "beta_weight",
        . == "beta[2]"  ~ "beta_pop",
        . == "beta[3]"  ~ "beta_tseas",
        . == "beta[4]"  ~ "beta_tvar",
        . == "beta[5]"  ~ "beta_ttrend",
        . == "beta[6]"  ~ "beta_pseas",
        . == "beta[7]"  ~ "beta_pvar",
        . == "beta[8]"  ~ "beta_ptrend",
        . == "gamma[1]" ~ "beta_sigma",
        . == ".draw"    ~ "draw",
        TRUE ~ .
      )
    )
}

# Function to extract coefficients of best-performing models
get_active_betas <- function(clim_row) {
  beta_names <- c("beta_tseas", "beta_tvar", "beta_ttrend", 
                  "beta_pseas", "beta_pvar", "beta_ptrend")
  
  # Match toggle names to beta names
  toggle_map <- c("use_t_seas", "use_t_var", "use_t_trend", 
                  "use_p_seas", "use_p_var", "use_p_trend")
  
  # Get active ones (toggle == 1)
  keep <- beta_names[as.logical(clim_row[toggle_map])]
  
  # Always include base betas
  c("beta_weight", "beta_pop", keep)
}

# Function to join formatted posterior to climate lag data
joinPostClim <- function(posterior, lag_grid) {
  posterior %>%
    # Select model coefficients
    select(draw, process, contains(c("alpha", "beta", "lag_weights")) & !contains("demog") & !contains("clim")) %>%
    # Add temperature & precipitation lag values, matching by lag_weight naming
    pivot_longer(contains("lag_weights"), names_to = "lag_wgt_var", values_to = "lag_wgt_val") %>%
    left_join(lag_grid)
}

# Get best-fitting models for each process (removing non-converged models)
best_fit <- LL_dat %>% filter(is_converged == "Yes") %>% group_by(model) %>% filter(LL == max(LL))

best_mod <- as.list(rep(NA, 3))  # best-fitting model for each process
for (i in 1:nrow(best_fit)) {
  best_mod[[i]] <- readRDS(paste0("Outputs/stan_fits/", models[i], "/fits/", "mod_", best_fit$run[i], ".rds"))
}

rec_fit  <- best_mod[[which(best_fit$model == "Pregnancy")]]
surv_fit <- best_mod[[which(best_fit$model == "Survival")]]
grow_fit <- best_mod[[which(best_fit$model == "Growth")]]
rm(best_mod)  # free up storage

# Save best-fitting models
saveRDS(rec_fit,  "Outputs/stan_fits/surv_stan_bestfit.rds")
saveRDS(surv_fit, "Outputs/stan_fits/rec_stan_bestfit.rds")
saveRDS(grow_fit, "Outputs/stan_fits/grow_stan_bestfit.rds")

# Edit lag grid to join posteriors to climate lags by variable matching in `joinPostClim()`
lag_grid <- lag_grid %>%
  mutate(lag_wgt_var = paste0("lag_weights[", row_number(), "]"))

# Function to extract parameters from best fitting model & return tidy posterior
tidyPosterior <- function(fit, model_name, process_name, best_fit, clim_toggles, lag_grid) {
  
  run_id <- best_fit %>% filter(model == model_name) %>% pull(run)
  clim_row <- clim_toggles[run_id, ]
  
  keep_betas <- get_active_betas(clim_row)
  if (process_name == "grow") keep_betas <- c(keep_betas, "beta_sigma")
  
  posterior <- fit %>% 
    tidybayes::tidy_draws() %>% 
    mutate(process = process_name) %>%
    renameVars() %>% 
    select(draw, process, all_of(keep_betas), contains(c("alpha", "lag_weights"))) %>%
    joinPostClim(posterior = ., lag_grid = lag_grid)
  
  return(posterior)
}

# Get tidy posterior distributions of best-performing models for each demographic process
surv_post  <- tidyPosterior(fit = surv_fit,  model_name = "Survival", process_name = "surv", 
                            best_fit = best_fit, clim_toggles = clim_toggles, lag_grid = lag_grid)
rec_post  <- tidyPosterior(fit = rec_fit,  model_name = "Pregnancy", process_name = "rec", 
                            best_fit = best_fit, clim_toggles = clim_toggles, lag_grid = lag_grid)
grow_post  <- tidyPosterior(fit = grow_fit,  model_name = "Growth", process_name = "grow", 
                            best_fit = best_fit, clim_toggles = clim_toggles, lag_grid = lag_grid)

# Join all posteriors together
posts <- bind_rows(surv_post, rec_post, grow_post)

# Add missing covariate columns with NA values
beta_names <- c("beta_tseas", "beta_tvar", "beta_ttrend", "beta_pseas", "beta_pvar", "beta_ptrend")
missing_covars <- beta_names[!(beta_names %in% colnames(posts))]
for (var in missing_covars) {
  posts[[var]] <- NA_real_
}

# Format male population size posteriors
mf_post <- readRDS("Outputs/stan_fits/mf_pop.rds") %>%
  tidybayes::tidy_draws() %>%
  rename(alpha = b_Intercept,
         beta  = b_mna.female) %>%
  # Rename month random intercept parameters
  rename_with(., ~ paste0("alpha_month_", readr::parse_number(.)), .cols = starts_with("r_month")) %>%
  select(draw = .draw, sigma, contains(c("alpha", "beta"))) %>%
  mutate(process = "male_pop")

# Add male population size posterior to demographic posteriors
all_posts <- bind_rows(posts, mf_post)

# Save all nicely formatted posteriors as one .rds (MNA fit, 3 demographic model fits)
saveRDS(all_posts, "Outputs/stan_fits/neat_posteriors_alldat_28d.rds")
