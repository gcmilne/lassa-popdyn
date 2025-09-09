#--------------------------------#
# Calibrate survival adjustments #
#--------------------------------#

# 0. Script purpose ----
# - Loads data: environmental, rodent trapping, demographic process model parameter estimates
# - Simulates integral projection model for N posterior simulations to calibrate 
#   survival adjustment for Morogoro & Nigerian states
# - Saves outputs as .rds files


# 1. Load packages ----
library(dplyr)
library(purrr)
library(lubridate)
library(stringr)
library(tidyr)
library(truncnorm)


# 2. Source scripts & set directory ----

# Set directory & file paths
use_cluster <- TRUE  # Use the high performance cluster?
cluster <- "cropdiversity"
# cluster   <- "nhm"

if (use_cluster) {
  proj_dir <- ifelse(cluster == "cropdiversity", 
                     "/home/gmilne/projects/nhm/gmilne/mastomys_modelling", 
                     "/home/gregm/mastomys_modelling")
  setwd(proj_dir)
  script_path <- "scripts"
  data_path   <- output_path <- "data"
} else {
  script_path <- "R"
  data_path   <- "Data"
  output_path <- "Outputs"
}

# Source mechanistic model & auxiliary functions 
source(paste0(script_path, "/rodent_model.R"))

# Source script for creating lagged climatic data
source(paste0(script_path, "/03_lag_env.R"))


# 3. Load data ----

# Rodent capture-mark-recapture data
cmr.data <- read.csv(paste0(data_path, "/MOSA_CMR_extracted_20240820.csv")) %>% 
  mutate(DATE = ymd(DATE))

# Model posterior estimates of climate-demographic parameters (unscaled)
post <- readRDS(paste0(output_path, "/unscaled_posteriors.rds"))

# Climate data (unscaled)
env_df_master <- readRDS(paste0(data_path, "/monthly.climate.1990-2025_separate_trend.rds")) %>%
  ungroup()


# 4. Define model parameters ----

# Generate climate lag grid to join temperature & precipitation lags to model fits
lag_seq  <- seq(0, 168, by = 28)
lag_grid <- setNames(expand.grid(lag_seq, lag_seq), c("t_lag", "p_lag")) %>% 
  as.data.frame() %>% 
  # Use variable naming syntax that matches syntax in posteriors
  mutate(lag_wgt_var = paste0("lag_weights[", row_number(), "]"))

# Generate dates at 28d intervals to extract climate data for
dates      <- getDateVec(start_date = min(cmr.data$DATE) - max(lag_grid$t_lag), # Start of CMR study
                         end_sampling = max(env_df_master$date),  # Last climate data
                         dt = 28)
start_date <- min(dates)
end_date   <- max(dates)

# IPM demographic parameters
demog.pars <- list(
  # Set litter size of n females
  effective_litter = 5,
  # Set population size adjustment parameter
  pop_adj = 1,
  # Set survival adjustment parameter
  survival_adj = 0.8,  # estimated in this script
  # Define juvenile weight range
  juvenile_wgt_range = c(11, 13)
)

# IPM infection parameters (irrelevant for population dynamics)
inf.pars <- list(
  # Set initial prevalence
  inf_prev = 0.4,
  # Set vertical transmission rate
  vertical_trans = 0.5,
  # Set transmission rate (beta)
  trans_rate = 0.001,
  # Toggle for maternal Abs (1 = on, 0 = off)
  mat_ab_toggle = 1,
  dens_cont_toggle = 0 
)

# Set initial population size
n_init <- 1000


# 5. Calibrate survival adjustment ----

# Survival adjustment = Probability of individual captured @ time t with weight w
#                       still being alive at t+1 despite never being recaptured

# Goal: to achieve mean pop. size of around n_init, using linear interpolation

# Function to calibrate survival adjustment
CalibrateSurvival <- function(lag_clim, stat.model, cmr.data,
                              demog.pars, inf.pars, n_init, start_date, end_date,
                              seed, cutoff_date, break.while.loop.time, sens = 0.01) {
  
  # lag_clim: lagged climate data
  # stat.model: Bayesian statistical model estimates of demographic process parameters
  # cmr.data: Rodent CMR data
  # demog.pars: IPM demographic parameters
  # inf.pars: IPM transmission parameters
  # n_init: initial simulated population size
  # start_date & end_date: start & end dates of simulation
  # seed: ensure repeatable random sampling
  # break.while.loop.time: time (in seconds) after which the while() loop should exit
  # sens: sensitivity when comparing simulated N vs. initial N

  # Set range of survival adjustment values
  x_val <- c(0.6, 1)
  
  # Initialise storage for state-specific differences btwn simulated mean N & initial N
  y_val <- rep(NA, length(x_val))
  
  # For given state, run model for each survival value
  sim <- list()
  for (i in seq_along(x_val)) {
    
    # Change survival value
    demog.pars$survival_adj <- x_val[i]
    
    # Run the model & store output
    sim[[i]] <- RunRodentModel(model      = stat.model,
                               env        = lag_clim,
                               cmr.data   = cmr.data,
                               demog.pars = demog.pars,
                               inf.pars   = inf.pars,
                               n_init     = n_init,
                               start_date = start_date,
                               end_date   = end_date,
                               seed       = 1) %>%
      # Ensure model output contains later dates than cut-off
      filter(date > as.Date(cutoff_date, origin = "1970-01-01"))
    
    # Calculate difference btwn mean simulated pop. size & initial pop. size
    y_val[i] <-  mean(sim[[i]]$mod.pop) - n_init
    
  }
  
  # Get value with least difference to initial N
  dif <- y_val[which(abs(y_val) == min(abs(y_val)))][1]
  
  # Set time before while() loop
  t.before.loop <- Sys.time()
  
  # Interpolate new survival adjustment value based on difference from initial N
  while (abs(dif) > sens * n_init) {
    
    # New interpolated survival adjustment value
    demog.pars$survival_adj <- x_val[1] - (x_val[2] - x_val[1]) * y_val[1] / (0.5 * y_val[2] - y_val[1])
    
    t.prerun <- Sys.time()  # store time to calculate model run time
    
    # Simulate model with new estimated survival adjustment value
    D <- RunRodentModel(model = stat.model, env = lag_clim, cmr.data = cmr.data,
                        demog.pars = demog.pars, inf.pars = inf.pars,
                        n_init = n_init, start_date = start_date, end_date = end_date,
                        seed = seed) %>%
      # Ensure model output contains later dates than cut-off
      filter(date > as.Date(cutoff_date, origin = "1970-01-01"))
    
    # Calculate difference between mean simulated population size & initial N
    dif <- mean(D$mod.pop) - n_init
    
    # Set new values based on initial N vs. simulated N difference
    if (dif < 0 & dif >= y_val[1]) {
      y_val[1] <- dif
      x_val[1] <- demog.pars$survival_adj
      
    } else if (dif > 0 & dif <= y_val[2]) {
      y_val[2] <- dif
      x_val[2] <- demog.pars$survival_adj
    }
    
    # Set current system time
    t.in.loop <- Sys.time()
    
    # Leave while loop if can't find appropriate survival adjustment given time
    if (as.numeric(t.in.loop - t.before.loop) > break.while.loop.time) break
  }
  
  # Store adjustment value that minimises this difference
  survival.est <- x_val[which.min(abs(y_val))]
  
  # Return the updated survival estimates
  return(survival.est)
  
}

# Set state index (j) & iteration number (i)
if (use_cluster) {
  i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  j <- as.integer(Sys.getenv("STATE_INDEX"))
} else {
  i <- 1; j <- 1
}

# Select state & get dates
selected_state <- unique(env_df_master$state)[j]
end_date <- if_else(selected_state == "Morogoro", max(cmr.data$DATE), max(env_df_master$date))
dates_cropped <- dates[dates <= end_date]

# Create dataframe of state-specific lagged climatic data using function in R/03_lag_env.R
lag_clim <- getLagClim(dates = dates_cropped, 
                       clim_dat = env_df_master %>% filter(state == selected_state,
                                                           date <= end_date) %>% 
                         dplyr::select(date, contains(c("seas", "var", "trend"))), 
                       lag_grid = lag_grid)

# Sample climate-demographic parameter values from posteriors
set.seed(i)
par.sam <- post %>%
  group_by(process) %>%
  mutate(sampled_draw = sample(unique(draw), size = 1)) %>%
  filter(draw == sampled_draw) %>%
  ungroup()

# Calibrate survival adjustment
surv_adj <- CalibrateSurvival(lag_clim = lag_clim,
                              stat.model =  par.sam,
                              cmr.data = cmr.data,
                              demog.pars = demog.pars,
                              inf.pars = inf.pars,
                              n_init = n_init,
                              start_date = start_date,
                              end_date = end_date,
                              seed = 1,
                              cutoff_date = start_date,
                              break.while.loop.time = 25)

# Save
if (use_cluster) saveRDS(surv_adj, file = paste0("out/surv_adj/", selected_state, "/", i, ".rds"))
