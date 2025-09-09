#-------------------------------------#
# Script to set up IPM for simulation #
#-------------------------------------#

# 0. Script purpose ----
# - Sets up environment needed to simulate the integral projection model (sourced 
#   by other scripts)
# - E.g. data, model parameters


# 1. Load packages ----
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(truncnorm)
library(matrixStats)


# 2. Source scripts ----

# Requires definition of use_cluster (FALSE or TRUE before sourcing of script)

# Set file paths
if (use_cluster) {
  script_path <- "scripts"
  data_path   <- output_path <- "data"
} else {
  # Set script & data paths
  script_path <- "R"
  data_path   <- "Data"
  output_path <- "Outputs"
}

# Source scripts
source(paste0(script_path, "/rodent_model.R"))  # rodent model
source(paste0(script_path, "/03_lag_env.R"))    # for creating lagged climatic data
source(paste0(script_path, "/fit_funcs.R"))     # fitting functions


# 3. Load data ----

# Morogoro climate data (unscaled)
morogoro.clim <- readRDS(paste0(data_path, "/monthly.climate.1990-2024_separate_trend.rds")) %>%
  filter(state == "Morogoro") %>% 
  ungroup() %>% 
  dplyr::select(date, contains(c("seas", "var", "trend")))

# Load data to infer priors for IFA sensitivity and specificity
sens_spec <- read.csv(paste0(data_path, "/IFA_se_sp.csv"))
shape2 <- 8 # set shape parameter 2 for deriving priors
set.seed(001)

# Set priors for IFA sensitivity
mean_sens <- weighted.mean(sens_spec$sensitivity, w = sens_spec$n_tot) # desired mean sensitivity
shape1_est <- mean_sens / (1 - mean_sens) * shape2 # calculate shape1 given set shape2 & mean se
sens_prior <- rbeta(1000, shape1 = shape1_est, shape2 = shape2)

# Set priors for IFA specificity
mean_spec <- weighted.mean(sens_spec$specificity, w = sens_spec$n_tot) # desired mean specificity
shape1_est <- mean_spec / (1 - mean_spec) * shape2  # calculate shape1 given set shape2 & mean sp
spec_prior <- rbeta(1000, shape1 = shape1_est, shape2 = shape2)

# Sensitivity & specificity priors (one row per "posterior" sample)
sero_posts <- data.frame(sens = sens_prior, spec = spec_prior)

# Load posteriors of demographic model fitting (unscaled)
posts <- readRDS(paste0(output_path, "/unscaled_posteriors.rds"))

# Load MOSA CMR data
cmr.data <- read.csv(paste0(data_path, "/MOSA_CMR_extracted_20240820.csv"), stringsAsFactors = FALSE) %>%
  # Convert date to correct format
  mutate(DATE = ymd(DATE)) %>%
  # Remove rows with NA date
  filter(!is.na(DATE))

# Load formatted MORV seroprevalence data
seroprev.dat <- readRDS(paste0(data_path, "/MORV_seroprev_formatted.rds")) %>% 
  # Leave out data after 2017 (use for validation)
  filter(year(date_mid) < 2017)


# 4. Format climate data ----

# Generate climate lag grid to join temperature & precipitation lags to model fits
lag_seq  <- seq(0, 168, by = 28)
lag_grid <- setNames(expand.grid(lag_seq, lag_seq), c("t_lag", "p_lag")) %>% 
  as.data.frame() %>% 
  # Use variable naming syntax that matches syntax in posteriors
  mutate(lag_wgt_var = paste0("lag_weights[", row_number(), "]"))

# Generate dates at 28d intervals to extract climate data for
dates      <- getDateVec(start_date = min(seroprev.dat$session.start) - years(5), # burn-in of 5 years
                         end_sampling = max(seroprev.dat$session.start),  # Run model past final serosurvey sampling date, 
                         dt = 28)
start_date <- min(dates)
end_date   <- max(dates)

# Create dataframe of lagged climatic data using function in R/03_lag_env.R
lag_clim <- getLagClim(dates = dates, clim_dat = morogoro.clim, lag_grid = lag_grid)


# 5. Set model parameters ----

# Estimated survival adjustments for Morogoro
surv.adj <- readRDS(paste0(output_path, "/Morogoro_survadj_N=1000.rds"))

# Demographic parameters
demog.pars <- list(
  # Set litter size of n females
  effective_litter = 5,
  # Set population size adjustment parameter
  pop_adj = 1,
  # Set survival adjustment parameter using median from estimated values
  survival_adj = median(surv.adj),
  # Define juvenile weight range
  juvenile_wgt_range = c(11, 13)
)

# Infection parameters
inf.pars <- data.frame(
  inf_prev = mean(seroprev.dat$seroprev),
  vertical_trans = 0.9,
  trans_rate = 0.001,
  mat_ab_toggle = 1,    # Maternal Ab toggled on
  dens_cont_toggle = 0  # no density-dependent transmission 
)

# Set initial population size
n_init <- 1000
