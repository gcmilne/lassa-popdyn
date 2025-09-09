#------------------------------------------------#
# Simulating Lassa infection dynamics in Nigeria #
#------------------------------------------------#

# 0. Script purpose ----
# - Uses fitted demographic process models & transmission posteriors
# - Does posterior predictions of rodent infection dynamics in Nigerian states
# - Saves output as .rds files


# 1. Set up IPM for simulation (load packages, data, set parameters) ----
use_cluster <- TRUE  # Use the high performance cluster?
cluster <- "cropdiversity"
# cluster   <- "nhm"
if (use_cluster) {  # Set directory for HPC
  proj_dir <- ifelse(cluster == "cropdiversity", 
                     "/home/gmilne/projects/nhm/gmilne/mastomys_modelling", 
                     "/home/gregm/mastomys_modelling")
  setwd(proj_dir)
}
script_path <- ifelse(use_cluster, "scripts", "R")
source(paste0(script_path, "/setup_IPM.R"))

# Remove unnecessary files
rm(seroprev.dat, surv.adj, morogoro.clim)

# Account for uncertainty in IFA sensitivity & specificity during fitting?
use_se_sp <- FALSE

# Set state index (j) & iteration number (i)
if (use_cluster) {
  i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  j <- as.integer(Sys.getenv("STATE_INDEX"))
} else {
  i <- 1; j <- 2
}


# 2. Load data ----

# Nigeria climate data (env_formatting.R)
nigeria.clim <- readRDS(paste0(data_path, "/monthly.climate.1990-2025_separate_trend.rds")) %>%
  filter(state != "Morogoro") %>% 
  ungroup()

# Select state
selected_state <- unique(nigeria.clim$state)[j]

# Lassa confirmed case data (format_lassa_case_data.R)
formatted.case.data <- readRDS(paste0(data_path, "/combined_case_data.rds"))

# Posterior parameter samples for posterior predictions
post_samples <- readRDS(paste0(output_path, "/posterior_predict_samples_N=1000.rds"))
sampled_inf_post <- post_samples[i, ]  # subset iteration number
rm(post_samples)

# Generate climate lag grid to join temperature & precipitation lags to model fits
lag_seq  <- seq(0, 168, by = 28)
lag_grid <- setNames(expand.grid(lag_seq, lag_seq), c("t_lag", "p_lag")) %>% 
  as.data.frame() %>% 
  # Use variable naming syntax that matches syntax in posteriors
  mutate(lag_wgt_var = paste0("lag_weights[", row_number(), "]"))

# Generate dates at 28d intervals to extract climate data for
dates   <- getDateVec(start_date = min(formatted.case.data$date_mid) - years(5),  # Burn-in
                      end_sampling = max(formatted.case.data$date_mid),  # Run model past final case data 
                      dt = 28)
start_date <- min(dates)
end_date   <- max(dates)

# Create dataframe of lagged climatic data using function in R/03_lag_env.R
lag_clim <- getLagClim(dates = dates, 
                       clim_dat = nigeria.clim %>% filter(state %in% selected_state), 
                       lag_grid = lag_grid)
lag_clim$state <- selected_state

# Clear environment
rm(nigeria.clim, formatted.case.data)
gc()

# Read state-specific survival adjustment values & demographic parameter
demog.pars$survival_adj <- median(readRDS(paste0(output_path, "/", selected_state, "_survadj_N=1000.rds")))


# 3. Posterior predictions of rodent infection dynamics for Nigerian state ----

# Define function for posterior prediction
sim_nigeria <- function(iter, inf_post, demog_posts, inf_pars, lagged_clim) {
  
  # Get matching demographic posterior samples
  sampled_demog_posts <- demog_posts %>%
    filter(
      (process == "surv"     & draw == inf_post$surv.draw.id) |
        (process == "rec"      & draw == inf_post$rec.draw.id)  |
        (process == "grow"     & draw == inf_post$grow.draw.id) |
        (process == "male_pop" & draw == inf_post$malepop.draw.id)
    )
  
  rm(demog_posts)
  
  # Run IPM with climate-demographic & transmission inputs & store
  inf_pars$trans_rate     <- inf_post$trans_rate
  inf_pars$vertical_trans <- inf_post$vertical_trans
  
  mod <- RunRodentModel(model = sampled_demog_posts, env = lag_clim, cmr.data = cmr.data,
                        demog.pars = demog.pars, inf.pars = inf_pars,
                        n_init = n_init, start_date = start_date, end_date = end_date,
                        seed = 1, summarised = TRUE) %>%
    mutate(iter = iter, state = lagged_clim$state[1]) %>% 
    filter(date >= ymd("2017-12-29")) %>% # date before 1st case data (2018)
    reframe(date, state, iter, inf.adult, inf.sub,
            inf.preg = J_h + J_v,
            inf = inf.adult + inf.sub, 
            rel_inf = inf / max(inf)) %>% 
    distinct()
  
  return(mod)
}

# Simulate model for given state for given iteration number
sim <- sim_nigeria(iter = i, inf_post = sampled_inf_post, demog_posts = posts, inf_pars = inf.pars, lagged_clim = lag_clim)

# Save output
if (use_cluster) saveRDS(sim, file = paste0("out/nigeria_sims/", selected_state, "/", i, ".rds"))
