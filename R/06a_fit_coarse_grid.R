#-----------------------------------------#
# Fit coarse grid of infection parameters #
#-----------------------------------------#

# 0. Script purpose ----
# - Uses Morogoroarenavirus serooprevalence data to derive coarse grid estimates
#   of horizontal & vertical transmission parameters for N demographic process 
#   model posterior samples
# - These coarse grids are saved & subsequently used to derive finer posterior 
#    distributions in the subsequent script


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

# Account for uncertainty in IFA sensitivity & specificity during fitting?
use_se_sp <- FALSE


# 2. Estimate transmission parameters on coarse grid ----

## Define coarse grid search ----

# Set probability threshold for keeping posterior samples
prob_thresh <- 1 - (1/1e7)

# Set grid parameter bounds
alpha_min <- 0.01
alpha_max <- 0.99
beta_min  <- 0.00001
beta_max  <- 0.0015

# Define no. grid points for coarse grid
n_grid <- 20

# Define max permitted noise (ensures bounds respected)
max_alpha_noise <- (alpha_max - alpha_min) / (n_grid) / 2
max_beta_noise  <- (beta_max - beta_min)  / (n_grid) / 2

# Define coarse grid, leaving space for noise at both ends
alpha_base <- seq(alpha_min + max_alpha_noise, alpha_max - max_alpha_noise, length.out = n_grid)
beta_base  <- seq(beta_min  + max_beta_noise, beta_max  - max_beta_noise,  length.out = n_grid)
# expand.grid(alpha_base, beta_base) %>% ggplot(aes(Var1, Var2)) + geom_tile(col="white")

# Function to run grid search for given posterior sample
run_grid_for_posterior <- function(j) {
  
  set.seed(j)
  
  # Sample posterior draw for demographic processes
  par.sam <- posts %>%
    group_by(process) %>%
    mutate(sampled_draw = sample(unique(draw), size = 1)) %>%
    filter(draw == sampled_draw) %>%
    ungroup()
  
  # Sample sensitivity and specificity value from priors
  if (use_se_sp) {
    sero.draw.id <- sample(1:nrow(sero_posts), size = 1)
    sampled_sens <- sero_posts$sens[sero.draw.id]
    sampled_spec <- sero_posts$spec[sero.draw.id]
  } else {
    sero.draw.id <- NA
    sampled_sens <- 1
    sampled_spec <- 1
  }
  
  # Store sample IDs for demographic posteriors
  surv.draw.id    <- unique(par.sam$draw[par.sam$process == "surv"])
  rec.draw.id     <- unique(par.sam$draw[par.sam$process == "rec"])
  grow.draw.id    <- unique(par.sam$draw[par.sam$process == "grow"])
  malepop.draw.id <- unique(par.sam$draw[par.sam$process == "male_pop"])
  
  # Add noise
  alpha_noise <- runif(1, min = -max_alpha_noise, max = max_alpha_noise)
  beta_noise  <- runif(1, min = -max_beta_noise,  max = max_beta_noise)
  
  alpha_new <- alpha_base + alpha_noise
  beta_new  <- beta_base  + beta_noise
  
  # Run grid search
  grid_result <- GridSearch(alpha_vals = alpha_new, beta_vals = beta_new,
                            seroprev.dat = seroprev.dat, demog.posts = par.sam,
                            morogoro.clim = lag_clim, cmr.data = cmr.data,
                            demog.pars = demog.pars, inf.pars = inf.pars,
                            n_init = n_init, start_date = start_date, end_date = end_date,
                            seed = 1, 
                            sens = sampled_sens, spec = sampled_spec,
                            loglik_func = calcLogLik,
                            alpha_prior_func = prior_vertical_trans,
                            beta_prior_func = prior_trans_rate,
                            make_grid = TRUE)
  
  bind_rows(grid_result) %>%
    mutate(post_iter = j, surv.draw.id = surv.draw.id, rec.draw.id = rec.draw.id,
           grow.draw.id = grow.draw.id, malepop.draw.id = malepop.draw.id,
           sero.draw.id = sero.draw.id) %>%
    dplyr::select(vertical_trans, trans_rate, joint_prior, joint_loglik, joint_posterior,
                  i, post_iter, contains(".id")) %>%
    distinct()
}

# Set iteration number
if (use_cluster) j <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) else j <- 1

# Run the grid search for this posterior sample
result <- run_grid_for_posterior(j)

# Calculate which parts of coarse grid to keep based on weight threshold
result_keep <- result %>%
  mutate(weight = exp(joint_posterior - logSumExp(joint_posterior))) %>% 
  arrange(desc(weight)) %>% 
  mutate(cumsum_weight = cumsum(weight), 
         # captures edge parameters that cause cumulative sum to cross prob. threshold
         within_bound = cumsum_weight - weight < prob_thresh)

# Save the result
if (use_cluster) saveRDS(result_keep, file = paste0("out/coarse_grid/grid_", j, ".rds"))
