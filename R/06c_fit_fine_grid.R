#---------------------------------------#
# Fit fine grid of infection parameters #
#---------------------------------------#

# 0. Script purpose ----
# - Uses Morogoroarenavirus serooprevalence data & fine grid of transmission 
#   parameters estimated in previous script to derive fine grid estimates of 
#   horizontal & vertical transmission parameters for N demographic process 
#   model posterior samples
# - These fine grids are saved & subsequently used in posterior predictions


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


# 2. Define fine grid search workflow ----

# Load fine grid of parameter values to fit model for
grid_fine <- readRDS(paste0(output_path, "/fine_grid.rds"))

# Function to run finer grid search
run_fine_grid <- function(j, grid_fine) {
  
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
  
  # Set parameter values
  alpha <- grid_fine$vertical_trans
  beta  <- grid_fine$trans_rate
  
  # Run grid search
  grid_result <- GridSearch(alpha_vals = alpha, beta_vals = beta,
                            seroprev.dat = seroprev.dat, demog.posts = par.sam,
                            morogoro.clim = lag_clim, cmr.data = cmr.data,
                            demog.pars = demog.pars, inf.pars = inf.pars,
                            n_init = n_init, start_date = start_date, end_date = end_date,
                            seed = 1, 
                            sens = sampled_sens, spec = sampled_spec,
                            loglik_func = calcLogLik,
                            alpha_prior_func = prior_vertical_trans,
                            beta_prior_func = prior_trans_rate,
                            make_grid = FALSE)
  
  result <- bind_rows(grid_result) %>%
    mutate(post_iter = j, surv.draw.id = surv.draw.id, rec.draw.id = rec.draw.id,
           grow.draw.id = grow.draw.id, malepop.draw.id = malepop.draw.id,
           sero.draw.id = sero.draw.id) %>%
    dplyr::select(vertical_trans, trans_rate, joint_prior, joint_loglik, joint_posterior, 
                  i, post_iter, contains(".id")) %>%
    distinct()
  
  return(result)
  
}

# Get SLURM task ID
if (use_cluster) j <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) else j <- 1

# Run the grid search for this posterior sample
result <- run_fine_grid(j = j, grid_fine = grid_fine)

# Calculate posterior weights
result_weight <- result %>%
  mutate(weight = exp(joint_posterior - logSumExp(joint_posterior))) %>% 
  arrange(desc(weight))

# Save the result
if (use_cluster) saveRDS(result_weight, file = paste0("out/fine_grid/grid_", j, ".rds"))
