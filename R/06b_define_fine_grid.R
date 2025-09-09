#----------------------------------------------------------#
# Process coarse grid transmission fitting output from HPC #
#----------------------------------------------------------#

# 0. Script purpose ----
# - Load coarse grid fits .rds files
# - Combine & only keep those within a given probability threshold
# - Define new fine grid over feasible parameter space
# - Save fine grid for use in subsequent script to estimate approximate transmission posteriors


# 1. Load packages ----
library(tidyverse)


# 2. Load & combine coarse grid estimates ----

# List file paths of all .rds files in the directory (case-insensitive)
file_list <- list.files(path = "Outputs/model_fits/coarse_grid/", 
                        pattern = "\\.rds$", 
                        full.names = TRUE, 
                        ignore.case = TRUE)

# Read & combine .rds files
coarse_grid <- file_list %>%
  map(readRDS) %>%
  bind_rows()

# Filter to get parameter sets within probability threshold
coarse_subset <- coarse_grid %>% filter(within_bound)


# 3. Generate finer grid based on the noisy coarse grid search ----

# Set grid parameter bounds used in coarse grid
alpha_min <- 0.01
alpha_max <- 0.99
beta_min  <- 0.00001
beta_max  <- 0.0015

# Define no. grid points used for coarse grid
n_grid <- 20

# Define max permitted noise (ensures bounds respected)
max_alpha_noise <- (alpha_max - alpha_min) / (n_grid) / 2
max_beta_noise  <- (beta_max - beta_min)  / (n_grid) / 2

# Define coarse grid, leaving space for noise at both ends
alpha_base <- seq(alpha_min + max_alpha_noise, alpha_max - max_alpha_noise, length.out = n_grid)
beta_base  <- seq(beta_min  + max_beta_noise, beta_max  - max_beta_noise,  length.out = n_grid)

# Finding a bounding box for the feasible parameter space to make fine grid generation computationally cheaper
realised_alpha_min <- min(coarse_subset$vertical_trans) - max_alpha_noise
realised_alpha_max <- max(coarse_subset$vertical_trans) + max_alpha_noise
realised_beta_min <- min(coarse_subset$trans_rate) - max_beta_noise
realised_beta_max <- max(coarse_subset$trans_rate) + max_beta_noise

# Find elements of the base grid in the feasible region
base_subset <- coarse_subset %>%
  group_by(vertical_trans, trans_rate) %>%
  mutate(alpha_base = alpha_base[which.min(abs(vertical_trans-alpha_base))],
         beta_base = beta_base[which.min(abs(trans_rate-beta_base))]) %>% # pull the subsetted values to their "base" value on the coarse grid
  ungroup() %>%
  dplyr::select(alpha_base, beta_base) %>%
  unique() # remove duplicates

# Get coarse grid step sizes for both axes
tmp <- base_subset %>% arrange(alpha_base)
alpha_step <- diff(unique(tmp$alpha_base))[1]
tmp <- base_subset %>% arrange(beta_base)
beta_step  <- diff(unique(tmp$beta_base))[1]
rm(tmp)

# Define how much finer to make fine grid vs. coarse grid
n_finer <- 8  # applies to both axes so overall increase in no. param sets = n_finer^2

# Create a fine grid by uniformly subdividing each retained coarse grid cell
grid_finer_list <- purrr::pmap(base_subset, function(alpha_base, beta_base) {
  
  alpha_seq <- seq(alpha_base - alpha_step / 2 + alpha_step / (2 * n_finer),
                   alpha_base + alpha_step / 2 - alpha_step / (2 * n_finer),
                   length.out = n_finer)
  
  beta_seq <- seq(beta_base - beta_step / 2 + beta_step / (2 * n_finer),
                  beta_base + beta_step / 2 - beta_step / (2 * n_finer),
                  length.out = n_finer)
  
  expand.grid(vertical_trans = alpha_seq, trans_rate = beta_seq)
})

# Combine all grid cells to generate full fine grid
grid_finer <- bind_rows(grid_finer_list) 
sqrt(nrow(grid_finer))

# Save fine grid for use in next HPC script
saveRDS(grid_finer, "Outputs/fine_grid.rds")
