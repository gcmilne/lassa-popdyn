#-------------------------------------#
# Process survival calibration output #
#-------------------------------------#

# 0. Script purpose ----
# - Read in multiple survival adjustments files calculated in calibrate_survival.R
# - Combine & save estimates for each state, each as a single .rds


# 1. Load packages ----
pacman::p_load(dplyr, purrr)


# 2. Process & save survival adjustments ----

# Read in climate data to get state names
states <- readRDS("Data/monthly.climate.1990-2025_separate_trend.rds") %>% 
  pull(state) %>% unique()

for (i in 1:length(states)) {
  
  # Get state
  selected_state <- states[i]
  
  # File path
  filepath <- paste0("Outputs/surv_adj/", selected_state, "/")
  
  if (!file.exists(filepath)) {
    warning(paste0("Files for ", selected_state, " not found"))
    next
    
  } else {
    
    # Read in & combine survival adjustment values for given state
    file_list <- list.files(path = paste0("Outputs/surv_adj/", selected_state, "/"),
                            pattern = "\\.rds$",
                            full.names = TRUE, 
                            ignore.case = TRUE)
    surv_adj <- file_list %>% map(readRDS) %>% unlist()
    
    # Save as one .rds
    saveRDS(surv_adj, paste0("Outputs/", selected_state, "_survadj_N=", length(surv_adj), ".rds"))
  }
}
