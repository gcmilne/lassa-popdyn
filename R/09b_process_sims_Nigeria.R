#--------------------
# Process simulations
#--------------------

# 0. Script purpose ----
# - Loads posterior predictions of Nigerian state rodent infection dynamics
# - Saves the N simulations for each state as an .rds file


# 1. Load packages ----
pacman::p_load(dplyr, purrr)


# 2. Read in & process Nigeria simulations ----

# Read in climate data to get state names
states <- readRDS("Data/monthly.climate.1990-2025_separate_trend.rds") %>% 
  pull(state) %>% unique()

for (i in 1:length(states)) {
  
  # Get state
  selected_state <- states[i]
  
  # File path
  filepath <- paste0("Outputs/nigeria_sims/", selected_state, "/")
  
  if (!file.exists(filepath)) {
    warning(paste0("Files for ", selected_state, " not found"))
    next
    
  } else {
    
    # Read in & combine simulations for given state
    file_list <- list.files(path = paste0("Outputs/nigeria_sims/", selected_state, "/"),
                            pattern = "\\.rds$",
                            full.names = TRUE, 
                            ignore.case = TRUE)
    sims <- file_list %>% map(readRDS) %>% bind_rows()
    
    # Save as one .rds
    saveRDS(sims, paste0("Outputs/", selected_state, "_sims_N=", n_distinct(sims$iter), ".rds"))
  }
}
