#---------------------------------------------------------------#
# Prepare fitted model outputs for input into mechanistic model #
#---------------------------------------------------------------#

# 0. Script purpose ----
# Defines functions to lag climate data for input into rodent integral projection model


# 1. Load packages ----
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)


# 2. Define functions ----

# Function to generate lags in climate data in 2 dimensions
getBothLags <- function(clim_dat, lags) {
  
  # clim_dat: Climate data with columns: date, t_seas, t_var, p_seas, p_var
  # lags: climate lags (in days) in order: (1) temperature, (2) precipitation
  
  # Generate lagged temp data
  clim_lag_t <- clim_dat %>% 
    reframe(date) %>% 
    mutate(lag_date_t = ymd(date) - lags[1]) %>% 
    left_join(clim_dat, by = c("lag_date_t" = "date")) %>% 
    dplyr::select(date, t_seas, t_var, t_trend)
  
  # Generate lagged precipitation data
  clim_lag_p <- clim_dat %>% 
    reframe(date) %>% 
    mutate(lag_date_p = ymd(date) - lags[2]) %>% 
    left_join(clim_dat, by = c("lag_date_p" = "date")) %>% 
    dplyr::select(date, p_seas, p_var, p_trend)
  
  # Join both climate data together
  both_lag_clim <- full_join(clim_lag_t, clim_lag_p, by = "date")
  
  # Return
  return(both_lag_clim)
  
}

# Function to generate lag array for all combination of lags for each climatic variable
getLagClim <- function(dates, clim_dat, lag_grid) {
  
  # dates: dates to keep climate data for (ideally in 28d intervals for simple input into model)
  # clim_dat: Climate data with columns: date, t_seas, t_var, p_seas, p_var
  # lag_grid: Grid of climate lags (in days) for temperature and precipitation

  lag_list <- vector("list", nrow(lag_grid))
  
  # Populate list for all variables & all lags
  for (i in 1:nrow(lag_grid)) {
    
    # Get lagged data & keep only relevant dates
    lagged_data <- getBothLags(clim_dat = clim_dat, lags = as.numeric(lag_grid[i, c("t_lag", "p_lag")])) %>% 
      filter(date %in% dates)
    
    # Store dataframe in list element
    lag_list[[i]] <- data.frame(
      date  = lagged_data$date,
      ttrend = lagged_data$t_trend,
      tseas = lagged_data$t_seas,
      tvar  = lagged_data$t_var,
      ptrend = lagged_data$p_trend,
      pseas = lagged_data$p_seas,
      pvar  = lagged_data$p_var,
      t_lag = lag_grid[i, "t_lag"],
      p_lag = lag_grid[i, "p_lag"]
    )
    
  }
  
  bound_lag_list <- bind_rows(lag_list)
  
  return(bound_lag_list)
  
}

# Function to create a date vector given a start date & target end date
getDateVec <- function(start_date, end_sampling, dt) {
  
  # start_date: date to start vector from, in ymd format
  # end_sampling: target date to end vector, in ymd format
  # dt: no. days between each date element (numeric input)
  
  # Create date vector given inputs
  date_vec <- seq(start_date, end_sampling, by = dt)
  
  # Append extra dates to end of vector if sampling end date not within vector range
  while (max(date_vec) < end_sampling) {
    date_vec <- c(date_vec, max(date_vec) + dt)
  }
  
  # Return date vector
  return(date_vec)
  
}
