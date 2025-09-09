#--------------------------------------------#
# Extract precipitation and temperature data #
#--------------------------------------------#

# 0. Session details ----
# - Get Nigerian state centroids & Morogoro field site coordinates
# - Extract rainfall & temperature data from rasters at these centroids/coordinates
# - Format data to account for NAs & erroneous values (eg sub-zero temps in Nigeria)
# - Calculate rolling averages at same temporal resolution as rodent population model
# - Decompose climatic variability into seasonal & variable components
# - Make plots of outputs
# - Save decomposed climate variables for analyses in other scripts


# 1. Load libraries ----
pacman::p_load(dplyr, ggplot2, terra, sf, ncdf4, raster, stringr, purrr, tidyr,
               lubridate, patchwork, zoo)

# Set toggle (0/1; off/on) for plotting throughout the script 
doplot <- 0


# 2. Define functions ----

# Function to extract climatic value at specific coordinates (state centroids)
get.daily.values  <- function(lon, lat, raster) {
  
  # lon: longitude coordinates
  # lat: latitude coordinates
  # raster: SpatRaster object as loaded by terra::rast()
  
  # Create dataframe of lon/lat coordinates
  points <- data.frame(lon = lon, lat = lat)
  
  # Convert points to spatial object
  coords <- vect(points, crs = crs(raster))
  
  # Extract all time steps at once (terra::extract to avoid use of tidyr::extract)
  daily.values <- terra::extract(raster, coords, method = "bilinear") %>% 
    # Remove ID column
    dplyr::select(!(contains("ID"))) %>% 
    # Return as numeric vector
    as.numeric()
  
  # Handle NA values in raster data
  if (any(is.na(daily.values))) {
    
    # Indices of NA values in daily values vector
    na.positions <- which(is.na(daily.values))
    
    # Print informative warning message giving state & date
    for (i in seq_along(na.positions)) {
      warning(paste0("NA detected for ", clim.df$state[clim.df$lon == lon], 
                     " on ", time(raster[[na.positions[i]]]), 
                     "; performing linear interpolation to replace missing values"))
    }
    
    # Replace NAs with values using linear interpolation
    daily.values[na.positions] <- zoo::na.approx(daily.values, na.rm = FALSE)[na.positions]
    
  }
  
  # Return daily values
  return(daily.values)
}


# Function to correct anomalous daily climate data by comparing to long-term mean
correct.daily.anomalies <- function(daily.values, lt.mean, threshold.prop) {
  
  # daily.values: daily values across given year
  # lt.mean: long-term mean
  # threshold.prop: proportional difference to long-term mean that classifies value as anomalous
  
  # Check if daily.values are from a leap year
  if (length(daily.values) == 366) {

    # Interpolate long-term mean values for Feb 29th to make vector lengths equal
    lt.mean <- c(lt.mean[1:59], NA, lt.mean[60:length(lt.mean)])
    lt.mean[is.na(lt.mean)] <- zoo::na.approx(lt.mean, na.rm = FALSE)[is.na(lt.mean)]

    # If data for <365 days (e.g. data for ongoing year)
  } else if (length(daily.values) < 365) {

    # Chop long-term mean vector to be same length
    lt.mean <- lt.mean[1:length(daily.values)]

  }
  
  # Calculate relative threshold as proportion of long-term mean
  rel.threshold <- lt.mean * threshold.prop
  
  # Calculate absolute difference between each daily value and matched LT mean value
  daily.diff <- abs(daily.values - lt.mean)
  
  # Get indices of anomalous daily values where difference from LT mean > threshold
  anom.boolean <- daily.diff > rel.threshold
  
  # Replace anomalous values with NA
  daily.values[anom.boolean] <- NA
  
  # Estimate NA value(s) via linear interpolation
  daily.values[anom.boolean] <- zoo::na.approx(daily.values, na.rm = FALSE)[anom.boolean]
  
  # Return a list containing corrected values and logical indicating whether data are interpolated
  result <- list(values = daily.values, interpolated = anom.boolean)
  
  return(result)  
  
}


# Function to append NAs to a vector such that its final size equals the year's
pad.na <- function(vec, year) {
  
  # vec: A vector of any class
  # year: The year, as a character or numeric input
  
  # Is the year a leap year?
  leap.year <- lubridate::leap_year(as.numeric(year))
  
  # Calculate no. days in year
  days.in.year <- ifelse(test = leap.year, yes = 366, no = 365)
  
  # Calculate no. NAs to add
  num.na <- days.in.year - length(vec)
  
  # Append NAs to vector & return
  vec <- append(vec, rep(NA, num.na))
  
  return(vec)
}


# Function to format decomposed time series to return dataframe
FormatDecomposedTimeSeries <- function(decomposed_ts, year, day) {
  formatted_data <- data.frame(year  = year,
                               day   = day,
                               seas  = decomposed_ts$time.series[, "seasonal"], 
                               var   = decomposed_ts$time.series[, "remainder"], 
                               trend = decomposed_ts$time.series[, "trend"])
  return(formatted_data)
}


# 3. Define coordinates for raster data extraction ----

# Define lon/lat coordinates for Morogoro field site
clim.df <- data.frame(
  country = "Tanzania",
  state   = "Morogoro",
  lon     = 37.6333, # 37°38′E = 37 + (38/60) = 37.6333 (since East is positive)
  lat     = -6.85    # 6°51′S = -6 - (51/60) = -6.85 (since South is negative)
)

# Load Nigeria state shape files
states_sf <- st_as_sf(  # Convert shapefile to sf object
  rgeoboundaries::geoboundaries(  # NB requires internet connection
    country = "Nigeria",
    adm_lvl = "adm1",  # Get state borders 
    type = NULL,
    quiet = TRUE
  )
)

# Specify Nigerian states of interest
states.vec <- c("Edo", "Ondo", "Bauchi", "Taraba", "Ebonyi")

# Format shape files
states_sf <- states_sf %>% 
  # Define binary variable for states of interest
  mutate(selected.state = as.factor(ifelse(shapeName %in% states.vec, 1, 0))) %>% 
  # Calculate state centroids
  mutate(state.centroid = st_centroid(geometry)) 

# Plot Nigeria on map with relevant states coloured & state centroids plotted
if (doplot) {
  states_sf %>% 
    ggplot() + 
    geom_sf(aes(fill = selected.state)) + 
    geom_sf(aes(geometry = state.centroid), color = "white", size = 2) +
    scale_fill_viridis_d() + 
    theme_minimal()
}

# Create dataframe of Nigerian state names & centroids
clim.nigeria.df <- states_sf %>% 
  filter(selected.state == 1) %>% 
  reframe(
    country = "Nigeria",
    state = shapeName,
    # Put centroids in long/lat format
    lon = st_coordinates(state.centroid)[ ,1],
    lat = st_coordinates(state.centroid)[ ,2]
  ) %>% 
  # Arrange alphabetically by state
  arrange(state)

# Combine Morogoro & Nigerian climate data into one dataframe
clim.df <- bind_rows(clim.nigeria.df, clim.df)


# 4. Load & stack daily climate rasters (1990-2025) ----

# Stack daily minimum temperature rasters (https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/)
tmin.raster.filepaths <- list.files("Data/rasters/min_temp", pattern = ".nc", full.names = TRUE)

# Stack daily maximum temperature rasters (https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/)
tmax.raster.filepaths <- list.files("Data/rasters/max_temp", pattern = ".nc", full.names = TRUE)

# Stack daily precipitation rasters (https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/)
precip.raster.filepaths <- list.files("Data/rasters/rainfall", pattern = ".nc", full.names = TRUE)

# Read all the raster files into separate lists
tmin.rasters   <- lapply(tmin.raster.filepaths, rast)
tmax.rasters   <- lapply(tmax.raster.filepaths, rast)
precip.rasters <- lapply(precip.raster.filepaths, rast)

# Visually check that all rasters contain same sampling years
sapply(X = list(tmin.raster.filepaths, tmax.raster.filepaths, precip.raster.filepaths),
       FUN = str_match_all, pattern = "[0-9]+")

# If so, name raster list elements with sampling year, using file path names
names(tmin.rasters) <- names(tmax.rasters) <- names(precip.rasters) <- precip.raster.filepaths %>% 
  str_match_all(., "[0-9]+") %>% 
  as.numeric()
# data.frame(names(precip.rasters), precip.raster.filepaths)  # check years match


# 5. Load & format long-term daily mean rasters (1991-2020) ----

# Load min temp (https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/)
tmin.ltmean <- rast("Data/rasters/lt_means/tmin.day.ltm.1991-2020.nc")

# Load max temp (https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/)
tmax.ltmean <- rast("Data/rasters/lt_means/tmax.day.ltm.1991-2020.nc")

# Load precipitation (https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/)
precip.ltmean <- rast("Data/rasters/lt_means/precip.day.1991-2020.ltm.nc")

# Extract all daily long-term means at Nigerian state centroids & Morogoro site
tmin.ltmean.list   <- pmap(clim.df, function(lon, lat, ...) get.daily.values(lon, lat, tmin.ltmean))
tmax.ltmean.list   <- pmap(clim.df, function(lon, lat, ...) get.daily.values(lon, lat, tmax.ltmean))
precip.ltmean.list <- pmap(clim.df, function(lon, lat, ...) get.daily.values(lon, lat, precip.ltmean))

# Put in list
ltmeans <- setNames(list(tmin.ltmean.list, tmax.ltmean.list, precip.ltmean.list), c("tmin", "tmax", "precip"))


# 6. Extract daily climate data for Nigeria & Morogoro ----

# Initialise list storage for climate data 
clim.list <- vector("list", length(precip.rasters))

# For each year, extract daily climatic values for all coordinates (Nigeria & Morogoro)
for (i in seq_along(clim.list)) {
  
  # Store current year of data extraction
  current.year <- names(precip.rasters)[i]
  
  # Print year 
  print(paste("Extracting climate data for", current.year))
  
  # Extract daily climate values at all coordinates (Nigeria states & Morogoro)
  precip.values <- pmap(clim.df, function(lon, lat, ...) get.daily.values(lon, lat, precip.rasters[[i]]))
  tmin.values   <- pmap(clim.df, function(lon, lat, ...) get.daily.values(lon, lat, tmin.rasters[[i]]))
  tmax.values   <- pmap(clim.df, function(lon, lat, ...) get.daily.values(lon, lat, tmax.rasters[[i]]))
  
  # Find & correct anomalous min temp data (by comparing to long-term mean)
  tmin.corrected <- pmap(.l = list(tmin.values, ltmeans[["tmin"]]),
                         .f = ~ correct.daily.anomalies(.x, .y, threshold.prop = 0.6))
  
  # Get daily values (values interpolated if found to be anomalous)
  tmin.values <- map(tmin.corrected, "values")
  
  # Logical flag indicating if interpolation (T or F)
  tmin.interpolated <- map(tmin.corrected, "interpolated")
  
  # Find & correct anomalous max temp data (by comparing to long-term mean)
  tmax.corrected <- pmap(.l = list(tmax.values, ltmeans[["tmax"]]),
                         .f = ~ correct.daily.anomalies(.x, .y, threshold.prop = 0.6))
  
  # Get daily values (values interpolated if found to be anomalous)
  tmax.values <- map(tmax.corrected, "values")
  
  # Logical flag indicating if interpolation (T or F)
  tmax.interpolated <- map(tmax.corrected, "interpolated")
  
  # For 2025, append NAs to missing data
  if (current.year == "2025") {
    precip.values <- lapply(precip.values, pad.na, year = current.year)
    tmin.values   <- lapply(tmin.values, pad.na, year = current.year)
    tmax.values   <- lapply(tmax.values, pad.na, year = current.year)
    tmin.interpolated <- lapply(tmin.interpolated, pad.na, year = current.year)
    tmax.interpolated <- lapply(tmax.interpolated, pad.na, year = current.year)
  }
  
  # Store daily climate data output in dataframe (1 list element/year)
  clim.list[[i]] <- clim.df %>%
    mutate(
      # Extract year of sampling (should be the same across all rasters)
      year = current.year,
      # Add data as list-columns, using I() to preserve structure
      daily.precip = I(precip.values), 
      daily.tmin   = I(tmin.values),
      daily.tmax   = I(tmax.values),
      tmin.interpolated = I(tmin.interpolated),
      tmax.interpolated = I(tmax.interpolated)
    )
  
  # Update progress
  print(paste(i, "year(s) of", length(clim.list), "completed"))
  
}

# Combine listed outputs into dataframes
daily.clim.df <- bind_rows(clim.list)

# Save dataframe if file doesn't already exist
filepath <- paste0("Data/daily.climate.data.", min(daily.clim.df$year), "-", 
                   max(daily.clim.df$year), ".rds")

if (!file.exists(filepath)) saveRDS(daily.clim.df, filepath)


# 7. Format extracted daily climate data ----

daily.clim.df <- readRDS("Data/daily.climate.data.1990-2025.rds")

# Get non leap years
years <- unique(daily.clim.df$year)
non_leapyears <- years[!leap_year(as.numeric(years))]

# Create imaginary Feb 29th for these years
feb29 <- daily.clim.df %>% 
  filter(year %in% non_leapyears) %>% 
  distinct(state, country, lon, lat, year) %>% 
  mutate(day = 60)

# Interpolate Feb 29th on non-leap years
daily.clim.adjusted <- daily.clim.df %>% 
  # Unnest lists into rows and columns
  tidyr::unnest(cols = c(daily.precip, daily.tmin, daily.tmax,
                         tmin.interpolated, tmax.interpolated)) %>% 
  mutate(is.leapyear = leap_year(as.numeric(year))) %>%
  # Get day in year
  group_by(state, year) %>%
  mutate(day = row_number(),
         day = case_when(
           !is.leapyear & day < 60  ~ day,
           !is.leapyear & day >= 60 ~ day + 1,
           .default = day)) %>% 
  ungroup() %>% 
  # Join imaginary Feb 29th days & interpolate daily climate data
  full_join(feb29) %>% 
  arrange(state, year, day) %>% 
  fill(country, state, lon, lat, is.leapyear) %>% 
  group_by(state) %>% 
  mutate(daily.precip = zoo::na.approx(daily.precip, na.rm = FALSE),
         daily.tmin = zoo::na.approx(daily.tmin, na.rm = FALSE),
         daily.tmax = zoo::na.approx(daily.tmax, na.rm = FALSE))

# Put daily climate data in long format (one row per day per state) & format
daily.clim.long <- daily.clim.adjusted %>%
  group_by(state) %>% 
  mutate(
    daily.temp = (daily.tmin + daily.tmax) / 2,
    # Calculate rolling means with k day windows (right-aligned so average over k previous values, not using values from future)
    roll.avg.temp   = zoo::rollmean(daily.temp,   k = 28, fill = NA, align = "right"),
    roll.avg.precip = zoo::rollmean(daily.precip, k = 28, fill = NA, align = "right")
  ) %>%
  ungroup() %>% 
  # Remove dates with any missing daily & rolling average climate data 
  filter(!is.na(daily.temp) & !is.na(daily.precip) & !is.na(roll.avg.temp) & !is.na(roll.avg.precip))


# 8. Decompose climate data into trend, seasonal & variable components ----

# Define all climate variables to decompose
all.climvars <- c("roll.avg.precip", "roll.avg.temp")  # Use rolling averages to smooth

all.states.decomp <- NULL

for (j in seq_along(unique(daily.clim.long$state))) {
  for (i in seq_along(all.climvars)) {
    
    # Get state-specific data
    daily.clim.long.filtered <- daily.clim.long %>% filter(state %in% unique(daily.clim.long$state)[j])
    
    # Create time series object
    time.series <- ts(coredata(daily.clim.long.filtered[[all.climvars[i]]]), frequency = 366)
    
    # Decompose into seasonal, trend & remainder components
    decomp <- stl(time.series, s.window = "periodic", t.window = 366*5, robust = TRUE)
    
    # Format decomposed data
    decomp.clean <- FormatDecomposedTimeSeries(decomposed_ts = decomp, year = daily.clim.long.filtered$year, 
                                               day = daily.clim.long.filtered$day)
    
    
    decomp.clean$clim <- all.climvars[i]  # Label with indexed climatic variable
    decomp.clean$state <- daily.clim.long.filtered$state # Label with state
    
    # Combine as dataframe
    all.states.decomp <- bind_rows(all.states.decomp, decomp.clean)
    
  }
}

# Decompose into components of trend, seasonality & variability
seasvar <- all.states.decomp %>%
  mutate(
    # Creating separate identifiers for seas, var & trend by clim
    t_seas  = ifelse(str_detect(clim, "temp"), seas, NA),
    t_var   = ifelse(str_detect(clim, "temp"), var, NA),
    t_trend = ifelse(str_detect(clim, "temp"), trend, NA),
    p_seas  = ifelse(str_detect(clim, "precip"), seas, NA),
    p_var   = ifelse(str_detect(clim, "precip"), var, NA),
    p_trend = ifelse(str_detect(clim, "precip"), trend, NA)
  ) %>%
  dplyr::select(!c(seas, var, clim, trend)) %>%  # Drop redundant columns
  distinct() %>% 
  group_by(state, year, day) %>%
  # Collapse rows to complete transition to wide format
  summarize(
    t_seas = max(t_seas, na.rm = TRUE),
    t_var = max(t_var, na.rm = TRUE),
    p_seas = max(p_seas, na.rm = TRUE),
    p_var = max(p_var, na.rm = TRUE),
    t_trend = max(t_trend, na.rm = TRUE),
    p_trend = max(p_trend, na.rm = TRUE)
  ) %>%
  ungroup()

if (doplot) {
  patchwork::wrap_plots(
    seasvar %>%
      ggplot(aes(day, t_seas, group = year)) +
      geom_line(aes(y = t_seas + t_var), alpha = 0.4, linewidth = .2, col = "indianred") +
      geom_line() +
      scale_colour_viridis_c() +
      labs(x = "Day of year", y = "Detrended temperature (°C)") +
      theme_minimal() + 
      facet_wrap(~state),
    seasvar %>%
      ggplot(aes(day, p_seas, group = year)) +
      geom_line(aes(y = p_seas + p_var), alpha = 0.4, linewidth = .2, col = "lightblue") +
      geom_line() +
      labs(x = "Day of year", y = "Detrended precipitaion (mm)") +
      theme_minimal() + 
      facet_wrap(~state), 
    nrow = 2
  )
  ggsave("Plots/climatic_data/seasvar.png", dpi = 1200, 
         width = 6, height = 6, units = "in", bg= "white")
}

# Remove imaginary Feb 29th & calculate date
seasvar_formatted <- seasvar %>% 
  anti_join(feb29) %>% 
  mutate(is.leapyear = leap_year(as.numeric(year))) %>% 
  group_by(year) %>%
  mutate(day = case_when(
    !is.leapyear & day < 60  ~ day,
    !is.leapyear & day >= 60 ~ day - 1,
    .default = day),
    date = ymd(paste0(year, "-01-01")) + days(day - 1))

if (doplot) {
  patchwork::wrap_plots(
    seasvar_formatted %>% ggplot(aes(date, t_trend)) + geom_line() + 
      labs(x = "Date", y = "Temperature trend (°C)") + 
      theme_minimal() +
      facet_wrap(~state),
    seasvar_formatted %>% ggplot(aes(date, p_trend)) + geom_line() + 
      labs(x = "Date", y = "Precipitaion trend (mm)") + 
      theme_minimal() +
      facet_wrap(~state),
    nrow = 2)
  ggsave("Plots/climatic_data/trends.png", dpi = 1200, 
         width = 6, height = 6, units = "in", bg= "white")
}


# 9. Save decomposed climatic variables ----

# Save climate data
saveRDS(seasvar_formatted, paste0("Data/monthly.climate.", min(daily.clim.df$year), "-", max(daily.clim.df$year), 
                                  "_separate_trend.rds"))

# 10. Make plots ----

# Define shading regions for the seasons
season_shading <- data.frame(
  xmin = c(60, 244),  # Example: March 1 and Sept 1 (yday)
  xmax = c(151, 334), # May 31 and Nov 30 (yday)
  ymin = -Inf,
  ymax = Inf
)

# Plot for temperature
temp_plot <- seasvar %>%
  ggplot(aes(x = day)) +
  # Shaded seasons
  geom_rect(
    data = season_shading,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "gray90"
  ) +
  # faint lines: seasonal + variable per year
  geom_line(aes(y = t_seas + t_var, group = year),
            alpha = 0.2, colour = "indianred", linewidth = 0.2) +
  # Bold seasonal mean line
  geom_line(aes(y = t_seas), colour = "indianred", linewidth = 0.8) +
  labs(x = "Day of year", y = "Detrended temperature (°C)", col = "Component") +
  facet_wrap(~state) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.35, "cm"),
    text = element_text(family = "serif", size = 10),
    strip.text.x = element_text()
  )

# Plot for precipitation
precip_plot <- seasvar %>%
  ggplot(aes(x = day)) +
  # Shaded seasons
  geom_rect(
    data = season_shading,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "gray90"
  ) +
  # faint lines: seasonal + variable per year
  geom_line(aes(y = p_seas + p_var, group = year),
            alpha = 0.2, colour = "steelblue", linewidth = 0.2) +
  # Bold seasonal mean line
  geom_line(aes(y = p_seas), colour = "steelblue", linewidth = 0.8) +
  labs(x = "Day of year", y = "Detrended precipitation (mm)", col = "Component") +
  facet_wrap(~state) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.35, "cm"),
    text = element_text(family = "serif", size = 10),
    strip.text.x = element_text()
  )

# Combine with patchwork
if (doplot) {
  patchwork::wrap_plots(
    temp_plot,
    precip_plot,
    nrow = 2
  )
  
  # Save plot
  ggsave("Plots/climatic_data/seasvar.all.states.pdf",
         width = 8, height = 8, units = "in", device = cairo_pdf)
  ggsave("Plots/climatic_data/seasvar.all.states.png",
         width = 8, height = 8, units = "in", bg = "white", dpi = 1200)
}
