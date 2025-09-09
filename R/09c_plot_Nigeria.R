#------------------------------------#
# Plot Nigeria posterior predictions #
#------------------------------------#

# 0. Script purpose ----
# - Load Nigeria simulations and plot posterior predictions of infection dynamics
#   for each state


# 1. Load packages ----
pacman::p_load(dplyr, ggplot2, tidyr, foreach, patchwork, lubridate)

# Source custom ggplot theme
source("R/ggplot_theme.R")

dodge_width <- 0.8 # set aesthetics


# 2. Load data ----

# Lassa confirmed case data (04_format_lassa_case_data.R)
formatted.case.data <- readRDS("Data/combined_case_data.rds")

# Load MCMC chain estimates of incubation period (remove burn-in iterations <10k)
chain_incub <- read.csv("Data/incubation_chain_lognormal.csv")[10000:100000, ] %>% 
  # Get parameter values: mean & SD of log-normal distribution
  dplyr::select(meanlog, sdlog)

# Read in climate data to get state names
states <- readRDS("Data/monthly.climate.1990-2025_separate_trend.rds") %>% 
  filter(state != "Morogoro") %>% pull(state) %>% unique()


# 3. Define functions ----

# Function to generate poisson-distributed case data from confirmed cases
GetPoisCases <- function(cases, seed) {
  
  # Set seed to ensure repeatable sampling
  set.seed(seed)
  
  # Generate & return poisson-distributed cases (assuming confirmed cases are the means)
  pois.cases <- rpois(n = length(cases), lambda = cases)
  return(pois.cases)
  
}


# Function to generate a random integer vector (for missing case data)
GenRandomIntVec <- function(n, k, seed) {
  
  # seed: to ensure repeatability
  # n: number of dates for which cases are missing
  # k: total number of missing confirmed cases
  
  # Set seed to ensure repeatable sampling
  set.seed(seed)
  
  # Generate n random numbers between 0 and k
  random_values <- sample(0:k, n, replace = TRUE)
  
  # Edge case where values are all zeros - assign 1 to a random position
  if (sum(random_values) == 0) {
    scaled_values <- rep(0, n)
    scaled_values[sample(1:n, 1)] <- k
    
  } else {
    
    # Scale the vector so that the sum is exactly k
    scaled_values <- floor(random_values / sum(random_values) * k)
    
    # Adjust the result to make sure the sum is exactly k
    diff <- k - sum(scaled_values)
    
    # Randomly distribute the remaining difference (positive or negative)
    for (i in seq_len(abs(diff))) {
      index <- sample(1:n, 1)
      scaled_values[index] <- scaled_values[index] + sign(diff)
    }
    
  }
  
  return(scaled_values)
}


# Function to get incubation period density distribution given mean & sd
GenIncubDist <- function(params) {
  
  # params: numeric vector of length 2 in order c(log_mean, log_sd)
  
  # Set the mean & standard deviation given input parameters
  meanlog <- as.numeric(params[1])
  sdlog   <- as.numeric(params[2])
  
  # Derive log-normal density distribution for incubation period of 1-30 days
  dd <- dlnorm(x = 1:30, meanlog = meanlog, sdlog = sdlog)
  
  # Normalise so sums to 1
  norm.dd <- dd / sum(dd)
  
  # Return the normalised density distribution
  return(norm.dd)
  
}

# Function to estimate time-specific no. infections given incubation period & dated case numbers
# NB: Wrapper around `GenIncubDist()`
CasesToInfections <- function(chain_incub, dates, cases, seed) {
  
  # chain_incub: MCMC samples of incubation period with columns in order (mean, sd)
  # dates: Dates of recorded confirmed case numbers
  # cases: Confirmed case numbers for each of the dates (length(cases) == length(dates))
  # seed: numeric seed value to ensure repeatability using `sample()`
  
  # Set seed so sampling is repeatable
  set.seed(seed)
  
  # Get number of days to sample incubation period over
  num_days <- length(dates)
  
  # Error handling for mismatch in case & date vector length
  if (!(num_days == length(cases))) stop("dates & cases must be equal length vectors")
  
  # Sample a row (mean & sd values) from incubation period MCMC chain
  sampled_row <- as.numeric(chain_incub[sample(1:num_days, 1), ])
  
  # Derive the incubation period density distribution using the mean & sd
  inc_dist <- GenIncubDist(sampled_row)
  
  # Transform daily cases into infection incidence
  infections <- rep(NA, num_days)
  
  for (d in 1:num_days) {
    
    # Create 1-30 day range in which infection may have occurred (incubation period)
    date_range <- seq(d + 1, min(d + 30, num_days))
    
    # Get cases confirmed over these days
    cases_range <- cases[date_range]
    
    # Calculate current no. infections by integrating (future) no. cases over incubation period dens dist
    infections[d] <- sum(cases_range * inc_dist[1:length(date_range)])
  }
  
  # Create dataframe of sampled incubation period parameters & derived no. infections
  results.df <- data.frame(infections, meanlog = sampled_row[1], sdlog = sampled_row[2])
  
  # Return data frame
  return(results.df)
}

# Define wrapper function to transform cases to infections
getInfections <- function(s, i, case_data) {
  
  seed.value <- s * 1e6 + i * 1e3  # Unique seed for each iteration (unique for i & j values ≤1000)
  
  # Get case data for selected state
  selected_state <- states[s]
  state_cases <- case_data %>%
    filter(state %in% selected_state)
  
  state_cases$con_cases[is.na(state_cases$con_cases)] <- 0 # set NA to 0
  
  # Generate random poisson-distributed counts from case data (to reflect uncertainty in true cases)
  state_cases <- state_cases %>%
    mutate(pois_cases = GetPoisCases(cases = con_cases, seed = seed.value))
  
  # Transform cases to infections using sampled values from incubation period
  trans_cases.df <- CasesToInfections(chain_incub = chain_incub,
                                      dates = state_cases$date_mid,
                                      cases = state_cases$pois_cases,
                                      seed  = seed.value)
  
  # Format & return output
  out <- data.frame(date = state_cases$date_mid, inf = trans_cases.df$infections, i = i, state = selected_state)
  return(out)
}


# 4. Summarise posterior predictions for each state ----

# Read in combined simulations
sims <- vector("list", length(states))

for(i in seq_along(sims)) {
  sims[[i]] <- readRDS(paste0("Outputs/", states[i], "_sims_N=1000.rds"))
}
sims <- bind_rows(sims)

## i. Calculate summary statistics for each state across all iterations ----

# Create summarised model output
out_summary <- sims %>%
  group_by(state, iter) %>% 
  pivot_longer(
    cols = starts_with("rel"),
    names_to = "var",
    values_to = "value"
  ) %>%
  group_by(state, date, var) %>%
  reframe(
    state,
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>% 
  distinct()

# Transform cases to infections
infections <- vector("list", length(states))
for (s in seq_along(states)) {
  
  # Initialise parallel backend before foreach() loop
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores)
  
  # Simulate model for state for all posterior samples
  infections[[s]] <- foreach(i = 1:n_distinct(sims$iter), .packages = c("dplyr")) %dopar% {
    getInfections(i = i, s = s, case_data = formatted.case.data)
  }
  
  # Stop and clean up the parallel cluster
  doParallel::stopImplicitCluster()
}

# Summarise transformed data
summary_infections <- bind_rows(infections) %>%
  filter(!is.na(inf)) %>% 
  group_by(state, i,) %>% 
  mutate(rel_inf = inf / max(inf, na.rm = TRUE)) %>% 
  group_by(state, date) %>% 
  reframe(
    state,
    date,
    median_inf = median(rel_inf),
    lower_inf = quantile(rel_inf, 0.025),
    upper_inf = quantile(rel_inf, 0.975)
  ) %>%
  distinct()


## ii. Calculate summary statistics for each state and each iteration ----

# Interpolate observations to same dates as predictions
obs <- bind_rows(infections)
obs_int <- obs %>% 
  rename(iter = i) %>% 
  group_by(state, iter) %>% 
  complete(date=full_seq(date, 1)) %>% # fill missing dates in 1 day intervals
  mutate(inf_int = zoo::na.approx(inf, na.rm = FALSE),  # interpolate infections
         rel_inf_dat = inf_int / max(inf_int, na.rm = TRUE)) %>% # calculate relative infections 
  reframe(date, state, iter, rel_inf_dat) %>% 
  filter(date %in% unique(sims$date))

# Calculate peak dates in data for each iteration
data_peaks <- obs_int %>% 
  mutate(season_year = if_else(month(date) >= 7, year(date), year(date) - 1)) %>%
  group_by(state, iter, season_year) %>%
  slice_max(order_by = rel_inf_dat, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  rename(data_peak_date = date, data_peak_val = rel_inf_dat)

# Calculate peak dates in predictions for each iteration
model_peaks <- sims %>% 
  reframe(date, state, iter, rel_inf_mod = rel_inf) %>% 
  mutate(season_year = if_else(month(date) >= 7, year(date), year(date) - 1)) %>%
  group_by(state, iter, season_year) %>%
  slice_max(order_by = rel_inf_mod, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  rename(model_peak_date = date, model_peak_val = rel_inf_mod)

# Combine both datasets
combined_peaks <- left_join(model_peaks, data_peaks, by = c("state", "iter", "season_year"))

# Calculate differences in peak dates & peak values
peak_diffs <- combined_peaks %>%
  filter(season_year != 2017, season_year != 2020) %>%  # remove 2017 peak & COVID year
  group_by(state, iter, season_year) %>% 
  mutate(peak_diff_days   = as.numeric(model_peak_date - data_peak_date),
         peak_diff_weeks  = peak_diff_days / 7,
         peak_diff_months = peak_diff_days / 30.44,
         peak_diff_val    = model_peak_val - data_peak_val) %>% 
  ungroup() %>% 
  mutate(year_label = paste0(substr(season_year, 3, 4), "/", substr(season_year + 1, 3, 4))) %>% 
  filter(year_label != "25/26")

# Difference in peak date: summarise median + 95% credible intervals
peak_summary_dates <- peak_diffs %>% 
  group_by(state, year_label) %>%
  reframe(med_diff = median(peak_diff_months),
          lower_diff = quantile(peak_diff_months, 0.025),
          upper_diff = quantile(peak_diff_months, 0.975))

peak_summary_dates$state <- factor(peak_summary_dates$state, levels = states)

# Difference in peak magnitude: summarise median + 95% credible intervals
peak_summary_vals <- peak_diffs %>% 
  group_by(state, year_label) %>%
  reframe(med_diff = median(peak_diff_val),
          lower_diff = quantile(peak_diff_val, 0.025),
          upper_diff = quantile(peak_diff_val, 0.975))

peak_summary_vals$state <- factor(peak_summary_vals$state, levels = states)


# 5. Make plots ----

## (i). Plot predicted vs. observed time series ----

# Generate palette making state names
palette <- RColorBrewer::brewer.pal(n = length(states), name = "Set2")
names(palette) <- states
palette <- c(palette, Overall = "black") # add black to the palette for overall summary statistic

# Combine datasets with an identifier for data vs. model
model_long <- out_summary %>%
  mutate(source = "Model") %>%
  rename(median_val = median, lower_val = lower, upper_val = upper)

data_long <- summary_infections %>%
  mutate(source = "Observed") %>%
  rename(median_val = median_inf, lower_val = lower_inf, upper_val = upper_inf)

plot_dat <- bind_rows(model_long, data_long)
plot_dat$state <- factor(plot_dat$state, levels = states)

# Make facet plot
alpha_val <- 0.5

# Create shading for 2020 (not included in analyses)
exclude_df <- data.frame(xmin = ymd("2021-01-01") - months(6), 
                       xmax = ymd("2021-01-01") + months(6), 
                       ymin = -Inf, ymax = Inf )

tseries_plot <- ggplot() +
  # Shade excluded year
  geom_rect(data = exclude_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "darkgrey", alpha = 0.8) +
  # Observed data (grey)
  geom_ribbon(data = filter(plot_dat, source == "Observed"),
              aes(x = date, ymin = lower_val, ymax = upper_val, group = state),
              fill = "lightgrey", alpha = alpha_val) +
  geom_line(data = filter(plot_dat, source == "Observed"),
            aes(x = date, y = median_val, group = state),
            colour = "grey40") +
  # Model predictions (state colours)
  geom_ribbon(data = filter(plot_dat, source == "Model"),
              aes(x = date, ymin = lower_val, ymax = upper_val, fill = state),
              alpha = alpha_val-.2, colour = NA) +
  geom_line(data = filter(plot_dat, source == "Model"),
            aes(x = date, y = median_val, colour = state)) +
  labs(x = "Year", y = "Relative size", col = "State", fill = "State") +
  scale_x_date(breaks = seq(
    from = as.Date(paste0(lubridate::year(min(plot_dat$date)), "-01-01")),
    to   = max(plot_dat$date),
    by   = "1 year"), labels = scales::date_format("%Y")) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 3)) +
  facet_wrap(~state, nrow = 5) +
  theme_Publication() +
  scale_color_manual(values = palette, guide = "legend") +
  scale_fill_manual(values = palette, guide = "legend")


## (ii). Plot Nigerian states on map ----

# Load Nigeria state shape files
states_sf <- sf::st_as_sf(  # Convert shapefile to sf object
  rgeoboundaries::geoboundaries(  # NB requires internet connection
    country = "Nigeria",
    adm_lvl = "adm1",  # Get state borders
    type = NULL,
    quiet = TRUE
  )
)

state_map <- states_sf %>% 
  ggplot(aes(fill = shapeName)) +
  geom_sf() +
  scale_fill_manual(values = palette, name = "State") + 
  theme_Publication() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))


## (iii). Plot predicted vs. observed peak timing ----

# Create shading for years on graph
odd_labels <- unique(peak_diffs$year_label)[c(FALSE, TRUE)]
shade_df <- data.frame(year_label = odd_labels)

# Summarise overall metric of difference in peak date
total_peak_summary_dates <- peak_diffs %>% 
  reframe(year_label = "Overall",
          med_diff = median(peak_diff_months),
          lower_diff = quantile(peak_diff_months, 0.025),
          upper_diff = quantile(peak_diff_months, 0.975),
          state = "Overall")
# range(peak_summary_dates$med_diff)
# which(peak_summary_dates$lower_diff <= 0) %>% length() / nrow(peak_summary_dates)

# Plot difference in peak dates (median + 95% CrI, dot + whisker)
peak_date_plot <- full_join(peak_summary_dates, total_peak_summary_dates) %>%
  ggplot(aes(year_label, med_diff, col = state, fill = state)) +
  geom_tile(data = shade_df,
            aes(x = year_label, y = 0),
            inherit.aes = FALSE,
            fill = "gray90", alpha = 0.5,
            width = 1, height = Inf) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_errorbar(aes(ymin = lower_diff, ymax = upper_diff), width = .2, position = position_dodge(width = dodge_width)) + 
  geom_point(position = position_dodge(width = dodge_width)) + 
  labs(x = "Year", y = "Months between peaks (predicted - observed)", col = "State", fill = "State") + 
  theme_Publication() + 
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-3, 3, 1.5), limits = c(-3, 3), expand = c(0, 0)) + 
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette)


## (iv). Plot predicted vs. observed peak magnitude ----

# Summarise overall metric of difference in peak magnitude
total_peak_summary_val <- peak_diffs %>% 
  reframe(year_label = "Overall",
          med_diff = median(peak_diff_val),
          lower_diff = quantile(peak_diff_val, 0.025),
          upper_diff = quantile(peak_diff_val, 0.975),
          state = "Overall")

# Plot difference in peak dates (median + 95% CrI, dot + whisker)
peak_height_plot <- full_join(peak_summary_vals, total_peak_summary_val) %>% 
  ggplot(aes(year_label, med_diff, col = state, fill = state)) +
  geom_tile(data = shade_df,
            aes(x = year_label, y = 0),
            inherit.aes = FALSE,
            fill = "gray90", alpha = 0.5,
            width = 1, height = Inf) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_errorbar(aes(ymin = lower_diff, ymax = upper_diff), width = .2, position = position_dodge(width = dodge_width)) + 
  geom_point(position = position_dodge(width = dodge_width)) + 
  labs(x = "Year", y = "Peak size difference (predicted - observed)", col = "State", fill = "State") +
  theme_Publication() +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-1, 1, length.out = 5), limits = c(-1, 1), expand = c(0, 0)) +
  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette)


## (v). Make panel plot ----
design <- "AABB
           CCDD"

wrap_plots(A = state_map, B = tseries_plot, C = peak_date_plot, D = peak_height_plot, 
           design = design) &
  theme(legend.position = "none") & 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18), text = element_text(size = 11))

ggsave("Plots/posterior_predictions/nigeria_peaks.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("Plots/posterior_predictions/nigeria_peaks.png", bg = "white", width = 8, height = 8, dpi = 1200)
