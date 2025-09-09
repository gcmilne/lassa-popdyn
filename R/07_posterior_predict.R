#---------------------------------------------------------#
# Posterior predictions using fitted infection parameters #
#---------------------------------------------------------#

# 0. Script purpose ----
# - Loads climate-demographic & fine grid estimates of transmission parameter posteriors
# - Performs posterior predictions for dynamics of seroprevalence & demographics
#   using Bayesian melding
# - Makes & saves plots of outputs
# - Calculates predictive performance of seroprevalence model on hold-out data


# 1. Set up IPM for simulation (load packages, data, set parameters) ----
use_cluster <- FALSE  # Use the high performance cluster?
if (use_cluster) {  # Set directory for HPC
  proj_dir <- "/home/gregm/mastomys_modelling"
  setwd(proj_dir)
}
script_path <- ifelse(use_cluster, "scripts", "R")
source(paste0(script_path, "/setup_IPM.R"))

# Account for uncertainty in IFA sensitivity & specificity during fitting?
use_se_sp <- FALSE

# Load additional packages
pacman::p_load(purrr, foreach, patchwork, tidybayes, bayesplot, brms)

# Source custom ggplot theme
source("R/ggplot_theme.R")


# 2. Load data ----

# Load seroprevalence data (inc. validation year left out of fitting)
seroprev.dat <- readRDS(paste0(data_path, "/MORV_seroprev_formatted.rds"))
seroprev.dat.sub <- seroprev.dat %>% filter(FER == FALSE)
seroprev.dat.ad  <- seroprev.dat %>% filter(FER == TRUE)

# Re-generate climate lag grid ensuring model runs past final seroprevalence date
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

# Load & combine fine grid transmission posteriors
file_list <- list.files(path = "Outputs/model_fits/fine_grid/",
                        full.names = TRUE, 
                        ignore.case = TRUE)

# Join together, calculate weights for importance sampling
inf_post <- file_list %>%
  map(readRDS) %>%
  bind_rows() %>% 
  filter(joint_posterior != -Inf & !is.na(weight)) %>%
  group_by(post_iter) %>%
  mutate(weight = exp(joint_posterior - matrixStats::logSumExp(joint_posterior))) %>%
  arrange(desc(weight)) %>%
  ungroup()


# 3. Posterior predictions using Bayesian melding ----

## Refs: https://www.tandfonline.com/doi/abs/10.1080/01621459.2000.10474324 https://sti.bmj.com/content/84/Suppl_1/i11 

if (file.exists("Outputs/posterior_predict_samples_N=1000.rds")) {
  posterior_samples <- readRDS("Outputs/posterior_predict_samples_N=1000.rds")
} else {
  
  ## i. Compute marginal likelihood per φ sample (post_iter) ----
  # Sum over the ψ grid for that φ
  marginals <- inf_post %>%
    group_by(post_iter) %>%
    summarise(
      log_weights = log(sum(exp(joint_posterior - max(joint_posterior)))) + max(joint_posterior),
      .groups = "drop"
    ) %>%
    mutate(
      weight = exp(log_weights - max(log_weights))
    )
  
  # Normalize weights
  marginals <- marginals %>%
    mutate(
      weight_norm = weight / sum(weight)
    )
  
  # Check
  sum(marginals$weight_norm)  # should be 1
  
  ## ii. Resample φ indices ----
  set.seed(123)
  N_samples <- 1000
  phi_samples <- sample(
    marginals$post_iter,
    size = N_samples,
    replace = TRUE,
    prob = marginals$weight_norm
  )
  
  ## iii. For each φ sample, sample ψ given φ ----
  # Need conditional weights for that φ's grid
  sampled_phi <- vector("list", N_samples)
  
  for (n in seq_len(N_samples)) {
    
    phi_id <- phi_samples[n]
    
    # Grid rows for this φ sample
    grid_phi <- inf_post %>% filter(post_iter == phi_id)
    
    # Conditional weights over ψ grid
    log_post <- grid_phi$joint_posterior
    log_post <- log_post - max(log_post)  # For stability
    weights <- exp(log_post)
    
    if (sum(weights) == 0) {
      warning(sprintf("All weights zero for post_iter=%d", phi_id))
      next
    }
    
    # Normalise weights to sum to 1
    weights <- weights / sum(weights)
    
    # Sample one ψ grid point
    sampled_idx <- sample(1:nrow(grid_phi), size = 1, prob = weights)
    
    sampled_phi[[n]] <- grid_phi[sampled_idx, ]
  }
  
  ## iv. Combine samples into a single dataframe ----
  posterior_samples <- bind_rows(sampled_phi)
  # saveRDS(posterior_samples, "Outputs/posterior_predict_samples_N=1000.rds")
}

## v. Simulate model with sampled parameters ----

# Define function for posterior prediction
posterior_predict <- function(j, sampled_inf_post, demog_posts, inf_pars) {
  
  # Get infection posterior for input index
  inf_post <- sampled_inf_post[j, ]
  
  # Get matching demographic posterior samples
  sampled_demog_posts <- demog_posts %>%
    filter(
        (process == "surv"     & draw == inf_post$surv.draw.id) |
        (process == "rec"      & draw == inf_post$rec.draw.id)  |
        (process == "grow"     & draw == inf_post$grow.draw.id) |
        (process == "male_pop" & draw == inf_post$malepop.draw.id)
    )
  
  # Run IPM with climate-demographic & transmission inputs & store
  inf_pars$trans_rate     <- inf_post$trans_rate
  inf_pars$vertical_trans <- inf_post$vertical_trans
  
  mod <- RunRodentModel(model = sampled_demog_posts, env = lag_clim, cmr.data = cmr.data,
                        demog.pars = demog.pars, inf.pars = inf_pars,
                        n_init = n_init, start_date = start_date, end_date = end_date,
                        seed = 1, summarised = TRUE) %>%
    # Store time-specific modelled population size & seroprevalence (subadults & adults)
    group_by(date) %>%
    mutate(seroprev = (D_all - D_adult_s - Dsub_s1 - Dsub_s2 - J_s) / D_all) %>% 
    select(date, mean_adult_wgt, mod.pop, seroprev, seroprev.sub, seroprev.adult, contains(c("J_", "D_adult_", "Dsub_"))) %>% 
    distinct() %>%
    filter(date >= min(seroprev.dat$date_mid) & date <= max(seroprev.dat$date_mid)) %>% 
    mutate(trans_rate = inf_pars$trans_rate,
           vertical_trans = inf_pars$vertical_trans,
           iter = j,
           )
  
  return(mod)
}

# Posterior predictions
if (file.exists("Outputs/posterior_predictions_N=1000.rds")) {
  preds <- readRDS("Outputs/posterior_predictions_N=1000.rds")
} else {  
  
  # Initialise parallel backend before foreach() loop
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores)
  
  # Run posterior predictions for all samples
  preds <- foreach(i = 1:1000, .packages = c("dplyr", "truncnorm", "lubridate")) %dopar% {
    posterior_predict(j = i, sampled_inf_post = posterior_samples, demog_posts = posts, 
                      inf_pars = inf.pars)
  }
  
  # Stop and clean up the parallel cluster
  doParallel::stopImplicitCluster()
  
  # Save posterior predictions
  preds <- bind_rows(preds)
  saveRDS(preds, "Outputs/posterior_predictions_N=1000.rds")
}

# Load seroprevalence data (inc. validation year left out of fitting)
seroprev.dat <- readRDS(paste0(data_path, "/MORV_seroprev_formatted.rds"))
seroprev.dat.sub <- seroprev.dat %>% filter(FER == FALSE)
seroprev.dat.ad  <- seroprev.dat %>% filter(FER == TRUE)

# 4. Calculate posterior predictive medians & 95% credible intervals ----
out_summary <- bind_rows(preds) %>%
  reframe(date, seroprev.sub, seroprev.adult, iter) %>%
  distinct() %>%
  pivot_longer(
    cols = starts_with("seroprev"),
    names_to = "population",
    values_to = "seroprevalence"
  ) %>%
  group_by(date, population) %>%
  summarise(
    median = median(seroprevalence),
    lower = quantile(seroprevalence, 0.025),
    upper = quantile(seroprevalence, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    population = case_when(
      population == "seroprev.adult" ~ "Adults",
      population == "seroprev.sub"  ~ "Sub-adults"
    ),
    is.fitted = if_else(year(date) < 2017, "Fit", "Validation")
  )

# Interpolate model output to be daily
full_date_seq <- tibble(date = seq(from = min(out_summary$date), to = max(out_summary$date), by = "1 day"))

out_summary_interp <- out_summary %>%
  group_by(population) %>%
  group_modify(
    ~{
      full_dates <- full_date_seq$date
      tibble(
        date = full_dates,
        population = .x$population[1],
        median = approx(x = .x$date, y = .x$median, xout = full_dates, rule = 2)$y,
        lower  = approx(x = .x$date, y = .x$lower,  xout = full_dates, rule = 2)$y,
        upper  = approx(x = .x$date, y = .x$upper,  xout = full_dates, rule = 2)$y
      )
    }
  ) %>%
  ungroup() %>% 
  # Visualise which part of model was fitted vs. used for prediction
  mutate(is.fitted = if_else(year(date) < 2017, "Fit", "Predicted"))

# Join data together for plot
seroprev.dat.long <- bind_rows(seroprev.dat.sub %>% mutate(population = "Sub-adults"),
                               seroprev.dat.ad  %>% mutate(population = "Adults"))

# Plot modelled seroprevalence against data
seroprev_plot <- out_summary_interp %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100, fill = is.fitted), alpha = 0.3) +
  geom_line(aes(y = median*100, col = is.fitted), linewidth = 0.4) +
  geom_point(data = seroprev.dat.long, 
             aes(x = date_mid, y = mean.seroprev*100), 
             col = "darkgrey", size = 0.8, inherit.aes = FALSE) +
  geom_errorbar(data = seroprev.dat.long, 
                aes(x = date_mid, ymin = lower.seroprev*100, ymax = upper.seroprev*100),
                col = "darkgrey",
                width = 0,
                linewidth = 0.2,
                inherit.aes = FALSE) +
  facet_wrap(~factor(population, levels = c("Sub-adults", "Adults")), nrow = 1) +
  labs(x = "Date", y = "Seroprevalence (%)", fill = "", col = "") +
  scale_x_date(breaks = seq(from = as.Date(paste0(lubridate::year(min(seroprev.dat$date_mid)), "-01-01")),
                 to = max(seroprev.dat$date_mid)+days(32), by = "1 year"), labels = scales::date_format("%Y")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "RdBu")[c(3, 1)]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "RdBu")[c(3, 1)])

#### Plot no. pregnancies over time
preg_dat <- bind_rows(preds) %>%
  reframe(date, J_all, iter) %>%
  distinct() %>%
  pivot_longer(
    cols = starts_with("J_"),
    names_to = "pop",
    values_to = "num_preg"
  ) %>%
  group_by(date) %>%
  summarise(
    median = median(num_preg),
    lower = quantile(num_preg, 0.025),
    upper = quantile(num_preg, 0.975),
    .groups = "drop"
  )

preg_dat %>% 
  ggplot(aes(x = date, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  geom_line() +
  scale_x_date(breaks = seq(from = as.Date(paste0(lubridate::year(min(seroprev.dat$date_mid)), "-01-01")),
                            to = max(seroprev.dat$date_mid)+days(32), by = "1 year"), labels = scales::date_format("%Y")) +
  labs(x = "Date", y = "Modelled no. pregnancies") +
  theme_Publication()

ggsave("Plots/posterior_predictions/num_preg.pdf", height = 6, width = 6, units = "in", device = cairo_pdf)
ggsave("Plots/posterior_predictions/num_preg.png", height = 6, width = 6, units = "in", dpi = 600, bg = "white")

# Summarise relative change in peak no. pregnancies by year
preg_dat %>% 
  filter(date >= min(out_summary_interp$date) & date <= max(out_summary_interp$date)) %>% 
  group_by(year(date)) %>% 
  filter(median %in% max(median)) %>%
  ungroup() %>% 
  mutate(prev_med = dplyr::lag(median),
         prev_lower = dplyr::lag(lower), prev_upper = dplyr::lag(upper)) %>% 
  reframe(year(date), 
          perc_change = (median / prev_med)*100,
          perc_change_lower = (lower / prev_lower)*100,
          perc_change_upper = (upper / prev_upper)*100)

#### Plot infection compartments by demographic group over time
compartment_data <- bind_rows(preds) %>%
  group_by(date, iter) %>%
  # Combine sub-adults within same compartment
  mutate(
    Dsub_s = Dsub_s1 + Dsub_s2,
    Dsub_h = Dsub_h1 + Dsub_h2,
    Dsub_v = Dsub_v1 + Dsub_v2,
    Dsub_r = Dsub_r1 + Dsub_r2
  ) %>%
  select(!c(Dsub_s1, Dsub_s2, Dsub_h1, Dsub_h2, Dsub_v1, Dsub_v2, Dsub_r1, Dsub_r2,
            contains("_all"))) %>%
  pivot_longer(
    cols = contains(c("J_", "D_adult_", "Dsub_")),
    names_to = "population",
    values_to = "value"
  ) %>%
  mutate(
    compartment = case_when(
      str_detect(population, "_s") ~ "S",
      str_detect(population, "_h") ~ "H",
      str_detect(population, "_v") ~ "V",
      str_detect(population, "_r") ~ "R",
      str_detect(population, "_m") ~ "M"
    ),
    group = case_when(
      # str_detect(population, "J_") ~ "Pregnant",
      # str_detect(population, "D_adult_") ~ "Non-pregnant adults",
      str_detect(population, "J_") | str_detect(population, "D_adult_") ~ "Adults",
      str_detect(population, "Dsub_") ~ "Sub-adults"
    )
  ) %>%
  ungroup()

# Make compartment an ordered factor
compartment_data <- compartment_data %>%
  mutate(compartment = factor(compartment, levels = c("S", "H", "V", "R", "M")))

# Summarise data by median & 95% crIs
summarised_compartment_data <- compartment_data %>%
  group_by(date, group, compartment) %>%
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )

#### Plot population size by compartment
inf_plot <- summarised_compartment_data %>%
  filter(compartment %in% c("H", "V")) %>%   # just infected
  ggplot(aes(x = date, fill = compartment, color = compartment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.4) +
  facet_wrap(~factor(group, levels = c("Sub-adults", "Adults")), scales = "free", ncol = 2) +
  labs(x = "Date", y = "Population size", fill = "", color = "") +
  scale_x_date(breaks = seq(from = as.Date(paste0(lubridate::year(min(seroprev.dat$date_mid)), "-01-01")),
                            to = max(seroprev.dat$date_mid)+days(32), by = "1 year"), labels = scales::date_format("%Y")) +
  scale_fill_viridis_d(option = "H") + scale_color_viridis_d(option = "H")

#### Summarise change in no. infected by year
summarised_compartment_data %>%
  filter(compartment %in% c("H"), year(date) %in% 2015:2016) %>%
  group_by(compartment, group, year(date)) %>%
  filter(median %in% max(median))

#### Plot transmission parameter posteriors
posteriors <- bind_rows(preds) %>%
  ungroup() %>% 
  reframe(iter, phi = vertical_trans, beta = trans_rate) %>% 
  distinct() %>%
  dplyr::select(!iter) %>% 
  mutate(group = "Posterior")

set.seed(001)
priors <- data.frame(
  name  = c(rep("phi", 10000), rep("beta", 10000)),
  value = c(rbeta(10000, shape1 = 20, shape2 = 2), seq(0.0001, 0.0015, length.out = 10000)),
  group = "Prior"
)

# Join prior to posterior (in long format)
prior_post <- posteriors %>% 
  pivot_longer(cols = c(phi, beta)) %>% 
  bind_rows(priors) %>% 
  mutate(group = factor(group, levels = c("Prior", "Posterior")))

# Define custom parameter labels
par_labs <- c("phi" = "\u03C6", "beta" = "\u03B2")

# Summarise posteriors
posterior_summary <- prior_post %>%
  filter(group == "Posterior") %>%
  group_by(name) %>%
  summarise(
    median = median(value),
    lower80 = quantile(value, 0.10),
    upper80 = quantile(value, 0.90),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975),
    .groups = "drop"
  )

#### Plot transmission parameter posterior distributions
trans_post_plot <- prior_post %>%
  ggplot(aes(x = value, fill = group)) +
  geom_density(alpha = 0.8, col = NA) +
  geom_segment(  # Outer interval (thin)
    data = posterior_summary,
    aes(x = lower95, xend = upper95, y = 0, yend = 0), linewidth = 0.4, inherit.aes = FALSE
  ) +
  geom_segment(  # Inner interval (thicker)
    data = posterior_summary, 
    aes(x = lower80, xend = upper80, y = 0, yend = 0), linewidth = 0.8, inherit.aes = FALSE
  ) +
  geom_point(
    data = posterior_summary, 
    aes(x = median, y = 0), size = 1.5, inherit.aes = FALSE
  ) +
    facet_wrap(~name, scales = "free", labeller = labeller(name = as_labeller(par_labs, label_parsed))) +
  labs(x = "Parameter value", y = "Density", fill = "") + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(2, "PuOr")[c(1, 3)])
  

#### Panel plot: seroprevalence fits, transmission posteriors & posterior predictions
patchwork::wrap_plots(trans_post_plot, seroprev_plot, inf_plot, nrow = 3) +
  plot_annotation(tag_levels = 'A',
                  theme = theme(plot.tag = element_text(size = 12, face = "bold"))) & 
  theme_Publication()

ggsave("Plots/fits_postpreds_plot.pdf", height = 8, width = 8, units = "in", device = cairo_pdf)
ggsave("Plots/fits_postpreds_plot.png", height = 8, width = 8, units = "in", dpi = 600, bg = "white")


#### Plot proportion of population in each compartment
compartment_data %>%
  group_by(date, iter, group) %>%
  mutate(
    prop_compartment = value / sum(value)
  ) %>%
  group_by(date, group, compartment) %>%
  summarise(
    median = median(prop_compartment),
    lower = quantile(prop_compartment, 0.025),
    upper = quantile(prop_compartment, 0.975),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = date, fill = compartment, color = compartment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  facet_wrap(~group, ncol = 1) +
  scale_x_date(breaks = seq(from = as.Date(paste0(lubridate::year(min(seroprev.dat$date_mid)), "-01-01")),
                            to = max(seroprev.dat$date_mid)+days(32), by = "1 year"), labels = scales::date_format("%Y")) +
  scale_fill_viridis_d() + scale_color_viridis_d() +
  labs(x = "Date", y = "Proportion of population", fill = "Compartment", color = "Compartment") +
  theme_Publication()

ggsave("Plots/posterior_predictions/prop_compartments.pdf",
       height = 6, width = 6, units = "in", device = cairo_pdf)
ggsave("Plots/posterior_predictions/prop_compartments.png",
       height = 6, width = 6, units = "in", dpi = 1200, bg = "white")

#### Plot time series of prevalence vs seroprevalence
prev_sero_dat <- bind_rows(preds) %>%
  group_by(date, iter) %>%
  mutate(prev = (D_adult_h + D_adult_v +  J_h + J_v + Dsub_h1 + Dsub_h2 + Dsub_v1 + Dsub_v2) / (D_adult_all + J_all + Dsub_all)) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(
    median_prev = median(prev),
    lower_prev = quantile(prev, 0.025),
    upper_prev = quantile(prev, 0.975),
    median_sero = median(seroprev),
    lower_sero = quantile(seroprev, 0.025),
    upper_sero = quantile(seroprev, 0.975)
  ) %>%
  pivot_longer(
    cols = c(median_prev, median_sero),
    names_to = "variable",
    values_to = "median"
  ) %>%
  mutate(
    lower = ifelse(variable == "median_prev", lower_prev, lower_sero),
    upper = ifelse(variable == "median_prev", upper_prev, upper_sero),
    variable = recode(variable, median_prev = "Prevalence", median_sero = "Seroprevalence")
  )

prev_sero_dat %>%
  ggplot(aes(x = date, y = median, color = variable, fill = variable)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Date", y = "Proportion", color = "", fill = "") +
  scale_x_date(breaks = seq(from = as.Date(paste0(lubridate::year(min(seroprev.dat$date_mid)), "-01-01")),
                            to = max(seroprev.dat$date_mid)+days(32), by = "1 year"), labels = scales::date_format("%Y")) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
  scale_color_viridis_d() + scale_fill_viridis_d() + 
  theme_Publication()

ggsave("Plots/posterior_predictions/prev_vs_seroprev_tseries.pdf",
       height = 6, width = 6, units = "in", device = cairo_pdf)
ggsave("Plots/posterior_predictions/prev_vs_seroprev_tseries.png",
       height = 6, width = 6, units = "in", dpi = 1200, bg = "white")

# Dates of peak prevalence vs. peak seroprevalence
prev_sero_dat %>% 
  group_by(year(date), variable) %>% 
  filter(median == max(median)) %>% View()

# ymd("2015-03-25") - ymd("2014-09-10")  # median 196 days between peak seroprev & peak prev (2014/15 season)
# ymd("2017-03-22") - ymd("2016-09-07")  # median 196 days between peak seroprev & peak prev (2016/17 season)


# 5. Calculate agreement between predicted and observed seroprevalence ----

# Get modelled output for 2017 (2017 data not used to fit model)
out <- bind_rows(preds) %>%
  reframe(date, seroprev.sub, seroprev.adult, iter) %>%
  distinct() %>% 
  filter(date >= ymd("2016-12-28"))

# Interpolate model output to be daily
full_date_seq <- tibble(date = seq(from = min(out$date), to = max(out$date), by = "1 day"))

out_interp <- out %>%
  group_by(iter) %>%
  group_modify(
    ~{
      full_dates <- full_date_seq$date
      tibble(
        date = full_dates,
        mod.seroprev.sub = approx(x = .x$date, y = .x$seroprev.sub, xout = full_dates, rule = 2)$y,
        mod.seroprev.adult = approx(x = .x$date, y = .x$seroprev.adult, xout = full_dates, rule = 2)$y
      )
    }
  ) %>%
  ungroup()

# Calculate no. seropositive implied by modelled seroprevalence
subs <- out_interp %>%   # Sub-adults
  left_join(seroprev.dat.sub, by = c("date" = "date_mid")) %>% 
  filter(!is.na(mean.seroprev)) %>% 
  reframe(date, mod_seroprev = mod.seroprev.sub, obs_seropos = seropos, obs_tested = tested, group = "Sub-adults")
ads <- out_interp %>%   # Adults
  left_join(seroprev.dat.ad, by = c("date" = "date_mid")) %>% 
  filter(!is.na(mean.seroprev)) %>% 
  reframe(date, mod_seroprev = mod.seroprev.adult, obs_seropos = seropos, obs_tested = tested, group = "Adults")

both <- bind_rows(subs, ads)

# Get median model predictions
df_med <- both %>%
  group_by(date,group) %>%
  summarise(mod_seroprev=median(mod_seroprev),obs_tested=sum(obs_tested)/1000,obs_seropos=sum(obs_seropos)/1000,.groups="keep") %>%
  ungroup()

# Predict observed seroprevalence with median modelled seroprevalence
mod <- brm(data=df_med,
           formula=obs_seropos|trials(obs_tested)~mod_seroprev,
           family=binomial("identity"),
           chains=4,
           iter=2000,
           warmup=1000,
           thin=1)

# Null model (only fit intercept)
mod_null <- brm(data=df_med,
                formula=obs_seropos|trials(obs_tested)~1,
                family=binomial("identity"),
                chains=4,
                iter=2000,
                warmup=1000,
                thin=1)

# Compare model to null model
loo_mod  <- loo(mod)
loo_null <- loo(mod_null)
loo_compare(loo_mod, loo_null)
