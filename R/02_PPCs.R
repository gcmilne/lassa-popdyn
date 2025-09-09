#-----------------------------------------------------------------------------#
# Do posterior predictive checks (PPCs) & plot demographic process model fits #
#-----------------------------------------------------------------------------#

# 0. Session details ----
# - Read in fitted demographic process models
# - Perform posterior predictive checks & make plots


# 1. Load packages ----
pacman::p_load(ggplot2, lubridate, dplyr, stringr, tidyr, rstan, bayesplot, tidybayes, purrr,
               ggridges, patchwork, brms, loo, bridgesampling, boot, doParallel, posterior, forcats)

# Source custom ggplot theme
source("R/ggplot_theme.R")
color_scheme_set("purple") # set colour scheme


# 2. Load data ----

# Function to apply unit scaling to covariates (mean of 0, SD of 1)
unitScale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

## Demographic data (from demo_formatting_Morogoro.R)
survival <- read.csv("Data/survival_data.csv", stringsAsFactors = FALSE) %>% 
  mutate(DATE = ymd(DATE))

recruitment <- read.csv("Data/recruitment_data.csv", stringsAsFactors = FALSE) %>%
  mutate(DATE = ymd(DATE))

growth <- read.csv("Data/growth_data.csv", stringsAsFactors = FALSE) %>%  # apply same scaling as before model fitting
  mutate(DATE1 = ymd(DATE1), DATE2 = ymd(DATE2),
         mean_W1 = mean(W1),
         sd_W1 = sd(W1),
         mean_W2 = mean(W2),
         sd_W2 = sd(W2),
         across(c(W1, W2), unitScale))

# Minimum number alive (MNA) data (from demo_formatting_Morogoro.R)
MNA <- read.csv("Data/MNA.by.sex.csv") %>% 
  mutate(month = month(trap.date))

## Stan model fits
rec_fit  <- readRDS("Outputs/stan_fits/surv_stan_bestfit.rds")
surv_fit <- readRDS("Outputs/stan_fits/rec_stan_bestfit.rds")
grow_fit <- readRDS("Outputs/stan_fits/grow_stan_bestfit.rds")
mf_fit   <- readRDS("Outputs/stan_fits/mf_pop.rds")


# 3. Posterior predictive checks ----

# Aim: Compare distribution of observations to simulated data 
# Ref: https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html

## (a) Recruitment ----

# Generate simulated data (yrep) for PPCs using posterior values for eta
rec_post <- rec_fit %>% tidy_draws()  # get posterior of best-fitting model
eta_mat   <- rec_post %>%  # Convert posterior eta values to matrix
  select(starts_with("eta[")) %>%
  as.matrix()  # [draws x N]

# Simulate yrep ~ bernoulli_logit(eta)
set.seed(1234)
yrep_rec <- matrix(
  rbinom(length(eta_mat), size = 1, prob = plogis(eta_mat)),
  nrow = nrow(eta_mat)
)
rm(eta_mat)  # clear environment

# Do PPCs
y <- as.numeric(recruitment$REC)

ppc_tot_mean   <- ppc_stat(y = y, yrep = yrep_rec, stat = mean, binwidth = 0.005) &
  scale_x_continuous(n.breaks = 5) & 
  labs(x = "Mean proportion pregnant")

ppc_month_mean <- ppc_stat_grouped(y = y, yrep = yrep_rec, group = month(recruitment$DATE, label = TRUE), stat = mean) & 
  scale_x_continuous(limits = c(0, 0.6)) &
  labs(x = "Mean proportion pregnant")

ppc_tot_sd   <- ppc_stat(y = y, yrep = yrep_rec, stat = sd, binwidth = 0.005) &
  scale_x_continuous(n.breaks = 5) & 
  labs(x = "SD proportion pregnant")

ppc_month_sd <- ppc_stat_grouped(y = y, yrep = yrep_rec, group = month(recruitment$DATE, label = TRUE), stat = sd) & 
  scale_x_continuous(limits = c(0, 0.6)) &
  labs(x = "SD proportion pregnant")

wrap_plots(list(ppc_tot_mean, ppc_month_mean, ppc_tot_sd, ppc_month_sd), nrow = 2) & 
  plot_annotation(tag_levels = 'A') & 
  theme_Publication()

# Save figure
ggsave("Plots/model_fits/PPC_rec.pdf", width = 10, height = 10, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/PPC_rec.png", width = 10, height = 10, units = "in", dpi = 1200, bg = "white")


## (b) Survival ----

# Generate yrep for PPCs using posterior values for eta
surv_post <- surv_fit %>% tidy_draws()  # get posterior of best-fitting model
eta_mat   <- surv_post %>%  # Convert posterior eta values to matrix
  select(starts_with("eta[")) %>%
  as.matrix()  # [draws x N]

# Simulate yrep ~ bernoulli_logit(eta)
set.seed(1234)
yrep_surv <- matrix(
  rbinom(length(eta_mat), size = 1, prob = plogis(eta_mat)),
  nrow = nrow(eta_mat)
)
rm(eta_mat)  # clear environment

# Do PPCs
y <- as.numeric(survival$SURV)

ppc_tot_mean   <- ppc_stat(y = y, yrep = yrep_surv, stat = mean, binwidth = 0.005) &
  scale_x_continuous(n.breaks = 5) &
  labs(x = "Mean proportion recaptured")

ppc_month_mean <- ppc_stat_grouped(y = y, yrep = yrep_surv, group = month(survival$DATE, label = TRUE), stat = mean) & 
  scale_x_continuous(limits = c(0.25, 0.6)) & 
  labs(x = "Mean proportion recaptured")

ppc_tot_sd   <- ppc_stat(y = y, yrep = yrep_surv, stat = sd, binwidth = 0.005) &
  scale_x_continuous(n.breaks = 5) &
  labs(x = "SD proportion recaptured")

ppc_month_sd <- ppc_stat_grouped(y = y, yrep = yrep_surv, group = month(survival$DATE, label = TRUE), stat = sd) & 
  scale_x_continuous(limits = c(0.46, 0.51), n.breaks = 3) & 
  labs(x = "SD proportion recaptured")

wrap_plots(list(ppc_tot_mean, ppc_month_mean, ppc_tot_sd, ppc_month_sd), nrow = 2) & 
  plot_annotation(tag_levels = 'A') & 
  theme_Publication()

# Save figure
ggsave("Plots/model_fits/PPC_surv.pdf", width = 10, height = 10, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/PPC_surv.png", width = 10, height = 10, units = "in", dpi = 1200, bg = "white")


## (c) Body weight change (growth) ----

# Extract posterior draws
grow_post <- tidy_draws(grow_fit)

# Generate yrep for PPCs using posterior values for mean & variance
mu_mat <- grow_post %>% select(starts_with("mu[")) %>% as.matrix()
gamma_mat <- grow_post %>% select(starts_with("gamma[")) %>% as.matrix()
alpha_sigma_vec <- grow_post$alpha_sigma  

# Construct sigma matrix: sigma = exp(alpha_sigma + Z %*% gamma)
Z <- growth$W1 %>% as.matrix()
eta_sd_mat <- gamma_mat %*% t(Z)  # [draws x N] <- t(Z): [K x N]
alpha_mat <- matrix(alpha_sigma_vec, nrow = nrow(grow_post), ncol = nrow(growth))  # replicate across columns
sigma_mat <- exp(alpha_mat + eta_sd_mat)

# Simulate yrep ~ Normal(mu, sigma)
set.seed(1234)
yrep <- matrix(
  rnorm(length(mu_mat), mean = mu_mat, sd = sigma_mat),
  nrow = nrow(mu_mat)
)

# Function to undo unit scaling & transform back to original scale
deUnitScale <- function(scaled_x, sd_x, mean_x) {
  scaled_x * sd_x + mean_x
}

unscaled_y    <- deUnitScale(scaled_x = growth$W2, sd_x = growth$sd_W2, mean_x = growth$mean_W2)
unscaled_yrep <- apply(X = yrep, MARGIN = 1, FUN = deUnitScale, sd_x = growth$sd_W2, mean_x = growth$mean_W2) %>% t()

ppc_tot_mean   <- ppc_stat(y = unscaled_y, yrep = unscaled_yrep, stat = mean) &
  scale_x_continuous(n.breaks = 5) & 
  labs(x = "Mean weight at recapture (g)")

ppc_month_mean <- ppc_stat_grouped(y = unscaled_y, yrep = unscaled_yrep, group = month(growth$DATE1, label = TRUE), stat = mean) & 
  scale_x_continuous(limits = c(25, 50), n.breaks = 3) &
  labs(x = "Mean weight at recapture (g)")

ppc_tot_sd   <- ppc_stat(y = unscaled_y, yrep = unscaled_yrep, stat = sd) &
  scale_x_continuous(n.breaks = 5) & 
  labs(x = "SD weight at recapture (g)")

ppc_month_sd <- ppc_stat_grouped(y = unscaled_y, yrep = unscaled_yrep, group = month(growth$DATE1, label = TRUE), stat = sd) & 
  scale_x_continuous(limits = c(5, 18), n.breaks = 3) &
  labs(x = "SD weight at recapture (g)")

wrap_plots(list(ppc_tot_mean, ppc_month_mean, ppc_tot_sd, ppc_month_sd), nrow = 2) & 
  plot_annotation(tag_levels = 'A') & 
  theme_Publication()

# Save figure
ggsave("Plots/model_fits/PPC_growth.pdf", width = 10, height = 10, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/PPC_growth.png", width = 10, height = 10, units = "in", dpi = 1200, bg = "white")

# Compare predicted growth vs. observations by month
ppc_intervals_grouped(
  y = unscaled_y,
  yrep = unscaled_yrep,
  x = deUnitScale(scaled_x = growth$W1, sd_x = growth$sd_W1, mean_x = growth$mean_W1),
  group = month(growth$DATE1, label = TRUE),
  prob = 0.8,
  prob_outer = 0.95
) + 
  labs(x = "Weight at t (g)", y = "Weight at t+1 (g)") & 
  scale_x_continuous(limits = c(0, 85)) &
  scale_y_continuous(limits = c(0, 85)) &
  theme_Publication()

# Save figure
ggsave("Plots/model_fits/PPC_growth_month.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/PPC_growth_month.png", width = 8, height = 8, units = "in", dpi = 1200, bg = "white")


## (d) Male population size ----
y    <- MNA$mna.male
yrep <- posterior_predict(mf_fit)

ppc_tot_mean   <- ppc_stat(y = y, yrep = yrep, stat = mean) &
  scale_x_continuous(n.breaks = 5) & 
  labs(x = "Mean male population size")

ppc_month_mean <- ppc_stat_grouped(y = y, yrep = yrep, group = month(MNA$trap.date, label = TRUE), stat = mean) & 
  scale_x_continuous(limits = c(0, 170), n.breaks = 3) &
  labs(x = "Mean male population size")

ppc_tot_sd   <- ppc_stat(y = y, yrep = yrep, stat = sd) &
  scale_x_continuous(n.breaks = 5) & 
  labs(x = "SD male population size")

ppc_month_sd <- ppc_stat_grouped(y = y, yrep = yrep, group = month(MNA$trap.date, label = TRUE), stat = sd) & 
  scale_x_continuous(limits = c(0, 110), n.breaks = 3) &
  labs(x = "SD male population size")

wrap_plots(list(ppc_tot_mean, ppc_month_mean, ppc_tot_sd, ppc_month_sd), nrow = 2) & 
  plot_annotation(tag_levels = 'A') & 
  theme_Publication()

# Save figure
ggsave("Plots/model_fits/PPC_malepop.pdf", width = 10, height = 10, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/PPC_malepop.png", width = 10, height = 10, units = "in", dpi = 1200, bg = "white")

# Compare predicted male population size vs. observations by month
ppc_intervals_grouped(
  y = MNA$mna.male,
  yrep = yrep,
  x = MNA$mna.female,
  group = month(MNA$month, label = TRUE),
  prob = 0.8,
  prob_outer = 0.95
) + 
  labs(x = "Female population size", y = "Male population size") & 
  scale_y_continuous(limits = c(0, 500)) & 
  scale_x_continuous(limits = c(0, 400)) & 
  theme_Publication()

# Save figure
ggsave("Plots/model_fits/PPC_mf_month.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/PPC_mf_month.png", width = 8, height = 8, units = "in", dpi = 1200, bg = "white")


# 4. Plot posterior distributions ----

# Read in neat posterior distributions generated in model_fit_bayesian.R
all_posts <-readRDS("Outputs/stan_fits/neat_posteriors_alldat_28d.rds")

## (a) All posteriors ----

# Define custom parameter labels
par_labs <- c(
  "alpha_mu" = expression(alpha[mu]),
  "alpha_sigma" = expression(alpha[sigma]),
  "alpha" = expression(alpha),
  "beta_weight" = expression(beta[w]),
  "beta_pop" = expression(beta[pop]),
  "beta_tseas" = expression(beta[tseas]),
  "beta_tvar" = expression(beta[tvar]),
  "beta_ttrend" = expression(beta[ttrend]),
  "beta_pseas" = expression(beta[pseas]),
  "beta_pvar" = expression(beta[pvar]),
  "beta_ptrend" = expression(beta[ptrend]),
  "beta_sigma"= expression(beta[sigma])
)

# Define the desired parameter order (must match names in your data)
param_order <- c(
  "alpha",
  "alpha_mu",
  "alpha_sigma",
  "beta_sigma",
  "beta_weight",
  "beta_pop",
  "beta_ttrend",
  "beta_tseas",
  "beta_tvar",
  "beta_ptrend",
  "beta_pseas",
  "beta_pvar"
)

# Format posteriors
all_posts_plot <- all_posts %>%
  filter(process != "male_pop") %>% 
  select(!contains("lag")) %>% 
  distinct() %>% 
  pivot_longer(cols = c(contains("alpha"), starts_with("beta")), names_to = "parameter") %>%
  filter(!is.na(value)) %>%
  mutate(parameter = factor(parameter, levels = param_order)) %>% 
  # Create nice groupings for plot colours
  mutate(
    param_group = 
      case_when(
        str_detect(parameter, "_t") ~ "Temperature",
        str_detect(parameter, "pseas|pvar|ptrend") ~ "Precipitation",
        str_detect(parameter, "pop") ~ "Population size",
        str_detect(parameter, "weight") | str_detect(parameter, "beta_sigma") ~ "Weight",
        str_detect(parameter, "alpha") ~ "Intercept"
      ),
    process = case_when(
      str_detect(process, "rec")  ~ "Recruitment",
      str_detect(process, "surv") ~ "Survival",
      str_detect(process, "grow") ~ "Body weight change",
    )
  )

# Define manual colour palette
pal <- c(
  "Temperature"     = "darkred",
  "Precipitation"   = "lightblue", 
  "Population size" = "orange",
  "Weight"          = "forestgreen",
  "Intercept"       = "grey"
)

# Summarise posteriors
posterior_summary <- all_posts_plot %>%
  group_by(process, parameter) %>%
  summarise(
    median = median(value),
    lower80 = quantile(value, 0.10),
    upper80 = quantile(value, 0.90),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975),
    .groups = "drop"
  )

# Plot all posteriors
all_posts_plot %>%
  ggplot(aes(x = value, y = fct_rev(parameter))) +
  geom_density_ridges(aes(fill = param_group), scale = 1.5, alpha = .8, col = NA) +
  geom_segment(  # Outer interval (thin)
    data = posterior_summary,
    aes(x = lower95, xend = upper95, y = parameter, yend = parameter), linewidth = 0.4, inherit.aes = FALSE
  ) +
  geom_segment(  # Inner interval (thicker)
    data = posterior_summary,
    aes(x = lower80, xend = upper80, y = parameter, yend = parameter), linewidth = 0.8, inherit.aes = FALSE
  ) +
  geom_point(
    data = posterior_summary,
    aes(x = median, y = parameter), size = 1.5, inherit.aes = FALSE
  ) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  facet_wrap(~process, nrow = 1, scales = "free") +
  scale_y_discrete(labels = par_labs) +
  labs(x = "Scaled posterior parameter value",
       y = "Parameter",
       fill = "Parameter group") +
  scale_fill_manual(values = pal) +
  theme_Publication()

# Save plot
ggsave("Plots/model_fits/all_posteriors.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)
ggsave("Plots/model_fits/all_posteriors.png", width = 8, height = 8, units = "in", dpi = 1200, bg = "white")


## (b) Climate lags ----

# Function to extract lag weights from stanfit object
extractLagWeights <- function(model_fit, lag_grid) {
  
  # Extract posterior
  posterior <- tidy_draws(model_fit)
  
  # Name lag grid so can join to lag weight posteriors
  lag_grid <- lag_grid %>% mutate(lag_name = paste0("lag_weights[", row_number(), "]"))
  
  # Extract lag weights
  lag_weights <- posterior %>%
    select(draw = .draw, starts_with("lag_weights")) %>%
    pivot_longer(cols = !draw, names_to = "lag", values_to = "weight") %>%
    # Join lag values
    left_join(lag_grid, by = c("lag" = "lag_name")) %>%
    mutate(model = model_fit@model_name)
  
  return(lag_weights)
}

# Generate grid of climate lags
lag_int  <- 28
lag_seq  <- seq(0, 168, by = lag_int)
lag_grid <- setNames(expand.grid(lag_seq, lag_seq), c("t_lag", "p_lag")) %>% as.data.frame()

# Extract lag weights from all demographic process models
surv_fit@model_name <- "Survival"; rec_fit@model_name <- "Pregnancy"; grow_fit@model_name <- "Growth"
lag_weights <- lapply(X = list(surv_fit, rec_fit, grow_fit), FUN = extractLagWeights, lag_grid = lag_grid) %>% bind_rows()

# Calculate medians
lag_medians <- lag_weights %>%
  group_by(model, t_lag, p_lag) %>% 
  mutate(med_weight = median(weight)) %>% 
  ungroup() %>% 
  distinct()

# Plot heat map & bar plot of posterior weight for each model
wrap_plots(
  lag_medians %>% 
    filter(model %in% c("Growth", "Pregnancy")) %>% 
    select(model, p_lag, med_weight) %>% distinct() %>% 
    ggplot() + 
    geom_col(aes(x = p_lag, y = med_weight)) + 
    facet_wrap(~model) + 
    labs(x = "Precipitation lag (days)", y = "Median posterior weight"), 
  lag_medians %>% 
    filter(model == "Survival") %>% 
    ggplot(aes(x = p_lag, y = t_lag, fill = med_weight)) +
    geom_tile() +
    facet_wrap(~model) +
    scale_fill_viridis_c() +
    labs(x = "Precipitation lag (days)", y = "Temperature lag (days)", fill = "Median posterior weight"), 
  nrow = 2
) & 
  theme_Publication() & 
  theme(legend.position = "bottom", legend.key.width = unit(1, "null"))

# Save plot
ggsave("Plots/model_fits/lag_posteriors_allmods.pdf", width = 8, height = 10, units = "in")
ggsave("Plots/model_fits/lag_posteriors_allmods.png", width = 8, height = 10, units = "in", dpi = 1200, bg = "white")

# Plot posterior lag weights for survival
lag_weights %>%
  filter(model == "Survival") %>% 
  # Calculate median weight for each lag combo
  group_by(model, lag) %>%
  mutate(med_weight = median(weight)) %>%
  ggplot() +
  geom_density_ridges(aes(x = weight, y = as.factor(t_lag), scale = 3,  fill = (med_weight)), alpha = .8, col = NA) +
  facet_grid(p_lag ~ model) +
  labs(x = "Posterior weight", y = "Temperature lag (days)", fill = "Median posterior weight") +
  scale_fill_viridis_c() +
  theme(
    axis.title.y  = element_text(vjust = 2.5),
    plot.margin   = unit(rep(0.5, 4), "cm"),
    panel.spacing = unit(0.5, "lines"),
    strip.text.x  = element_text(face = "bold"),
    strip.text.y  = element_text(face = "bold", angle = 0),
    legend.position = "bottom",
    legend.key.width = unit(1, "null")
  )

# Save plot
ggsave(filename = "Plots/model_fits/lag_posteriors_surv.pdf", width = 8, height = 10, units = "in")
ggsave(filename = "Plots/model_fits/lag_posteriors_surv.png", width = 8, height = 10, units = "in", dpi = 1200, bg = "white")
