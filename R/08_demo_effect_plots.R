#----------------------------------------------------------------#
# Plot climate and demographic effects on demographic parameters #
#----------------------------------------------------------------#

# 0. Script purpose ----
# - Read in posterior predictions
# - Calculate contributions of climate & demographic effects to demographic processes
# - Make & save plots


# 1. Load packages ----
pacman::p_load(purrr, patchwork, ggplot2, ggridges, tidybayes, bayesplot, forcats, foreach)
source("R/ggplot_theme.R")  # Source custom ggplot theme
use_cluster <- FALSE  # Use the high performance cluster?
source("R/setup_IPM.R")


# 2. Load simulated output ----

# Read in & combine .rds files
preds <- readRDS("Outputs/posterior_predictions_N=1000.rds")
posterior_samples <- readRDS("Outputs/posterior_predict_samples_N=1000.rds")
n_samp <- n_distinct(preds$iter)
subsample_indices <- 1:nrow(posterior_samples)

if(max(subsample_indices) > n_samp) subsample_indices <- subsample_indices[1:n_samp]

pop.ajd  <- mean(cmr.data$POP) / (n_init * 2)
mean.wgt <- mean(preds$mean_adult_wgt)


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


# 3. Compute effect contributions from different covariates ----

# Effects of:
# - temperature (seas,var,trend)
# - precipitation (seas,var,trend)
# - populatin size
# - weight
# On:
# - growth
# - recruitment
# - survival

# Define function for getting effects
get_effects_over_time <- function(j, sampled_inf_post, demog_posts, preds){
  
  simulation <- preds %>% filter(iter==subsample_indices[j])
  
  dates <- lag_clim %>%
    filter(date %in% simulation$date) %>%
    pull(date)
  
  simulation <- simulation %>% filter(date %in% dates) %>% arrange(date)
  clim_df <- lag_clim %>% filter(date %in% dates) %>% arrange(date)
  
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
  
  posts_surv <- sampled_demog_posts %>% filter(process=="surv")
  posts_rec <- sampled_demog_posts %>% filter(process=="rec")
  posts_grow <- sampled_demog_posts %>% filter(process=="grow")
  
  # fn to compute lagged clim
  lagged_clim <- function(posts){
    lagged_env <- inner_join(clim_df,
                             posts %>% select(lag_wgt_val, t_lag, p_lag),
                             by = join_by(t_lag,p_lag)) %>%
      group_by(date) %>%
      summarise(ttrend = sum(ttrend*lag_wgt_val),
                tseas = sum(tseas*lag_wgt_val),
                tvar = sum(tvar*lag_wgt_val),
                ptrend = sum(ptrend*lag_wgt_val),
                pseas = sum(pseas*lag_wgt_val),
                pvar = sum(pvar*lag_wgt_val)) %>%
      as.data.frame() %>%
      arrange(date)
    
    return(lagged_env)
  }
  
  # Function to compute effects for a given process
  get_process_effects <- function(chosen_process){
    
    process_posts <- sampled_demog_posts %>% filter(process==chosen_process)
    process_lagged_env <- lagged_clim(process_posts)
    
    process_pars <- sampled_demog_posts %>%
      filter(process==chosen_process) %>%
      select(process,alpha,alpha_mu,beta_weight,beta_pop,beta_ttrend,beta_tseas,beta_tvar,beta_ptrend,beta_pseas,beta_pvar) %>%
      unique()
    
    if(chosen_process=="grow"){
      process_pars <- process_pars %>%
        mutate(alpha = alpha_mu)
    }
    
    process_pars <- process_pars %>% select(-alpha_mu)
    
    if(chosen_process=="grow"){
      # divide body weight response by mean simulated body weight
      effect.mod <- 1/mean.wgt
    }else{
      effect.mod <- 1
    }
    
    effect <- data.frame(
      iter = j,
      date = simulation$date,
      process = chosen_process,
      alpha = effect.mod * process_pars$alpha,
      env_seas = effect.mod * (process_lagged_env$tseas * ifelse(is.na(process_pars$beta_tseas), 0, process_pars$beta_tseas) +
                                 process_lagged_env$pseas * ifelse(is.na(process_pars$beta_pseas), 0, process_pars$beta_pseas)),
      env_var = effect.mod * (process_lagged_env$tvar * ifelse(is.na(process_pars$beta_tvar), 0, process_pars$beta_tvar) +
                                process_lagged_env$pvar * ifelse(is.na(process_pars$beta_pvar), 0, process_pars$beta_pvar)),
      env_trend = effect.mod * (process_lagged_env$ttrend * ifelse(is.na(process_pars$beta_ttrend), 0, process_pars$beta_ttrend) +
                                  process_lagged_env$ptrend * ifelse(is.na(process_pars$beta_ptrend), 0, process_pars$beta_ptrend)),
      pop = effect.mod * simulation$mod.pop * pop.ajd * ifelse(is.na(process_pars$beta_pop), 0, process_pars$beta_pop),
      weight = effect.mod * simulation$mean_adult_wgt * ifelse(is.na(process_pars$beta_weight), 0, process_pars$beta_weight)
    )
    
    return(effect)
  }
  
  effects <- rbind(get_process_effects("surv"),
                   get_process_effects("rec"),
                   get_process_effects("grow"))
  
  lagged_env_surv <- lagged_clim(posts_surv)
  lagged_env_rec <- lagged_clim(posts_rec)
  lagged_env_grow <- lagged_clim(posts_grow)
  
  demog_pars <- sampled_demog_posts %>% 
    filter(process %in% c("surv","rec","grow")) %>%
    select(process,
           alpha, alpha_mu, alpha_sigma,
           beta_weight, beta_pop, beta_sigma,
           beta_ttrend, beta_tseas, beta_tvar, beta_ptrend, beta_pseas, beta_pvar) %>%
    unique()
  
  pars_surv <- demog_pars %>% filter(process=="surv")
  pars_rec <- demog_pars %>% filter(process=="rec")
  pars_grow <- demog_pars %>% filter(process=="grow")
  
  return(effects)
}

if (!file.exists("Outputs/demographic_effects_N=1000.rds")) {
  
  # Initialise parallel backend before foreach() loop
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores)
  
  # Estimate climate effects for all samples
  effects <- foreach(j = 1:1000, .packages = c("dplyr", "truncnorm", "lubridate")) %dopar% {
    get_effects_over_time(j = j, sampled_inf_post = posterior_samples, demog_posts = posts,
                          preds = preds)
  }
  
  # Stop and clean up the parallel cluster
  doParallel::stopImplicitCluster()
  
  # Save predictions
  effects <- bind_rows(effects)
  saveRDS(effects, "Outputs/demographic_effects_N=1000.rds")
  
} else {
  effects <- readRDS("Outputs/demographic_effects_N=1000.rds")
}


# 4. Wrangle effects ----

effect_names <- c("alpha","env_trend","env_seas","env_var","pop","weight")

## ACROSS ALL SAMPLES

# get median values across samples
effects_median <- effects %>%
  pivot_longer(cols = -c("date","iter","process"), names_to = "effect") %>%
  group_by(date,process,effect) %>%
  summarise(value=median(value), .groups="keep") %>%
  ungroup() %>%
  pivot_wider(names_from = "effect")

# get marginal effects (all other covariates held at their mean value)
effects_marginal <- effects_median %>%
  pivot_longer(cols = -c("date","process","alpha"), names_to = "beta_effect") %>%
  group_by(process,beta_effect) %>%
  mutate(meanvalue = mean(value),
         value = value-meanvalue) %>%
  ungroup() %>%
  group_by(date,process) %>%
  mutate(alpha = alpha + sum(meanvalue)) %>%
  ungroup() %>%
  select(-meanvalue)

# Get effects relative to some baseline value?
effects_baselined <- effects_median %>%
  pivot_longer(cols = -c("date","process","alpha"), names_to="beta_effect") %>%
  group_by(process,beta_effect) %>%
  mutate(minvalue = min(value),
         value = value-minvalue) %>%
  ungroup() %>%
  group_by(date,process) %>%
  mutate(alpha = alpha + sum(minvalue)) %>%
  ungroup() %>%
  select(-minvalue)

kappa <- median(surv.adj)  # plot survival probability
# kappa <- 0  # plot recapture probability

# Transform logit-space into probabilities for recruitment & survival
effects_prob <- effects_baselined %>%
  filter(process%in%c("rec","surv")) %>%
  group_by(date,process) %>%
  mutate(effects_total = sum(value),
         p_baseline = 1/(1+exp(-alpha)),
         p_total = 1/(1+exp(-alpha-effects_total)),
         p_effects = p_total-p_baseline,
         alpha = p_baseline,
         value = (value / effects_total) * p_effects,
         alpha = ifelse(process=="surv", kappa + (1-kappa)*alpha, alpha),
         value = ifelse(process=="surv", (1-kappa)*value, value)) %>%
  ungroup() %>%
  select(date,process,alpha,beta_effect,value)

# Combine all output
effects_tf <- effects_baselined %>%
  filter(process=="grow") %>%
  rbind(effects_prob)


## GROUPED BY POSTERIOR SAMPLE

# baseline per sample as well as per time step
effects_baselined_per_sample <- effects %>%
  pivot_longer(cols = -c("iter","date","process","alpha"), names_to = "beta_effect") %>%
  group_by(iter,process,beta_effect) %>%
  mutate(minvalue = min(value),
         value = value-minvalue) %>%
  ungroup() %>%
  group_by(iter,date,process) %>%
  mutate(alpha = alpha + sum(minvalue)) %>%
  ungroup() %>%
  select(-minvalue)

# get probabilities for the above baselined effects
effects_prob_per_sample <- effects_baselined_per_sample %>%
  filter(process %in% c("surv","rec")) %>%
  group_by(iter,date,process) %>%
  mutate(effects_total = sum(value),
         p_baseline = 1/(1+exp(-alpha)),
         p_total = 1/(1+exp(-alpha-effects_total)),
         p_effects = p_total-p_baseline,
         alpha = p_baseline,
         value = (value / effects_total) * p_total,
         alpha = ifelse(process=="surv", kappa + (1-kappa)*alpha, alpha),
         value = ifelse(process=="surv", (1-kappa)*value, value)) %>%
  ungroup() %>%
  select(iter,date,process,alpha,beta_effect,value)

# combine probabilities (surv and rec) with growth
effects_tf_per_sample <- effects_baselined_per_sample %>%
  filter(process=="grow") %>%
  rbind(effects_prob_per_sample)

# converting the effects into a proportion of overall response (still per sample and per time step)
effects_proportion_per_sample <- effects_tf_per_sample %>%
  pivot_wider(names_from = "beta_effect") %>%
  pivot_longer(cols = -c("iter","date","process"), names_to="effect") %>%
  group_by(iter,date,process) %>%
  mutate(value=value/sum(value))


# 5. Visualise effects ----

# Define effect names in display order
effect_names <- c(
  "Baseline",
  "Climate trend",
  "Climate seasonality",
  "Climate variability",
  "Population size",
  "Weight"
)

# Generate palette
palette <- viridis::viridis(length(effect_names), option = "D")
names(palette) <- effect_names

## MEDIAN EFFECTS

# Plot contributions of each covariate to each demographic process
stacked_plot <- effects_tf %>%
  pivot_wider(names_from="beta_effect") %>%
  pivot_longer(cols = -c("date","process"), names_to="effect") %>%
  mutate(process = case_when(process == "grow" ~ "Body weight change",
                             process == "rec"  ~ "Recruitment",
                             process == "surv" ~ "Survival"), 
         effect = case_when(effect == "alpha"     ~ "Baseline",
                            effect == "env_trend" ~ "Climate trend",
                            effect == "env_seas"  ~ "Climate seasonality",
                            effect == "env_var"   ~ "Climate variability",
                            effect == "pop"       ~ "Population size",
                            effect == "weight"    ~ "Weight")) %>% 
  mutate(effect = factor(effect,levels=rev(effect_names))) %>%
  ggplot() +
  geom_col(aes(x=date,y=value,fill=effect), width = 28, col = NA, position="stack") +
  geom_hline(data = data.frame(y=1,process="Body weight change"), aes(yintercept=y), linetype="dashed") +
  facet_wrap(~process, scales="free") + 
  theme_Publication() + 
  scale_fill_manual(values = palette) +
  labs(x = "Date", y = "Estimate", fill = "Effect", col = "Effect")


## EFFECTS OVER WHOLE POSTERIOR

# Summarise data for whole posterior
whole_post_dat <- effects_proportion_per_sample %>% 
  filter(effect != "alpha") %>% 
  pivot_wider(names_from="effect") %>%
  pivot_longer(cols = -c("date","process"), names_to="effect") %>%
  mutate(process = case_when(process == "grow" ~ "Body weight change",
                             process == "rec"  ~ "Recruitment",
                             process == "surv" ~ "Survival"), 
         effect = case_when(effect == "alpha"     ~ "Baseline",
                            effect == "env_trend" ~ "Climate trend",
                            effect == "env_seas"  ~ "Climate seasonality",
                            effect == "env_var"   ~ "Climate variability",
                            effect == "pop"       ~ "Population size",
                            effect == "weight"    ~ "Weight")) %>% 
  mutate(effect = factor(effect,levels=rev(effect_names))) %>%
  filter(!is.na(effect)) %>% # remove effects w/ non-matching names (e.g. baseline)
  group_by(process, effect) %>% 
  filter(!all(value == 0)) %>%  # remove estimates for effects not present in given model
  mutate(year = year(date))

# Calculate credible intervals
posterior_summary <- whole_post_dat %>%
  group_by(process, effect) %>%
  summarise(
    median = median(value),
    lower80 = quantile(value, 0.10),
    upper80 = quantile(value, 0.90),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975)
  )

# Proportion of response explained by each covariate
prop_plot <- posterior_summary %>% 
  ggplot(aes(col = effect)) +
  geom_linerange(  # Outer interval (thin)
    aes(x = effect, ymin = lower95, ymax = upper95),
    size = 0.4,
  ) +
  geom_linerange(  # Inner interval (thicker)
    aes(x = effect, ymin = lower80, ymax = upper80),
    size = 0.6,
  ) +
  geom_point(aes(x = effect, y = median), size = 1.5) +
  facet_wrap(~process,scales="free") + 
  theme_Publication() + 
  scale_color_manual(values = palette) +
  labs(x = "Effect", y = "Contribution to response (proportion)", col = "Effect") + 
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Plot as time series summarising by year
years <- unique(whole_post_dat$year)
odd_years <- years[c(TRUE, FALSE)]
shading <- data.frame(
  xmin = odd_years - 0.5,
  xmax = odd_years + 0.5,
  ymin = -Inf,
  ymax = Inf
)

posterior_summary <- whole_post_dat %>%  # Calculate credible intervals
  group_by(process, effect, year) %>%
  summarise(
    median = median(value),
    lower80 = quantile(value, 0.10),
    upper80 = quantile(value, 0.90),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975)
  )
# View(posterior_summary %>% select(!contains("80"))) # summary of predictions + 95% credible intervals

prop_plot_time <- posterior_summary %>% 
  ggplot(aes(group = interaction(year, effect), col = effect)) +
  geom_rect(data = shading,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "gray90") +
  geom_linerange(  # Outer interval (thin)
    aes(x = year, ymin = lower95, ymax = upper95),
    linewidth = 0.4, position = position_dodge(width = .8)
  ) +
  geom_linerange(  # Inner interval (thicker)
    aes(x = year, ymin = lower80, ymax = upper80),
    linewidth = 0.6, position = position_dodge(width = .8)
  ) +
  geom_point(aes(x = year, y = median), size = 1.5, position = position_dodge(width = .8)) +
  facet_wrap(~process, scales="free") + 
  theme_Publication() + 
  scale_color_manual(values = palette) +
  labs(x = "Year", y = "Contribution to response (proportion)", col = "Effect")

# Panel plot of proportion contributing to response
wrap_plots(prop_plot, prop_plot_time, nrow = 2) & plot_annotation(tag_levels = "A")

# Save plot
ggsave("Plots/model_fits/process_model_prop_panel.pdf", device = cairo_pdf, height = 8, width = 8, units = "in")
ggsave("Plots/model_fits/process_model_prop_panel.png", bg = "white", dpi = 600, height = 8, width = 8, units = "in")

# Calculate credible intervals (e.g. for a table)
effects_full_credible_interval <- effects_tf_per_sample %>%
  pivot_wider(names_from = "beta_effect") %>%
  pivot_longer(cols = -c("iter","date","process"), names_to="effect") %>%
  mutate(effect = factor(effect, levels=rev(effect_names))) %>%
  group_by(process,effect) %>%
  summarise(lower=quantile(value,0.025),
            median=quantile(value,0.5),
            upper=quantile(value,0.975),
            .groups="keep")

effects_proportion_full_credible_interval <- effects_proportion_per_sample %>%
  group_by(process, effect) %>% 
  filter(effect != "alpha", !all(value == 0)) %>%  # remove estimates for effects not present in given model
  mutate(
    median=round(quantile(value,0.5)*100, 2),
    lower=round(quantile(value,0.025)*100, 2),
    upper=round(quantile(value,0.975)*100, 2)
  ) %>% 
  reframe(estimate=paste0(median, " (", lower, "–", upper, ")")) %>% 
  distinct() %>% 
  arrange(effect)
print(effects_proportion_full_credible_interval)


## Plot posterior distributions

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

# Define the desired parameter order (must match names in the data)
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

# Read in & clean demographic process posteriors
all_posts <- readRDS("Outputs/stan_fits/neat_posteriors_alldat_28d.rds")

all_posts_plot <- all_posts %>%
  filter(process != "male_pop") %>%
  select(!contains("lag")) %>%
  distinct() %>%
  pivot_longer(cols = c(contains("alpha"), starts_with("beta")), names_to = "parameter") %>%
  filter(!is.na(value)) %>%
  mutate(parameter = factor(parameter, levels = param_order)) %>%
  mutate(
    process = case_when(
      str_detect(process, "rec")  ~ "Recruitment",
      str_detect(process, "surv") ~ "Survival",
      str_detect(process, "grow") ~ "Body weight change",
    )
  )

# Map to same categories for consistent colours
all_posts_plot <- all_posts_plot %>%
  mutate(effect = case_when(
    parameter == "beta_ttrend" | parameter == "beta_ptrend" ~ "Climate trend",
    parameter == "beta_tseas" | parameter == "beta_pseas"  ~ "Climate seasonality",
    parameter == "beta_tvar" | parameter == "beta_pvar"  ~ "Climate variability",
    parameter == "beta_pop"    ~ "Population size",
    parameter == "beta_weight" ~ "Weight"
  ))

# Summarise posteriors
posterior_summary <- all_posts_plot %>%
  filter(!str_detect(parameter, "alpha|beta_sigma")) %>%
  group_by(process, parameter) %>%
  summarise(
    median = median(value),
    lower80 = quantile(value, 0.10),
    upper80 = quantile(value, 0.90),
    lower95 = quantile(value, 0.025),
    upper95 = quantile(value, 0.975),
    .groups = "drop"
  )

# Plot posteriors
posterior_plot <- all_posts_plot %>%
  filter(!str_detect(parameter, "alpha|beta_sigma")) %>%
  ggplot(aes(x = value, y = fct_rev(parameter))) +
  geom_density_ridges(aes(fill = effect), scale = 1.5, alpha = .8, col = NA) +
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
  scale_fill_manual(values = palette) +
  theme_Publication() +
  theme(legend.position = "none")

# Panel plot of posteriors & demographic/climate effects
wrap_plots(posterior_plot, stacked_plot, nrow = 2) & plot_annotation(tag_levels = "A")

# Save plot
ggsave("Plots/model_fits/process_model_panel.pdf", device = cairo_pdf, height = 8, width = 8, units = "in")
ggsave("Plots/model_fits/process_model_panel.png", bg = "white", dpi = 600, height = 8, width = 8, units = "in")
