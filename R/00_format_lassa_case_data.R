#------------------------------------------#
# Load & format Lassa case data in Nigeria #
#------------------------------------------#

# 0. Script purpose ----
# - Loads Lassa case data compiled from situation reports by Nigerian CDC
# - Harmonises & formats data & saves as .rds file for subsequent analyses in other scripts


# 1. Load packages ----
pacman::p_load(dplyr, tidyr, ggplot2, stringr, lubridate, zoo, gganimate)


# 2. Load raw Lassa case data ----

# Jan 2018 to Apr 2022
case.data.early <- read.csv("Data/weekly_cases_bystate.csv") %>% mutate(week_commencing = dmy(week_commencing))

# Feb 2020 to Aug 2025
case.data.late <- read.csv("Data/lassa_data_latest.csv") %>% 
  mutate(year = 2000 + year) %>% 
  reframe(year, epiweek = week, state = states, con_cases = confirmed)

library(ISOweek)

# build ISO week string (e.g. "2025-W32-1" = Monday of week 32, 2025)
case.data.late$week_commencing <- ISOweek2date(
  paste0(case.data.late$year, "-W", sprintf("%02d", case.data.late$epiweek), "-1")
)

# # Remove data found in latest dataset
# case.data.early <- case.data.early %>% filter(week_commencing < min(case.data.late$week_commencing))


# 3. Format & harmonise data ----
  
# Define column names to keep
shared.cols <- c("year", "epiweek", "week_commencing", "state", "con_cases")

# Define relevant Nigerian states to keep
states <- c("Bauchi", "Ebonyi", "Edo", "Ondo", "Taraba")

# Format Jan '18 - Apr '22 case data
formatted.case.data.early <- case.data.early %>%
  # Select relevant states & put all confirmed cases in same column (long format)
  pivot_longer(cols = all_of(states), names_to = "state", values_to = "con_cases") %>% 
  mutate(
    # Make NA cases = 0
    con_cases = ifelse(is.na(con_cases), 0, con_cases)
    ) %>% 
  # select columns
  dplyr::select(all_of(shared.cols))

# Format Feb '20 - Aug '25 case data
formatted.case.data.late <- case.data.late %>%
  filter(state %in% states) %>% 
  arrange(state, week_commencing) %>% 
  # Ensure weeks that are missing are represented as no reported cases for given state
  group_by(year, state) %>%
  complete(epiweek=full_seq(epiweek, 1)) %>% 
  mutate(con_cases = if_else(is.na(con_cases), 0, con_cases))

# Recalculate week commencing for previously missing epiweeks
formatted.case.data.late$week_commencing <- ISOweek2date(
  paste0(formatted.case.data.late$year, "-W", sprintf("%02d", formatted.case.data.late$epiweek), "-1")
)

# Exclude data overlapping with earlier time series
formatted.case.data.late <- formatted.case.data.late %>% 
  filter(week_commencing > max(formatted.case.data.early$week_commencing))



# 4. Combine & save data ----

# Combine all case data (early + late)
combined.case.data <- full_join(formatted.case.data.early, formatted.case.data.late) %>%
  # Create date mid-point
  mutate(date_mid = week_commencing + 3) %>%
  arrange(state, week_commencing)

# Save combined case data
saveRDS(combined.case.data, "Data/combined_case_data.rds")


# 5. Make plots ----

doplot <- 1

# Plot confirmed cases
if (doplot) {
  
  combined.case.data %>% 
    ggplot(aes(date_mid, con_cases)) +
    geom_col() +
    facet_wrap(~state, ncol = 1) + 
    labs(x = "", y = "Weekly confirmed cases") +
    theme_minimal() + 
    theme(text = element_text(size = 12, family = "serif"))
  
  # Save plot
  ggsave("Plots/case_data/confirmed_cases.pdf", width = 6, height = 8, units = "in", device = cairo_pdf)
  ggsave("Plots/case_data/confirmed_cases.png", width = 6, height = 8, units = "in", dpi = 1200, bg = "white")
}

# Plot relative cases
if (doplot) {
  
  combined.case.data %>% 
    ggplot(aes(date_mid, rel_cases)) +
    geom_col() +
    facet_wrap(~state, nrow = 1) + 
    labs(x = "Date", y = "Weekly relative confirmed cases") +
    theme_minimal() + 
    theme(text = element_text(size = 14, family = "serif"))
  
}

# Plot confirmed & relative cases on same figure
if (doplot) {
  
  combined.case.data %>%
    pivot_longer(contains("cases")) %>% 
    mutate(name = case_when(str_detect(name, "con_") ~ "Confirmed cases", .default = "Relative cases")) %>% 
    ggplot(aes(date_mid, value)) +
    geom_col() +
    facet_grid(name ~ state, scales = "free_y", switch = "y", axis.labels = "all_y") +
    labs(x = "Date", y = "") +
    theme_minimal() + 
    theme(text = element_text(size = 14, family = "serif"))
  
  # Save plot
  ggsave("Plots/case_data/confirmed_and_relative_cases.pdf", width = 12, height = 6, units = "in", device = cairo_pdf)
}
