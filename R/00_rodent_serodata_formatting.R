#------------------------------------------------------#
# Load & format M. natalensis Morogoro arenavirus data #
#------------------------------------------------------#

# 0. Script purpose ----
# - Loads Mastmoys natalensis Morogoro arenavirus seroprevalence data from MOSA CMR study
# - Formats data & saves as .rds file for subsequent analyses in other scripts


# 1. Load packages ----
pacman::p_load(readr, dplyr, tidyr, ggplot2, stringr, lubridate, lhs, ggridges, 
               foreach, doParallel)


# 2. Load data ----

# Morogoro arenavirus data from MOSA (Mariën et al., 2020. J Animal Ecol, 2)
ind.seroprev.dat <- readr::read_delim("Data/MOSA FINAL data Nov_2017.txt", delim = "\t", show_col_types = FALSE)


# 3. Format seroprevalence data ----

ind.seroprev.dat <- ind.seroprev.dat %>%
  mutate(DATE = dmy(DATE)) %>%
  filter(
    year(DATE) >= 2010,
    # Get only Mastomys natalensis
    SPECIES == "MN", 
    # Remove 0g weights
    WEIGHT != 0,
    # Remove NAs:
    !is.na(DATE),   # dates
    !is.na(WEIGHT), # weights
    !is.na(TOE),    # ID
    !is.na(SEX_COND),  # sexual condition
    !is.na(SEX),
    !is.na(IFA)  # Ab test
  ) %>% 
  mutate(
    # Define sexually active
    FER = case_when(
      # Males: sexually active if both testes are scrotal ("S") and gubernacula are visible ("V")
      SEX == "M" & str_detect(SEX_COND, "S") & str_detect(SEX_COND, "V") ~ TRUE,
      SEX == "M" ~ FALSE,
      # Females: sexually active if vagina is perforated ("P"), nipples lactating ("L"), or visibly pregnant ("Y")
      SEX == "F" & (str_detect(SEX_COND, "P") | str_detect(SEX_COND, "L") | str_detect(SEX_COND, "Y")) ~ TRUE,
      SEX == "F" ~ FALSE
    ),
    # Define antibody-positive & negative (don't classify uncertain as either)
    ab_pos = case_when(
      str_detect(IFA, "neg") ~ 0,
      str_detect(IFA, "pos") ~ 1
    )
  ) %>% 
  filter(!is.na(ab_pos)) %>% 
  # Select relevant columns
  dplyr::select(DATE, SEX, WEIGHT, TOE, FER, ab_pos) %>% 
  ### Create a trapping sessions group column
  # Arrange data by date
  arrange(DATE) %>%
  # Group sessions: if lagged difference > 3, assume different date session
  mutate(session = cumsum(c(1, diff(DATE) > 3))) %>%
  # Create session start date as the first date in each session group
  group_by(session) %>%
  mutate(session.start = first(DATE),
         date_mid = mean(session.start)) %>%
  ungroup() %>%
  # Drop session group (only need session start date as grouping factor)
  dplyr::select(!session)

# Group seroprevalence data by trap session, calculate k, n, seroprevalence
seroprev.dat <- ind.seroprev.dat %>% 
  group_by(session.start, SEX, FER) %>% 
  reframe(date_mid,
          session.start,
          seropos = sum(ab_pos),
          tested = n(),
          seroprev = seropos / tested,
          # Calculate binomial 95% CIs
          mean.seroprev  = Hmisc::binconf(seropos, tested, method = "wilson")[,1],
          lower.seroprev = Hmisc::binconf(seropos, tested, method = "wilson")[,2],
          upper.seroprev = Hmisc::binconf(seropos, tested, method = "wilson")[,3]) %>% 
  ungroup() %>% 
  distinct() %>% 
  # Get year & month of sampling
  mutate(year = year(session.start), month = month(session.start))

## No difference in seroprevalence by sex when adjusting for year & month
# lm(seroprev ~ month + year + SEX + FER, data = seroprev.dat) %>% summary()

# So recalculate seroprevalence with sexes combined but still separate by subadult/adult
seroprev.dat <- ind.seroprev.dat %>% 
  group_by(session.start, FER) %>% 
  reframe(date_mid,
          session.start,
          seropos = sum(ab_pos),
          tested = n(),
          seroprev = seropos / tested,
          # Calculate binomial 95% CIs
          mean.seroprev  = Hmisc::binconf(seropos, tested, method = "wilson")[,1],
          lower.seroprev = Hmisc::binconf(seropos, tested, method = "wilson")[,2],
          upper.seroprev = Hmisc::binconf(seropos, tested, method = "wilson")[,3]) %>% 
  ungroup() %>% 
  distinct()


# 4. Save formatted seroprevalence data ----
saveRDS(seroprev.dat, "Data/MORV_seroprev_formatted.rds")
