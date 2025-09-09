#----------------------------------------------------------------------------#
# Convert capture data into inputs to infer survival, growth and recruitment #
#----------------------------------------------------------------------------#

# 0. Session details ----

# Script purpose:
# - Read in the raw Mastomys natalensis capture-mark-recapture data
# - Convert data into inputs to infer survival, growth and recruitment
# - Save these data for later analyses


# 1. Load packages ----
pacman::p_load(dplyr, ggplot2, stringr, lubridate, purrr, tidyr)


# Set toggle (0/1; off/on) for plotting throughout the script 
doplot <- 0


# 2. Load & format trapping data ----

# Read in & format CMR data (Leirs et al, 2023. Scientific Data, 10(1): 798)
mosa <- read.csv("Data/mosa_Feb2023_flags.csv", sep = ";") %>%
  mutate(
    DATE = dmy(DATE),
    # Replace odd ID number
    TOE = case_when(
      NUMBER == 55881 ~ str_replace(string = TOE, pattern = "\xa3", replacement = "£"),
      .default = TOE
    ),
    # Create flag for OK IDs (those that contain any numbers)
    id.ok = grepl("[[:digit:]]", TOE)
  ) %>%
  filter(
    # Get only Mastomys natalensis
    SPECIES == "MN", 
    # Only those sampled in MOSA field site
    GRID == "MOSA",
    # Remove 0g weights
    WEIGHT != 0,
    # Remove suspect IDs
    id.ok == TRUE,
    # Remove NAs:
    !is.na(DATE),   # dates
    !is.na(WEIGHT), # weights
    !is.na(TOE),    # ID
    !is.na(SEXCON)  # sexual condition
  ) %>% 
    mutate(
      # Define lactating (L) or not lactating (S; small nipples)
      LAC = case_when(
        str_detect(SEXCON, "L") ~ TRUE,
        str_detect(SEXCON, "S") ~ FALSE
      ),
      # Define pregnant (Y) or not pregnant (N)
      PRE = case_when(
        str_detect(SEXCON, "Y") ~ TRUE,
        str_detect(SEXCON, "N") ~ FALSE
      ),
      # # Define sexually active
      FER = case_when(
        # Males: sexually active if both testes are scrotal ("S") and gubernacula are visible ("V")
        SX == "M" & str_detect(SEXCON, "S") & str_detect(SEXCON, "V") ~ TRUE,
        SX == "M" ~ FALSE,
        # Females: sexually active if vagina is perforated ("P"), nipples lactating ("L"), or visibly pregnant ("Y")
        SX == "F" & (str_detect(SEXCON, "P") | str_detect(SEXCON, "L") | str_detect(SEXCON, "Y")) ~ TRUE,
        SX == "F" ~ FALSE
      )
    ) %>% 
  # Select relevant columns
  dplyr::select(DATE, WEIGHT, SEX = SX, TOE, FER, PRE, LAC)


# 3. Classify trapping sessions ----

# Identify separate study trapping sessions based on days since last trapping session

# Format data to create a trapping sessions group column
mosa <- mosa %>%
  # Arrange data by date
  arrange(DATE) %>%
  # Group sessions: if lagged difference > 3, assume different date session
  mutate(session = cumsum(c(1, diff(DATE) > 3))) %>%
  # Create session start date as the first date in each session group
  group_by(session) %>%
  mutate(session.start = first(DATE)) %>%
  # Count number of trapping days per session
  group_by(session.start) %>% 
  mutate(num.trap.days = length(unique(DATE))) %>% 
  ungroup() %>%
  # Drop session group (only need session start date as grouping factor)
  dplyr::select(!session)

# Define average no. days between trapping sessions
session.gap <- mosa %>% 
  dplyr::select(session.start) %>% 
  distinct() %>% 
  mutate(session.gap = c(0, diff(session.start))) %>% 
  summarise(mean(session.gap)) %>% 
  as.numeric() %>% 
  round()


# 4. Calculate population size ----

# Function for creating occurrence matrix: presence/absence matrix of individuals on each trapping session date
getOccMat <- function(cmr_data, date_col, choose_sex) {
  
  # cmr_data: with columns SEX, TOE (for ID), date named by date_col input
  # date_col: name of date column in cmr_data you want to search for presence/absence
  # choose_sex: name of SEX column value (e.g. "F" for female, "M" for male)
  
  # Filter data by input sex
  cmr_data_sex <- cmr_data %>% 
    filter(SEX == choose_sex)
  
  # Create presence/absence matrix of all individuals for each trapping session date
  occurrence.mat <- data.frame(matrix(nrow = cmr_data_sex %>% dplyr::select(all_of(date_col)) %>% unique() %>% nrow(),
                                      ncol = cmr_data_sex %>% dplyr::select(TOE) %>% unique() %>% nrow() + 1))
  colnames(occurrence.mat)[1] <- "trap.date"
  colnames(occurrence.mat)[2:ncol(occurrence.mat)] <- unique(cmr_data_sex$TOE)
  
  occurrence.mat$trap.date <- ymd(unique(cmr_data_sex[[date_col]]))

  # Fill matrix with presence/absence on each date
  ID.index  <- 2:ncol(occurrence.mat)
  individuals <- colnames(occurrence.mat)[ID.index]
  
  # Find which individuals are caught in this trapping session
  for (trap.date in occurrence.mat$trap.date) {
    
    # Get dataframe row index according to current trapping date
    row.index <- which(occurrence.mat$trap.date == trap.date)
    
    occurrence.mat[row.index, ID.index] <- individuals %in% cmr_data_sex$TOE[cmr_data_sex[[date_col]] == trap.date]
    
  }
  
  # Return occurrence matrix
  return(occurrence.mat)
  
}

# Get occurrence matrices for females & males
occ.mat.f <- getOccMat(cmr_data = mosa, date_col = "session.start", choose_sex = "F")
occ.mat.m <- getOccMat(cmr_data = mosa, date_col = "session.start", choose_sex = "M")

pop.f <- rowSums(occ.mat.f[2:ncol(occ.mat.f)])
plot(occ.mat.f$trap.date, pop.f, 'l')

pop.m <- rowSums(occ.mat.m[2:ncol(occ.mat.m)])
plot(occ.mat.m$trap.date, pop.m, 'l')

# Function to calculate minimum number alive (MNA) population size
calcMNA <- function(occ_mat, col_inds) {
  
  # occ_mat:  Occurrence matrix calculated by getOccMat()
  # col_inds: Column indices for individuals' IDs
  
  # individual considered present if found both before & after a given session
  for (inds in colnames(occ_mat)[col_inds]) {
    
    # Get column of presence/absence recordings for given individual
    ind.column <- occ_mat[, colnames(occ_mat) == inds]
    
    # If TRUE values surround FALSE value(s), set the FALSE to TRUE
    for (i in 1:(length(ind.column)-1)) {
      
      # If element is FALSE (individual not present)
      if (ind.column[i] == FALSE) {
        
        # Check if individual found before current session
        before_true <- any(ind.column[1:(i-1)] == TRUE)
        
        # Check if individual found after current session
        after_true <- any(ind.column[(i+1):length(ind.column)] == TRUE)
        
        # If both conditions are met, replace FALSE with TRUE
        if (before_true & after_true) {
          ind.column[i] <- TRUE
        }
        
      }
      
    }
    
    # Replace matrix column with updated presence/absence values
    occ_mat[, colnames(occ_mat) == inds] <- ind.column
    
  }
  
  # Calculate MNA by summing adjusted presences for each trapping session date
  occ_mat$mna <- rowSums(occ_mat[, col_inds])
  
  # Return dataframe of capture session date & MNA
  return(occ_mat %>% reframe(trap.date, mna))
  
}

# Calculate MNA for females & males 
MNA.f <- calcMNA(occ_mat = occ.mat.f, col_inds = 2:(ncol(occ.mat.f)-1))
MNA.m <- calcMNA(occ_mat = occ.mat.m, col_inds = 2:(ncol(occ.mat.m)-1))

# Join female & male MNAs into one dataframe
MNA <- full_join(
  MNA.f %>% rename(mna.female = mna),
  MNA.m %>% rename(mna.male = mna)) %>% 
  # 0 MNAs indicate no individuals of that sex present so set to 0
  mutate(mna.female = case_when(is.na(mna.female) ~ 0, TRUE ~ mna.female),
         mna.male   = case_when(is.na(mna.male)   ~ 0, TRUE ~ mna.male))

# Save MNA data for subsequent analyses
write.csv(MNA, "Data/MNA.by.sex.csv")

# Join calculated MNA (summed across both sexes) to MOSA data
mosa <- left_join(mosa, MNA %>% mutate(trap.date = ymd(trap.date)) %>% reframe(trap.date, mna = mna.female + mna.male),
                  by = c("session.start" = "trap.date")) %>%
  rename(POP = mna)


# 5. Classify sub-adults vs adults ----

# Only include females
mosa.f <- mosa %>% 
  filter(SEX == "F") %>% 
  rename(is.adult = FER) %>% 
  mutate(is.adult = as.numeric(is.adult))

# Save formatted CMR data (females only)
write.csv(mosa.f, "Data/MOSA_CMR_extracted_20240820.csv", row.names = FALSE)

# Save formatted CMR data (males & females)
write.csv(mosa, "Data/MOSA_CMR_bothsexes.csv", row.names = FALSE)


# 6. Calculate survival probability ----

# For each individual, find all captures & recaptures to infer survival:
# survival = 1 (TRUE) if recaptured in subsequent time step
# survival = 0 (FALSE) if not recaptured in subsequent time step

# Create array to store survival data in loop
mosa_colnames <- c("DATE", "session.start", "WEIGHT", "TOE", "FER", "POP", "is.adult")
survival <- as.data.frame(array(NA, c(0, length(mosa_colnames)+1)))
colnames(survival)[1:length(mosa_colnames)] <- mosa_colnames
colnames(survival)[length(mosa_colnames)+1] <- "SURV"

# Initialise iterator for within-trapping session observations of an individual
j <- 1

# Loop through all individuals
for (i in mosa.f$TOE %>% unique()){

  # Get data for individual
  subs <- mosa.f %>% filter(TOE == i) %>% mutate(DATE = as.Date(DATE))

  # Start with latest date individual captured (sorting dates in descending order)
  last_date <- TRUE
  sorted_dates <- sort(subs$DATE, decreasing = TRUE)

  # Loop through trap night dates, starting from the most recent
  for (d in seq_along(sorted_dates)) {

    # If this IS NOT not the last date the individual is captured
    if (!last_date){

      # If capture is in a different trapping session (>3d before subsequent capture)
      if (sorted_dates[d] < (sorted_dates[d-1] - 2)) {
        
        # Extract capture data for individual at this date
        date_subs <- subs[subs$DATE == sorted_dates[d], mosa_colnames]
        
        # If only 1 row in the capture data, store this
        if (nrow(date_subs) == 1) {
          survival[j, 1:length(mosa_colnames)] <- date_subs
          
          # Otherwise if >1 row in the capture data, store the 1st row
        } else {
          survival[j, 1:length(mosa_colnames)] <- date_subs[1, ]
        }
        
        # Set survival to TRUE (because individual is captured at a subsequent date)
        survival[j, "SURV"] <- TRUE
        
        # Increase row iterator by 1
        j <- j + 1
        
      }

      # If this IS the latest date the individual is captured
    } else {
      
      # Extract capture data for individual at this date
      date_subs <- subs[subs$DATE == sorted_dates[d], mosa_colnames]

      # If only 1 row in the capture data, store this
      if (nrow(date_subs) == 1) {
        survival[j, 1:length(mosa_colnames)] <- date_subs
        
        # Otherwise if >1 row in the capture data, store the 1st row
      } else {
        survival[j, 1:length(mosa_colnames)] <- date_subs[1, ]
      }

      # Set survival to FALSE (because individual is not recaptured after this date)
      survival[j, "SURV"] <- FALSE
      
      # Set last_date to FALSE (because looping to dates prior to last capture)
      last_date <- FALSE
      
      # Increase row iterator by 1
      j <- j + 1

    }

  }

}

# Format data before saving
survival.formatted <- survival %>% 
  # Convert date to date format
  mutate(DATE = as.Date(DATE),
         session.start = as.Date(session.start)) %>% 
  # Arrange individuals in chronological date order
  arrange(DATE, TOE)

# Save survival data
write.csv(survival.formatted, "Data/survival_data.csv", row.names = FALSE)


# 7. Calculate growth ----

# For each individual, identify all consecutive captures and record weights

# Initialise dataframe for storage
growth <- as.data.frame(array(NA, c(0, 8)))
colnames(growth) <- c("DATE1", "DATE2", "W1", "W2", "TOE", "FER1", "FER2", "POP")

# Prepare data for calculating growth
mosa_g <- mosa.f %>%
  filter(PRE == FALSE, # remove pregnant rodents
         LAC == FALSE) # remove lactating rodents

# Get unique ID of all individuals
inds_g <- unique(mosa_g$TOE)

# Initialise row iterator
j <- 1

# Loop through all individuals
for (i in inds_g) {

  # Get data for individual
  subs <- mosa_g %>% filter(TOE == i)

  # If individual has >1 capture record on different dates
  if (length(unique(subs$DATE)) > 1) {

    # Initialise date selector
    d1 <- 0

    # Loop through capture records for individual
    for (d in subs$DATE) {

      # If given date is >14 days after previous date
      if (d > (d1 + (session.gap / 2))) {

        # Update date selector
        d1 <- as.Date(d)

        # Are any capture records > 14d (0.5 trapping session) from selected date?
        condition1 <- subs$DATE > (d1 + (session.gap / 2))

        # Are any capture records < 42d (1.5 trapping sessions) from selected date?
        condition2 <- subs$DATE < (d1 + (session.gap * 1.5))

        # Find records where both conditions are TRUE
        if ( any(condition1 & condition2) ) {

          # If >1 record, get the date closest to selected date
          d2 <- min(subs$DATE[ subs$DATE > (d1 + (session.gap/2)) & subs$DATE < (d1 + (session.gap * 1.5)) ])

          # Store consecutive capture records for individual
          subs_date1 <- subs %>% filter(DATE == d1)  # First capture
          subs_date2 <- subs %>% filter(DATE == d2)  # Recapture

          # Store first capture record (if >1 on same date use first record)
          growth[j, c("DATE1","W1","FER1","POP")] <- subs_date1 %>%
            filter(row_number() == 1) %>%
            dplyr::select(DATE, WEIGHT, FER, POP)

          # Store recapture record (if >1 on same date use first record)
          growth[j, c("DATE2","W2","FER2")] <- subs_date2 %>%
            filter(row_number() == 1) %>%
            dplyr::select(DATE, WEIGHT, FER)

          # Record individual ID
          growth[j, "TOE"] <- i

          # Increase row iterator by 1
          j <- j + 1

        }

      }

    }

  }

}

# Format growth data
growth <- growth %>%
  mutate(
    DATE1 = as.Date(DATE1),
    DATE2 = as.Date(DATE2)
  )

# Save growth data
write.csv(growth, "Data/growth_data.csv", row.names = FALSE)


# 8. Calculate recruitment ----

# Define a function to process a single individual's data
process_individual <- function(subs, session.gap) {
  
  subs <- subs %>% arrange(DATE)  # data for individual
  results <- list()
  d1 <- as.Date("1900-01-01")  # Arbitrary old date
  
  for (d in subs$DATE) {
    # If given date is >14 days after previous date
    if (d > (d1 + days(session.gap / 2))) {
      
      d1 <- as.Date(d)
      
      # Define conditions
      condition1 <- subs$DATE >  (d1 + days(session.gap / 2))    # capture records > 14d (0.5 trapping session) from selected date
      condition2 <- subs$DATE <  (d1 + days(session.gap * 1.5))  # capture records < 42d (1.5 trapping sessions) from selected date
      condition3 <- subs$DATE >= (d1 + days(session.gap * 1.5))  # capture records >= 42d (1.5 trapping sessions) from selected date
      condition4 <- subs$DATE <  (d1 + days(session.gap * 2.5))  # capture records < 70d (2.5 trapping sessions) from selected date
      
      # 1 month recapture
      if (any(condition1 & condition2)) {
        
        # If >1 record, get the date closest to selected date
        d2 <- min(subs$DATE[condition1 & condition2])
        
        # Store consecutive capture records for individual
        rec1 <- subs %>% filter(DATE == d1) %>% slice(1)
        rec2 <- subs %>% filter(DATE == d2) %>% slice(1)
        
        # Store data for capture & recapture records
        PRE1 <- rec1$PRE
        LAC1 <- rec1$LAC
        PRE2 <- rec2$PRE
        LAC2 <- rec2$LAC
        TOE  <- rec1$TOE
        POP  <- rec1$POP
        WEIGHT <- rec1$WEIGHT
        
        # 2 month recapture
        if (any(condition3 & condition4)) {
          # If >1 record, get the date closest to selected date
          d3 <- min(subs$DATE[condition3 & condition4])
          LAC3 <- subs %>% filter(DATE == d3) %>% slice(1) %>% pull(LAC)
        } else {
          # If no recaptures at this time, set lactating status to FALSE
          LAC3 <- FALSE
        }
        
        # Store the record in next list element
        results[[length(results) + 1]] <- tibble(
          DATE = d1,
          WEIGHT = WEIGHT,
          TOE = TOE,
          POP = POP,
          PRE1 = PRE1,
          LAC1 = LAC1,
          PRE2 = PRE2,
          LAC2 = LAC2,
          LAC3 = LAC3
        )
      }
    }
  }
  
  bind_rows(results)
}

# Apply to all individuals
recr <- mosa.f %>% 
  group_split(TOE) %>% 
  map_dfr(~ process_individual(.x, session.gap))

# Final formatting and filtering
recr <- recr %>% 
  mutate(
    DATE = as.Date(DATE),
    LAC3 = if_else(is.na(LAC3), FALSE, LAC3),
    # Define recruitment: if individual is pregnant or lactating 1 month later, or lactating 2 months later
    REC = PRE2 | LAC2 | LAC3
  ) %>% 
  # Remove entries where individual started off pregnant or lactating
  filter(PRE1 != TRUE, LAC1 != TRUE) %>% 
  select(DATE, WEIGHT, TOE, POP, REC) %>% 
  filter(!is.na(WEIGHT), !is.na(POP), !is.na(REC))

# Double check looks OK 
if (doplot) {
  
  # Sample a year
  year.sample <- sample(min(year(recr$DATE)):max(year(recr$DATE)), size = 1)
  
  # Make plot of individuals in given year
  recr %>% 
    filter(year(DATE) == year.sample) %>%
    ggplot(aes(group = TOE)) + 
    geom_line(aes(DATE, REC)) + 
    labs(x = "Date", y = "Pregnant?") + 
    theme_minimal()
  
}

# Save recruitment data
write.csv(recr, "Data/recruitment_data.csv", row.names = FALSE)
