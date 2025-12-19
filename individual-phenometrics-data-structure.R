library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)

# Using the following data to look at how individual phenometrics data are 
# constructed at the beginning of each period (by default, period = calendar year)
  # McDowell (network_id = 622)
  # Jojoba (species_id = 769)
  # Flowers (500) or open flowers (501)

# Status-intensity data
si <- npn_download_status_data(
  request_source = "erinz",
  years = 2019:2025,
  network_ids = 622,
  species_ids = 769,
  phenophase_ids = c(500, 501)
) %>% data.frame()

# Individual phenometrics, calendar year
ipc <- npn_download_individual_phenometrics(
  request_source = "erinz",
  years = 2020:2025,
  network_ids = 622, 
  species_ids = 769,
  phenophase_ids = c(500, 501)
) %>% data.frame()

# Individual phenometrics, water year
ipw <- npn_download_individual_phenometrics(
  request_source = "erinz",
  years = 2019:2024,
  network_ids = 622, 
  species_ids = 769,
  phenophase_ids = c(500, 501),
  period_start = "10-01",
  period_end = "09-30"
) %>% data.frame()

# Notes on data structure when using water year: 

# Downloads data for water years, beginning with Oct in first year listed.
# So here, 6 water years, last one incomplete (Oct2019-Sep2020, ..., Oct 2024-Sep2025)

# First yes year, doy are related to the CALENDAR year, not the water year
# If we want to select the first yes in a water year, need to create that 
# variable first.

# Histogram of open flower dates
filter(si, phenophase_id == 501 & phenophase_status == 1) %>% 
  pull(day_of_year) %>%
  hist(breaks = 50)

# Histogram with flower dates
filter(si, phenophase_id == 500 & phenophase_status == 1) %>% 
  pull(day_of_year) %>%
  hist(breaks = 50)

# How much data for each plant?
count(si, individual_id)

# Look at some examples of status changes near the turn of the year to see
# how it's reported in the individual phenometrics datasets

# No at last obs of 2022, yes at first obs of 2023
si %>%
  filter(individual_id == 117051) %>%
  filter(phenophase_id == 500) %>%
  filter(observation_date > "2022-12-01" & observation_date < "2023-01-30") %>%
  select(phenophase_id, phenophase_status, observation_date, day_of_year, intensity_value)
  # No on 2022-12-31; Yes on 2023-01-02 

  # Individual phenometrics calculated for calendar year:
  ipc %>%
    filter(individual_id == 117051) %>%
    filter(phenophase_id == 500) %>%
    filter(first_yes_year == 2023) %>%
    select(phenophase_id, first_yes_year, first_yes_doy, numdays_since_prior_no)
    # first_yes_doy == 2; numdays_since_prior_no == NA
  
  # Individual phenometrics calculated for water year:
  ipw %>%
    filter(individual_id == 117051) %>%
    filter(phenophase_id == 500) %>%
    filter(first_yes_year == 2023) %>%
    select(phenophase_id, first_yes_year, first_yes_doy, numdays_since_prior_no)
    # first_yes_doy == 2; numdays_since_prior_no == 2

# Yes at last obs of 2021, yes at first obs of 2022
si %>%
  filter(individual_id == 117061) %>%
  filter(phenophase_id == 500) %>%
  filter(observation_date > "2021-12-15" & observation_date < "2022-01-30") %>%
  arrange(phenophase_id, observation_date) %>%
  select(phenophase_id, phenophase_status, observation_date, day_of_year, intensity_value)
  # No on 2021-12-27; Yes on 2021-12-31; Yes on 2022-01-06, 2022-01-13; No on 2022-01-18

  # Individual phenometrics calculated for calendar year:
  ipc %>%
    filter(individual_id == 117061) %>%
    filter(phenophase_id == 500) %>%
    filter(first_yes_year == 2022) %>%
    select(phenophase_id, first_yes_year, first_yes_doy, numdays_since_prior_no)
    # first_yes_doy == 6; numdays_since_prior_no == NA
  
  # Individual phenometrics calculated for water year:
  ipw %>%
    filter(individual_id == 117061) %>%
    filter(phenophase_id == 500) %>%
    filter(first_yes_year == 2022) %>%
    select(phenophase_id, first_yes_year, first_yes_doy, numdays_since_prior_no)
  # first_yes_doy == 19; numdays_since_prior_no == 1 (Jan 6th is in the middle
  # of a yes series that overlaps the new year, so it isn't a first yes here)

