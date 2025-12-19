# Compare site and individual phenometrics data structure
library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)

# Download data for red maple trees in TN and NC in 2024

# Individual phenometrics
ip_orig <- npn_download_individual_phenometrics(
  request_source = "erinz",
  years = 2024,
  states = c("TN", "NC"),
  species_ids = 3,
  phenophase_ids = c(371, 483, 498, 500, 501, 516)
) %>% data.frame()

# Site phenometrics
sp_orig <- npn_download_site_phenometrics(
  request_source = "erinz",
  years = 2024,
  states = c("TN", "NC"),
  species_ids = 3,
  phenophase_ids = c(371, 483, 498, 500, 501, 516)
) %>% data.frame()

head(ip_orig)
head(sp_orig)

# Simplify dataframes
sp <- sp_orig %>%
  select(-c(mean_first_yes_julian_date, mean_last_yes_julian_date,
            latitude, longitude, elevation_in_meters, state, species_id, genus,
            species, kingdom)) %>%
  rename(first_yes_n = first_yes_sample_size,
         first_yes_mn = mean_first_yes_doy,
         first_yes_se = se_first_yes_in_days,
         prior_no_mn = mean_numdays_since_prior_no,
         prior_no_se = se_numdays_since_prior_no, 
         last_yes_n = last_yes_sample_size,
         last_yes_mn = mean_last_yes_doy,
         last_yes_se = se_last_yes_in_days,
         next_no_mn = mean_numdays_until_next_no,
         next_no_se = se_numdays_until_next_no)

ip <- ip_orig %>%
  select(-c(first_yes_julian_date, last_yes_julian_date, first_yes_month, 
            first_yes_day, last_yes_month, last_yes_day,
            latitude, longitude, elevation_in_meters, state, species_id, genus,
            species, kingdom)) %>%
  rename(first_yes = first_yes_doy,
         prior_no = numdays_since_prior_no,
         last_yes = last_yes_doy,
         next_no = numdays_until_next_no)

# Number of trees at each site
sites <- ip %>%
  group_by(site_id) %>%
  summarize(ntrees = n_distinct(individual_id),
            nphp = n_distinct(phenophase_description),
            nfirsts = n_distinct(paste0(individual_id, "_", phenophase_id))) %>%
  data.frame()

# Look at one site, one phenophase
site <- 2772
php <- 500 # 516, 501, 483
filter(sp, site_id == site, phenophase_id == php)
filter(ip, site_id == site, phenophase_id == php) %>% 
  arrange(phenophase_id, individual_id)

# Looks like means/SEs in site phenoometrics are only calculated if the number
# of days since prior no (or number of days until next no for last yes) is
# <= 30 (or at least not NA)

# What if we had multiple first yeses for a plant and phenophase? Are both
# included in the means/SEs?

ip %>%
  group_by(individual_id, phenophase_id) %>%
  summarize(n_obs = n(), .groups = "keep") %>%
  data.frame() %>%
  filter(n_obs > 1)
# id = 23058, phenophase = 371

filter(ip, individual_id == 23058, phenophase_id == 371) # site_id == 8182

# First yes mean DOY: just used first observation in each year for each plant
# Last yes mean DOY: just used last observation in each year for each plant

filter(ip, site_id == 8182, phenophase_id == 371)
filter(sp, site_id == 8182, phenophase_id == 371)


# Which observation would they use if there were first yeses around turn of the
# year if using a water or summer year?
# Look at different datasets to figure this out... Jojoba

# Individual phenometrics, water year
ip_orig2 <- npn_download_individual_phenometrics(
  request_source = "erinz",
  years = 2020:2023,
  network_ids = 622, 
  species_ids = 769,
  phenophase_ids = c(500, 501),
  period_start = "10-01",
  period_end = "09-30"
) %>% data.frame()

# Site phenometrics, water year
sp_orig2 <- npn_download_site_phenometrics(
  request_source = "erinz",
  years = 2020:2023,
  network_ids = 622, 
  species_ids = 769,
  phenophase_ids = c(500, 501),
  period_start = "10-01",
  period_end = "09-30"
) %>% data.frame()

ip2 <- ip_orig2 %>%
  select(-c(first_yes_month, 
            first_yes_day, last_yes_month, last_yes_day,
            latitude, longitude, elevation_in_meters, state, species_id, genus,
            species, kingdom)) %>%
  rename(first_yes = first_yes_doy,
         prior_no = numdays_since_prior_no,
         last_yes = last_yes_doy,
         next_no = numdays_until_next_no)

sp2 <- sp_orig2 %>%
  select(-c(latitude, longitude, elevation_in_meters, state, species_id, genus,
            species, kingdom)) %>%
  rename(first_yes_n = first_yes_sample_size,
         first_yes_mn = mean_first_yes_doy,
         first_yes_se = se_first_yes_in_days,
         prior_no_mn = mean_numdays_since_prior_no,
         prior_no_se = se_numdays_since_prior_no, 
         last_yes_n = last_yes_sample_size,
         last_yes_mn = mean_last_yes_doy,
         last_yes_se = se_last_yes_in_days,
         next_no_mn = mean_numdays_until_next_no,
         next_no_se = se_numdays_until_next_no)

# Look at one site, one phenophase
ip2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500) %>%
  filter(first_yes_year == 2020 | (first_yes_year == 2021 & first_yes < 70)) %>%
  arrange(individual_id, first_yes_year, first_yes)

sp2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500)
# Note that there are only 4 rows, with mean_first_yes_year = 2021:2024 but 
# mean_last_yes_year = 2021, 2021, 2023, 2024??

# Julian dates are in there, which would allow us to calculate the mean first 
# date within each WATER year a phenophase occurred.
# HOWEVER, it look like within each water year, they've just averaged the 
# minimum day of calendar year a phenophase was observed for each individual:

# Calculating this incorrectly based on individual phenometrics
ip2_temp <- ip2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500) %>%
  group_by(first_yes_year, individual_id) %>%
  summarize(first_yes_doy = min(first_yes[!is.na(prior_no)]),
            last_yes_doy = max(last_yes))
ip2_temp %>%
  group_by(first_yes_year) %>%
  summarize(nobs = n(),
            first_yes_mn = round(mean(first_yes_doy)),
            first_yes_se = sd(first_yes_doy)/sqrt(nobs),
            last_yes_mn = round(mean(last_yes_doy)),
            last_yes_se = sd(last_yes_doy)/sqrt(nobs)) %>%
  data.frame()
# first_yes_mn and first_yes_se for 2021-2024 match up with what's in the 
# site phenometrics dataset.


count(sp2, mean_first_yes_year, mean_last_yes_year)
# Not sure how calculations are being done for last yeses... 
# Stuff at the very bottom is trying to figure that out....

# How often are people downloading site phenometrics data? And of those, do
# people every specify a period other than the default calendar year?

# -----------------------------------------------------------------------------#
# Made same downloads of jojoba data from POP just to make sure it's providing
# same data as rnpn

ip_pop_orig <- read.csv(file.choose())
sp_pop_orig <- read.csv(file.choose())

ip_pop <- ip_pop_orig %>%
  mutate(NumDays_Since_Prior_No = ifelse(NumDays_Since_Prior_No == -9999, NA, NumDays_Since_Prior_No)) %>%
  mutate(NumDays_Until_Next_No = ifelse(NumDays_Until_Next_No == -9999, NA, NumDays_Until_Next_No))
sp_pop <- sp_pop_orig

colnames(ip_pop) <- str_to_lower(colnames(ip_pop))
colnames(sp_pop) <- str_to_lower(colnames(sp_pop))

ip_pop2 <- ip_pop %>%
  select(-c(first_yes_month, 
            first_yes_day, last_yes_month, last_yes_day,
            latitude, longitude, elevation_in_meters, state, species_id, genus,
            species, kingdom)) %>%
  rename(first_yes = first_yes_doy,
         prior_no = numdays_since_prior_no,
         last_yes = last_yes_doy,
         next_no = numdays_until_next_no)

sp_pop2 <- sp_pop %>%
  select(-c(latitude, longitude, elevation_in_meters, state, species_id, genus,
            species, kingdom)) %>%
  rename(first_yes_n = first_yes_sample_size,
         first_yes_mn = mean_first_yes_doy,
         first_yes_se = se_first_yes_in_days,
         prior_no_mn = mean_numdays_since_prior_no,
         prior_no_se = se_numdays_since_prior_no, 
         last_yes_n = last_yes_sample_size,
         last_yes_mn = mean_last_yes_doy,
         last_yes_se = se_last_yes_in_days,
         next_no_mn = mean_numdays_until_next_no,
         next_no_se = se_numdays_until_next_no)

all.equal(ip_pop2[, -1], ip2)
filter(ip2, is.na(next_no))
ip_pop2[which(is.na(ip2$next_no)), -1]

# Looks like in the dataset from the POP, it provided the numdays_until_next_no
# for a couple rows when the yes ws on Sep 28 and the no was on Oct 1.
# In contrast, the dataset from rnpn had NA for numdays_until_next_no.

# Otherwise they're identical...

# -----------------------------------------------------------------------------#

ip2 <- ip2 %>%
  mutate(first_yes_date = parse_date_time(x = paste(first_yes_year, first_yes), orders = "yj"),
         last_yes_date = parse_date_time(x = paste(last_yes_year, last_yes), orders = "yj"),
         first_wateryr = ifelse(month(first_yes_date) %in% 10:12, first_yes_year + 1, first_yes_year),
         last_wateryr = ifelse(month(last_yes_date) %in% 10:12, last_yes_year + 1, last_yes_year))

count(ip2, first_yes_year, first_wateryr, last_yes_year, last_wateryr)
filter(ip2, phenophase_id == 500) %>%
  count(first_yes_year, first_wateryr, last_yes_year, last_wateryr)
filter(ip2, first_yes_year == 2021, first_wateryr == 2022, last_yes_year == 2022, last_wateryr == 2022)

ip2_tempF <- ip2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500) %>%
  group_by(first_wateryr, individual_id) %>%
  summarize(first_yes_doy = min(first_yes[!is.na(prior_no)]))

ip2_tempL <- ip2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500) %>%
  group_by(last_wateryr, individual_id) %>%
  summarize(last_yes_doy = max(last_yes))

# The first dates match up (but are incorrect -- using minimum day of calendar year)
ip2_tempF %>%
  group_by(first_wateryr) %>%
  summarize(nobs = n(),
            first_yes_mn = round(mean(first_yes_doy)),
            first_yes_se = sd(first_yes_doy)/sqrt(nobs)) %>%
  data.frame()
sp2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500)

ip2_temp <- ip2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500) %>%
  group_by(first_wateryr, individual_id) %>%
  summarize(first_yes_doy = min(first_yes[!is.na(prior_no)]),
            last_yes_doy = max(last_yes),
            last_yes_doy_noNA = max(last_yes[!is.na(next_no)]),
            .groups = "keep") %>%
  data.frame()

ip2_temp %>%
  group_by(first_wateryr) %>%
  summarize(nobs = n(),
            first_yes_mn = round(mean(first_yes_doy)),
            first_yes_se = sd(first_yes_doy)/sqrt(nobs),
            last_yes_mn = round(mean(last_yes_doy)),
            last_yes_se = sd(last_yes_doy)/sqrt(nobs),
            last_yes_mn2 = round(mean(last_yes_doy_noNA)),
            last_yes_se2 = sd(last_yes_doy_noNA)/sqrt(nobs)) %>%
  data.frame()
  data.frame()

ip2_temp <- ip2 %>%
  filter(site_id == 24702) %>%
  filter(phenophase_id == 500) %>%
  group_by(first_wateryr, last_wateryr, individual_id) %>%
  summarize(first_yes_doy = min(first_yes[!is.na(prior_no)]),
            last_yes_doy = max(last_yes),
            
            .groups = "keep") %>%
  data.frame()


