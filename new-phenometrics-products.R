# Using functions to create new phenometrics datasets

library(dplyr)
library(lubridate)
library(rnpn)

# Load functions
functions <- list.files("functions", full.names = TRUE)
for (f in functions) {
  source(f)
}

# Series ----------------------------------------------------------------------#

# User inputs
start_date <-  "2022-01-01"
end_date <- "2023-06-30"
max_yes_gap <- 90
request_source <- "erinz"
site_ids <- c(24702, 24703, 24705)
species_ids = c(16, 210, 317, 769, 945, 1022, 1170, 2133)
phenophase_ids <- c(500, 501)

# series <- create_series(
#   start_date = start_date,
#   end_date = end_date,
#   request_source = request_source,
#   max_yes_gap = max_yes_gap,
#   site_ids = site_ids,
#   species_ids = species_ids,
#   phenophase_ids = phenophase_ids
# )

# Individual phenometrics -----------------------------------------------------#

# User inputs
year_start <- 2017
year_end <- 2024
period_start <- "10-01" # "01-01"
period_end <- "09-30"   # "12-31"
onset_offset <- "both" # can be "onset", "offset", or "both"
max_prior_no <- 30
include_overlap_series <- TRUE
request_source <- "erinz"
max_yes_gap <- 90
site_ids <- c(24702, 24703, 24705)
# species_ids <- c(16, 210, 317, 769, 945, 1022, 1170, 2133)
species_ids <- 769
phenophase_ids <- c(500, 501)

ip <- create_individual_phenometrics(
  year_start = year_start,
  year_end = year_end,
  period_start = period_start,
  period_end = period_end,
  onset_offset = onset_offset,
  max_prior_no = max_prior_no,
  include_overlap_series = include_overlap_series,
  request_source = request_source,
  max_yes_gap = max_yes_gap,
  site_ids = site_ids,
  species_ids = species_ids,
  phenophase_ids = phenophase_ids
)

ip %>%
  filter(metric == "onset") %>%
  group_by(common_name, phenophase_id, period, period_start, period_end) %>%
  summarize(n_plants = n(),
            first_yes = mean(first_yes_date),
            prior_no = mean(days_prior_no, na.rm = TRUE),
            # last_yes = mean(last_yes_date),
            .groups = "keep") %>%
  data.frame()




# Site phenometrics -----------------------------------------------------------#

# User inputs
year_start <- 2021
year_end <- 2023
period_start <- "10-01" # "01-01"
period_end <- "09-30"   # "12-31"
max_prior_no <- NA
include_overlap_series <- TRUE
request_source <- "erinz"
max_yes_gap <- 90
site_ids <- c(24702, 24703, 24705)
species_ids <- c(16, 210, 317, 769, 945, 1022, 1170, 2133)
phenophase_ids <- c(500, 501)

sp <- create_site_phenometrics(
  year_start = year_start,
  year_end = year_end,
  period_start = period_start,
  period_end = period_end,
  max_prior_no = max_prior_no,
  include_overlap_series = include_overlap_series,
  request_source = request_source,
  max_yes_gap = max_yes_gap,
  site_ids = site_ids,
  phenophase_ids = phenophase_ids
)

