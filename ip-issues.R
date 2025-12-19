# User inputs
year_start <- 2017
year_end <- 2024
period_start <- "12-01" # "01-01"
period_end <- "11-30"   # "12-31"
max_prior_no <- 30
include_overlap_series <- FALSE
request_source <- "erinz"
max_yes_gap <- 90
species_ids <- 769
phenophase_ids <- c(500, 501)

# Run create_individual_phenometrics through line 214...

# Find overlap series
filter(series, first_yes_period != last_yes_period) %>%
  arrange(id, php ,first_yes_period) %>%
  select(id, php, first_yes_date, prior_no, last_yes_date, next_no,
         first_yes_period, last_yes_period)

# id = 117072; php == 500, 1-2 and 2-3 overlap

# Current workflow:
filter(ip, id == 117072, php == 500) %>%
  select(id, php, period_no, period_start, period_end, first_yes_date,
         prior_no, last_yes_date, next_no)
# Don't see any problems, but that's because there are lots of other series in 
# each period.

# If we remove all other series for period 2
series <- series %>%
  filter(!(id == 117072 & php == 500 & first_yes_period == 2 & last_yes_period == 2))
# Current workflow with include_overlap_series = FALSE:
filter(ip, id == 117072, php == 500) %>%
  select(id, php, period_no, period_start, period_end, first_yes_date,
         prior_no, last_yes_date, next_no)
# Problem evident now. For period 2, first_yes date is AFTER last_yes_date ######

# Current workflow with include_overlap_series = TRUE:
filter(ip, id == 117072, php == 500) %>%
  select(id, php, period_no, period_start, period_end, first_yes_date,
         prior_no, last_yes_date, next_no)
# No problem if including overlap series





filter(series, id == 117072, php == 500, first_yes_period %in% 1:2) %>%
  select(id, php, first_yes_date, prior_no, last_yes_date, next_no, 
         first_yes_period, last_yes_period, series_id)


