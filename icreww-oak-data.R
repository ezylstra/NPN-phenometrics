# Creating series data for oak species (to use for ICREWW)

library(dplyr)
library(lubridate)
library(rnpn)
library(ggplot2)
library(tidyr)

# Load functions
functions <- list.files("functions", full.names = TRUE)
for (f in functions) {
  source(f)
}

# Identify species, phenophases -----------------------------------------------#

# First, get species IDs
spp <- npn_species() %>%
  select(-species_type) %>%
  data.frame()
oaks <- spp %>%
  filter(genus == "Quercus") %>%
  filter(species %in% c("alba",
                        "bicolor",
                        "geminata",
                        "lyrata",
                        "macrocarpa",
                        "michauxii",
                        "montana",
                        "muehlenbergii",
                        "prinoides",
                        "stellata",
                        "virginiana"))

# Then, get phenophases
phps <- npn_phenophases_by_species(
  species_ids = oaks$species_id,
  date = "2026-01-01"
) %>% data.frame()
oak_phps <- phps %>% 
  filter(pheno_class_id %in% c(1, 7)) %>%
  arrange(pheno_class_id) %>%
  count(pheno_class_id, phenophase_id, phenophase_name)
oak_phps <- oak_phps %>%
  filter(phenophase_name != "Full flowering")

# Create series data ----------------------------------------------------------#

# User inputs
start_date <- "2009-01-01"
end_date <- "2026-12-31"
request_source <- "erinz-icreww"
max_yes_gap <- 90
species_ids <- oaks$species_id
site_ids <- NULL
phenophase_ids <- oak_phps$phenophase_id

oak_series <- create_series(
  start_date = start_date,
  end_date = end_date,
  request_source = request_source,
  max_yes_gap = max_yes_gap,
  species_ids = species_ids,
  phenophase_ids = phenophase_ids
)

# Look at a few summaries -----------------------------------------------------#

# Series counts for each species
count(oak_series, common_name)

# Look at series start dates, by year (2020-2026) & phenophase
ggplot(filter(oak_series, first_yes_year >= 2020 & days_prior_no <= 30)) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 1) +
  facet_grid(first_yes_year ~ phenophase_description)

# Look at BLB series starts, by species
ggplot(filter(oak_series, days_prior_no <= 30 & phenophase_description == "Breaking leaf buds")) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 1) +
  facet_wrap(~common_name)

# Look at flowers series starts, by species
ggplot(filter(oak_series, days_prior_no <= 30 & phenophase_description == "Open flowers")) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 1) +
  facet_wrap(~common_name)

# Look at number of days since prior no, by year (2020-2026) & phenophase
ggplot(filter(oak_series, days_prior_no <= 100 & first_yes_year >= 2020)) +
  geom_histogram(aes(x = days_prior_no), binwidth = 1) +
  geom_vline(xintercept = c(7, 14, 30), color = "blue") + 
  facet_grid(first_yes_year ~ phenophase_description)

# Create onsets ---------------------------------------------------------------#
# First yes of the calendar year. Restricting to those with prior no within 30 days

# User inputs (in addition to those used above for series dataset)
year_start <- 2009
year_end <- 2026
period_start <- "01-01"
period_end <- "12-31"
onset_offset <- "onset"
include_overlap_series <- FALSE
max_prior_no <- 30

oak_onset <- create_individual_phenometrics(
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

# Onset counts for each species
count(oak_onset, common_name)

# Look at first yes dates, by year (2020-2026) & phenophase
ggplot(filter(oak_onset, year(period_start) >= 2020 & days_prior_no <= 30)) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 1) +
  facet_grid(year(period_start) ~ phenophase_description)

# Look at BLB first yeses, by species
ggplot(filter(oak_onset, days_prior_no <= 30 & phenophase_description == "Breaking leaf buds")) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 1) +
  facet_wrap(~common_name)

# Look at flowers first yeses, by species
ggplot(filter(oak_onset, days_prior_no <= 30 & phenophase_description == "Open flowers")) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 1) +
  facet_wrap(~common_name)

# Create offsets --------------------------------------------------------------#
# Last yes of the calendar year

# User inputs (in addition to those used above for series dataset)
onset_offset <- "offset"

oak_offset <- create_individual_phenometrics(
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

# Offset counts for each species
count(oak_offset, common_name)

# Look at last yes dates, by year (2020-2026) & phenophase
ggplot(filter(oak_offset, year(period_start) >= 2020)) +
  geom_histogram(aes(x = last_yes_doy), binwidth = 1) +
  facet_grid(year(period_start) ~ phenophase_description)

# Look at BLB first yeses, by species
ggplot(filter(oak_offset, phenophase_description == "Breaking leaf buds")) +
  geom_histogram(aes(x = last_yes_doy), binwidth = 1) +
  facet_wrap(~common_name)

# Look at flowers first yeses, by species
ggplot(filter(oak_offset, phenophase_description == "Open flowers")) +
  geom_histogram(aes(x = last_yes_doy), binwidth = 1) +
  facet_wrap(~common_name)

# Combine onset and offset dates for each individual, phenophase, year --------#
# Only include those individuals with both a first and last yes that year

oak_combined <- oak_onset %>%
  select(-c(metric, began_prior)) %>%
  inner_join(select(oak_offset, individual_id, phenophase_id, period, 
                   last_yes_date, last_yes_doy, next_no_date, days_next_no),
            by = c("individual_id", "phenophase_id", "period")) %>%
  data.frame()

# Calculate time elapsed between first, last yes
oak_combined <- oak_combined %>%
  mutate(duration = last_yes_doy - first_yes_doy)

# Visualize phenophase duration, by year (2020-2026) & phenophase
ggplot(filter(oak_combined, year(period_start) >= 2020)) +
  geom_histogram(aes(x = duration), binwidth = 1) +
  facet_grid(year(period_start) ~ phenophase_description)
# Lots of same day (duration = 0)
ggplot(filter(oak_combined, year(period_start) >= 2020 &
                duration %in% 1:150)) +
  geom_histogram(aes(x = duration), binwidth = 1) +
  facet_grid(year(period_start) ~ phenophase_description)

# Visualize BLB duration, by species (with duration > 0, < 150)
ggplot(filter(oak_combined, phenophase_description == "Breaking leaf buds" &
                duration %in% 1:150)) +
  geom_histogram(aes(x = duration), binwidth = 1) +
  facet_wrap(~common_name)

# Visualize open flowering duration, by species (with duration > 0, < 150)
ggplot(filter(oak_combined, phenophase_description == "Open flowers" &
                duration %in% 1:150)) +
  geom_histogram(aes(x = duration), binwidth = 1) +
  facet_wrap(~common_name)

# Write datasets to file ------------------------------------------------------#

write.csv(oak_series, "output/icreww-oaks/oak-series.csv", row.names = FALSE)
write.csv(oak_onset, "output/icreww-oaks/oak-onset.csv", row.names = FALSE)
write.csv(oak_offset, "output/icreww-oaks/oak-offset.csv", row.names = FALSE)
write.csv(oak_combined, "output/icreww-oaks/oak-combined.csv", row.names = FALSE)
