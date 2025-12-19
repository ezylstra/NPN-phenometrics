#' Create site phenometrics dataset
#' 
#' @param year_start numeric year indicating the year that the first "period"
#' will begin
#' @param year_end numeric year indicating the year that the last "period"
#' will begin
#' @param period_start month-day, formatted as "mm-dd" indicating the month and
#' day that each annual period will begin (default = "01-01", indicating Jan 1st
#' for calendar year. Set to "10-01" for water year or "07-01" for summer year).
#' It is typical to specify that each period is 1 year in duration.
#' @param period_end month-day, formatted as "mm-dd" indicating the month and
#' day that each annual period will end (default = "12-31", indicating Dec 31st
#' for calendar year. Set to "09-30" for water year or "06-30" for summer year).
#' It is typical to specify that each period is 1 year in duration.
#' @param max_prior_no maximum number of days a prior no occurs before the first
#' yes (default = 30). If NA, then any first yes will be included, regardless of
#' if/when a prior no occurred.
#' @param include_overlap_series logical indicating whether to use series that 
#' overlap one or more period boundaries when identifiying the first and last 
#' yes for an individual in each period (default = FALSE). If TRUE, the first 
#' yes date will be set to the first date of the period and thus will 
#' automatically be the earliest first yes across all series for that individual 
#' and period, so use with caution. Note that setting max_prior_no to any 
#' numeric value will override this choice, as series that overlap period 
#' boundaries will have prior no set to NA.
#' @param request_source character string indicating who is requesting the data
#' (required)
#' @param max_yes_gap maximum number of days between consecutive yeses before
#' the 2nd yes is considered the start of a new series (default = 90).
#' @param site_ids 
#' @param species_ids
#' @param phenophase_ids
#'
#' @details 
#' 
#' @return A dataframe where each row contains the mean earliest and mean latest
#' date of a positive phenophase status observation ("yes") within a 
#' user-defined period of time for all plants of a given species at a given site

create_site_phenometrics <- function(year_start,
                                     year_end,
                                     period_start = "01-01",
                                     period_end = "12-31",
                                     max_prior_no = 30,
                                     include_overlap_series = FALSE,
                                     request_source,
                                     max_yes_gap = 90,
                                     site_ids = NULL, 
                                     species_ids = NULL, 
                                     phenophase_ids = NULL) {

  # Extract the name of this function for reporting
  function_name <- as.character(match.call())[1]
  
  # Libraries required for this function to work
  dependencies <- c("dplyr", "lubridate", "rnpn")
  if (!all(unlist(lapply(X = dependencies, FUN = require, character.only = TRUE)))) {
    stop("At least one package required by ", function_name, 
         " could not be loaded: ", paste(dependencies, collapse = ", "),
         " are required.")
  }
  
  # Make sure dates are appropriate
  if (!year_start %in% 2000:2050) {
    stop(function_name, " requires a valid start year on or after 2000")
  }
  if (!year_end %in% 2000:2050) {
    stop(function_name, " requires a valid end year on or after 2000")
  }
  pstart_test <- suppressWarnings(parse_date_time(paste0(2024, "-", period_start),
                                                  orders ="%Y-%m-%d"))
  if (is.na(pstart_test)) {
    stop(function_name, " requires a valid period start formatted as 'MM-DD'")
  }
  pend_test <- suppressWarnings(parse_date_time(paste0(2024, "-", period_end),
                                                  orders ="%Y-%m-%d"))
  if (is.na(pend_test)) {
    stop(function_name, " requires a valid period end formatted as 'MM-DD'")
  }

  # Make sure request source is specified
  if (is.na(request_source)) {
    stop(function_name, " requires name of request source")
  }
  
  # Make sure max_yes_gap is numeric
  if (!is.numeric(max_yes_gap)) {
    stop(function_name, " requires numeric value for max_yes_gap")
  }
  
  # Make sure max_prior_no is numeric (or NA)
  if (!(is.na(max_prior_no) | is.numeric(max_prior_no))) {
    stop(function_name, " requires numeric value or NA for max_prior_no")
  }

  ip <- create_individual_phenometrics(
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
  
  # We've already done all the filtering by max_yes_gap, max_prior_no, and 
  # series overlaps with the previous functions. Now it's just a matter of 
  # calulating means and standard errors...

  sp <- ip %>% 
    group_by(site, species_id, common_name, php, period_no, period_start, 
             period_end) %>%
    summarize(first_yes_n = sum(!is.na(first_yes_date)),
              mean_first_yes = mean_na(first_yes_date),
              se_first_yes = se_na(first_yes_julian),
              mean_prior_no = mean_na(prior_no),
              se_prior_no = se_na(prior_no),
              last_yes_n = sum(!is.na(last_yes_date)),
              mean_last_yes = mean_na(last_yes_date),
              se_last_yes = se_na(last_yes_julian),
              mean_next_no = mean_na(next_no),
              se_next_no = se_na(next_no),
              # Flag if at least one series overlapped the start of the period
              flag_beganprior = ifelse(sum(flag_beganprior) > 0, 1, 0),
              # Flag if at least one series overlapped the end of the period
              flag_endedlater = ifelse(sum(flag_endedlater) > 0, 1, 0),
              .groups = "keep") %>%
    mutate(mean_first_julian = julian_date(mean_first_yes),
           mean_last_julian = julian_date(mean_last_yes)) %>%
    arrange(php, species_id, period_no, site) %>%
    data.frame()
  
  return(sp)
}  

  
  
  