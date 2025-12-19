#' Create individual phenometrics dataset
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
#' @return A dataframe where each row contains the earliest and latest 
#' positive phenophase status observation ("yes") for an individual plant with  
#' a user-defined period of time.

create_individual_phenometrics <- function(year_start,
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
  
  # Create start/end dates for each period
  start_yrs <- year_start:year_end
  if (ymd(paste0(2025, "-", period_start)) > ymd(paste0(2025, "-", period_end))) {
    end_yrs <- start_yrs + 1
  } else {
    end_yrs <- start_yrs
  }
  period_starts <- paste0(start_yrs, "-", period_start)
  period_ends <- paste0(end_yrs, "-", period_end)
  
  # Get start/end date across all periods
  start_date <- first(period_starts)
  end_date <- last(period_ends)
  
  series <- create_series(
    start_date = start_date,
    end_date = end_date,
    request_source = request_source,
    max_yes_gap = max_yes_gap,
    site_ids = site_ids,
    species_ids = species_ids,
    phenophase_ids = phenophase_ids
  )
  
  # Append period number to series dataset
  periods <- data.frame(
    period_start = period_starts,
    period_end = period_ends
    ) %>%
    mutate(period_no = 1:length(start_yrs),
           julian_start = julian_date(period_start),
           julian_end = julian_date(period_end)) 
  n_periods <- nrow(periods)
  # If period is < 1 year, then rows for days between periods
  if (length(period_starts) > 1 & 
      as.numeric(as.Date(period_starts[2]) - as.Date(period_ends[1])) > 1) {
    period_starts <- as.Date(period_starts)
    period_ends <- as.Date(period_ends)
    periods_add <- data.frame(
      period_start = period_ends[-length(period_ends)] + 1,
      period_end = period_starts[-1] - 1,
      period_no = seq(1.5, max(periods$period_no) - 0.5, by = 1)
      ) %>%
      mutate(julian_start = julian_date(period_start),
             julian_end = julian_date(period_end)) %>%
      mutate(period_start = as.character(period_start),
             period_end = as.character(period_end))
    periods <- rbind(periods, periods_add) %>%
      arrange(period_start)
  }
  
  # Add unique series ID
  series <- series %>%
    mutate(series_id = paste(id, php, first_yes_date, sep = "_"))
  
  # Create columns that indicate which period the first yes and last yes
  # dates fall into
  series2 <- series %>%
    cross_join(select(periods, period_start, period_end, period_no)) %>%
    mutate(first_yes_period = ifelse(first_yes_date >= period_start &
                                       first_yes_date <= period_end, 
                                     period_no, 0),
           last_yes_period = ifelse(last_yes_date >= period_start &
                                      last_yes_date <= period_end, 
                                    period_no, 0)) %>%
    group_by(series_id) %>%
    summarize(first_yes_period = max(first_yes_period),
              last_yes_period = max(last_yes_period)) %>%
    mutate(last_yes_period = ifelse(last_yes_period == 0, 
                                    max(periods$period_no) + 1, 
                                    last_yes_period)) %>%
    data.frame()
  # check:
  # count(series2, first_yes_period, last_yes_period)
  series <- series %>%
    left_join(series2, by = "series_id")
  
  # Get dataframe with just periods (not in-between times if period < 1 year)
  periodso <- periods %>%
    filter(period_no %% 1 == 0)
  
  # Create columns in series dataframe to indicate whether any dates (inclusive)
  # between first and last yes fall into each period
  period_matrix <- matrix(NA, 
                          nrow = nrow(series), 
                          ncol = nrow(periodso))
  for (i in 1:nrow(series)) {
    for (j in 1:ncol(period_matrix)) {
      period_matrix[i, j] <- ifelse(
        length(intersect(series$first_yes_julian[i]:series$last_yes_julian[i],
                         periodso$julian_start[j]:periodso$julian_end[j])) == 0,
        0, 1)
    }
  }
  period_df <- as.data.frame(period_matrix)
  colnames(period_df) <- paste0("period", periodso$period_no)
  period_df$nperiods <- rowSums(period_df)
  # check:
  # cols <- colnames(period_df)
  # period_df %>% count(!!!rlang::syms(cols))
  series <- cbind(series, period_df)
  
  # Remove any series that only occurred between periods (if periods < 1 year)
  series <- series %>%
    filter(nperiods > 0)

  # Identify how many series do NOT start and begin in the same period, and warn 
  # user if the proportion of series is higher than a given threshold
  noverlap <- sum(!series$first_yes_period %in% 1:n_periods | 
                    !series$last_yes_period %in% 1:n_periods |
                    series$first_yes_period != series$last_yes_period)
  propoverlap <- noverlap/nrow(series)
  if (propoverlap > 0.10) {
    warning("More than 10% of yes series overlap period boundaries")
  }
  
  # Extract plant info to add back in later on
  plantinfo <- series %>%
    distinct(site, species_id, common_name, id)
  
  if (include_overlap_series) {
    # If we want to include information from series that span mulitple periods:
      
    # To identify the earliest first yes in each period for each individual, we 
    # need to gather series that had a first yes in that period AND series that 
    # had yeses in that period as well as yeses in any prior period. For the 
    # latter type of series that overlap period boundaries, we'll change the
    # first yes date to the first date of the period and we'll create a flag.
    
    # IMPORTANT NOTE FOR FIRST YESES: if max_prior_no is set to a numeric value
    # and not NA, then any series that overlaps period boundaries will not be 
    # included in an individual phenometrics dataset since days since prior no 
    # will be set to NA for series that cross period boundaries.
      
    # To identify the latest last yes in each period for each individual, we 
    # need to gather series that had a last yes in that period AND series that 
    # had yeses in that period as well as yeses in any later period. For the 
    # latter type of series that overlap period boundaries, we'll change the 
    # last yes date to the last date of the period and we'll create a flag.
    
    # First yeses -------------#
    
    # Extract the simple ones, with first yes in periods of interest
    ip_first <- series %>%
      filter(first_yes_period %in% 1:n_periods) %>%
      select(series_id, id, php, first_yes_date, prior_no_date, n_yeses, 
             # flag_multobs, flag_conflict, 
             prior_no, first_yes_julian, prior_no_julian, first_yes_period) %>%
      mutate(flag_beganprior = 0)
    
    # Find series that overlap period boundaries
    series_overlap <- series %>%
      filter(first_yes_period != last_yes_period) 
    
    for (i in 1:nrow(series_overlap)) {
      
      # Identify first_yes_period(s) that we need to add:
      new_periods <- ceiling(series_overlap$first_yes_period[i]):floor(series_overlap$last_yes_period[i])
      new_periods <- new_periods[new_periods != series_overlap$first_yes_period[i]]
      if (length(new_periods) == 0) {next}

      row_add <- data.frame(
        series_id = series_overlap$series_id[i],
        id = series_overlap$id[i], 
        php = series_overlap$php[i],
        prior_no_date = NA,
        prior_no_julian = NA,
        prior_no = NA,
        n_yeses = NA,
        # flag_multobs = series_overlap$flag_multobs[i],
        # flag_conflict = series_overlap$flag_conflict[i],
        first_yes_period = new_periods
      )
      row_add <- filter(row_add, first_yes_period <= max(periods$period_no))
      if (nrow(row_add) == 0) {next}
      row_add <- row_add %>%  
        mutate(first_yes_date = periods$period_start[periods$period_no == first_yes_period],
               first_yes_julian = julian_date(first_yes_date),
               flag_beganprior = 1)
      if (!exists("rows_add_first")) {
        rows_add_first <- row_add
      } else {
        rows_add_first <- rbind(rows_add_first, row_add)
      }
    }
    rows_add_first <- rows_add_first %>%
      select(colnames(ip_first))
    
    # Merge these series together
    ip_first <- rbind(ip_first, rows_add_first)

    # Filter by max_prior_no (if not NA)
    if (!is.na(max_prior_no)) {
      ip_first <- ip_first %>%
        filter(prior_no <= max_prior_no)
    }
    
    # Then find the earliest first yes for each individual, phenophase and
    # period
    ip_first <- ip_first %>%
      group_by(id, php, first_yes_period) %>%
      summarize(first_yes_date = first_yes_date[first_yes_julian == min(first_yes_julian)],
                prior_no_date = prior_no_date[first_yes_julian == min(first_yes_julian)],
                prior_no = prior_no[first_yes_julian == min(first_yes_julian)],
                n_first_series = n(),
                # flag_multobs = ifelse(sum(flag_multobs) > 0, 1, 0),
                # flag_conflict = ifelse(sum(flag_conflict) > 0, 1, 0),
                flag_beganprior = ifelse(sum(flag_beganprior) > 0, 1, 0),
                .groups = "keep") %>% 
      data.frame()
    
    # Last yeses -------------#
    
    # Extract the simple ones, with last yes in periods of interest
    ip_last <- series %>%
      filter(last_yes_period %in% 1:n_periods) %>%
      select(series_id, id, php, last_yes_date, next_no_date, n_yeses, 
             # flag_multobs, flag_conflict, 
             next_no, last_yes_julian, next_no_julian, last_yes_period) %>%
      mutate(flag_endedlater = 0)
    
    # Find series that overlap period boundaries
    series_overlap <- series %>%
      filter(first_yes_period != last_yes_period)

    for (i in 1:nrow(series_overlap)) {
      
      # Identify last_yes_period(s) that we need to add:
      new_periods <- ceiling(series_overlap$first_yes_period[i]):floor(series_overlap$last_yes_period[i])
      new_periods <- new_periods[new_periods != series_overlap$last_yes_period[i]]
      if (length(new_periods) == 0) {next}
      
      row_add <- data.frame(
        series_id = series_overlap$series_id[i],
        id = series_overlap$id[i], 
        php = series_overlap$php[i],
        next_no_date = NA,
        next_no_julian = NA,
        next_no = NA,
        n_yeses = NA,
        # flag_multobs = series_overlap$flag_multobs[i],
        # flag_conflict = series_overlap$flag_conflict[i],
        last_yes_period = new_periods
      )
      row_add <- filter(row_add, last_yes_period > 0)
      if (nrow(row_add) == 0) {next}
      row_add <- row_add %>%  
        mutate(last_yes_date = periods$period_end[periods$period_no == last_yes_period],
               last_yes_julian = julian_date(last_yes_date),
               flag_endedlater = 1)
      if (!exists("rows_add_last")) {
        rows_add_last <- row_add
      } else {
        rows_add_last <- rbind(rows_add_last, row_add)
      }
    }
    rows_add_last <- rows_add_last %>%
      select(colnames(ip_last))
    
    # Merge these series together
    ip_last <- rbind(ip_last, rows_add_last)
    
    # Find the latest last yes for each individual, phenophase and period
    ip_last <- ip_last %>%
      group_by(id, php, last_yes_period) %>%
      summarize(last_yes_date = last_yes_date[last_yes_julian == max(last_yes_julian)],
                next_no_date = next_no_date[last_yes_julian == max(last_yes_julian)],
                next_no = next_no[last_yes_julian == max(last_yes_julian)],
                n_last_series = n(),
                # flag_multobs = ifelse(sum(flag_multobs) > 0, 1, 0),
                # flag_conflict = ifelse(sum(flag_conflict) > 0, 1, 0),
                flag_endedlater = ifelse(sum(flag_endedlater) > 0, 1, 0),
                .groups = "keep") %>% 
      data.frame()
    
  } else {  # If ignoring series that span multiple periods
    
    # First yeses -------------#
    
    # Extract the simple ones, with first yes in periods of interest
    ip_first <- series %>%
      filter(first_yes_period %in% 1:n_periods) %>%
      select(series_id, id, php, first_yes_date, prior_no_date, n_yeses, 
             # flag_multobs, flag_conflict, 
             prior_no, first_yes_julian, prior_no_julian, first_yes_period) %>%
      mutate(flag_beganprior = 0)

    # Filter by max_prior_no (if not NA)
    if (!is.na(max_prior_no)) {
      ip_first <- ip_first %>%
        filter(prior_no <= max_prior_no)
    }
    
    # Then find the earliest first yes for each individual, phenophase and
    # period
    ip_first <- ip_first %>%
      group_by(id, php, first_yes_period) %>%
      summarize(first_yes_date = first_yes_date[first_yes_julian == min(first_yes_julian)],
                prior_no_date = prior_no_date[first_yes_julian == min(first_yes_julian)],
                prior_no = prior_no[first_yes_julian == min(first_yes_julian)],
                n_first_series = n(),
                # flag_multobs = ifelse(sum(flag_multobs) > 0, 1, 0),
                # flag_conflict = ifelse(sum(flag_conflict) > 0, 1, 0),
                flag_beganprior = ifelse(sum(flag_beganprior) > 0, 1, 0),
                .groups = "keep") %>% 
      data.frame()
    
    # Last yeses -------------#
    
    # Extract the simple ones, with last yes in periods of interest
    ip_last <- series %>%
      filter(last_yes_period %in% 1:n_periods) %>%
      select(series_id, site, species_id, common_name, id, php, last_yes_date,
             next_no_date, n_yeses, 
             # flag_multobs, flag_conflict, 
             next_no, last_yes_julian, next_no_julian, last_yes_period) %>%
      mutate(flag_endedlater = 0)
    
    # Find the latest last yes for each individual, phenophase and period
    ip_last <- ip_last %>%
      group_by(id, php, last_yes_period) %>%
      summarize(last_yes_date = last_yes_date[last_yes_julian == max(last_yes_julian)],
                next_no_date = next_no_date[last_yes_julian == max(last_yes_julian)],
                next_no = next_no[last_yes_julian == max(last_yes_julian)],
                n_last_series = n(),
                # flag_multobs = ifelse(sum(flag_multobs) > 0, 1, 0),
                # flag_conflict = ifelse(sum(flag_conflict) > 0, 1, 0),
                flag_endedlater = ifelse(sum(flag_endedlater) > 0, 1, 0),
                .groups = "keep") %>% 
      data.frame()

  }

  # Merge first/last information for each individual, phenophase, period
  ip_first <- ip_first %>%
    rename(period_no = first_yes_period)
  ip_last <- ip_last %>%
    rename(period_no = last_yes_period)

  ip <- ip_first %>%
    full_join(ip_last, by = c("id", "php", "period_no")) %>%
    left_join(select(periods, period_no, period_start, period_end), by = "period_no") %>%
    left_join(plantinfo, by = "id") %>%
    mutate(first_yes_julian = julian_date(first_yes_date),
           prior_no_julian = julian_date(prior_no_date),
           last_yes_julian = julian_date(last_yes_date),
           next_no_julian = julian_date(next_no_date)) %>%
    select(site, species_id, common_name, id, php, period_no, period_start,
           period_end, n_first_series, first_yes_date, first_yes_julian, 
           prior_no, prior_no_date, prior_no_julian, flag_beganprior,
           n_last_series, last_yes_date, last_yes_julian, 
           next_no, next_no_date, next_no_julian, flag_endedlater) %>%
    arrange(period_no, php, common_name, id)
  
  return(ip)
}  

  
  
  