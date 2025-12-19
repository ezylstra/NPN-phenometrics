#' Create series dataset
#' 
#' @param start_date date, provided as "YYYY-MM-DD", that indicates the earliest
#' date an output series begins (default = "2009-01-01")
#' @param end_date date, provided as "YYYY-MM-DD", that indicates the latest
#' date an output series begins (default = today's date)
#' @param requestor character string indicating who is requesting the data
#' (required)
#' @param max_yes_gap maximum number of days between consecutive yeses before
#' the 2nd yes is considered the start of a new series (default = 90).
#' @param site_ids 
#' @param species_ids
#' @param phenophase_ids
#'
#' @details 
#' 
#' @return A dataframe, where each row provides information about a series of
#' positive phenophase status observations ("yeses")

create_series <- function(start_date = "2009-01-01",
                          end_date = Sys.Date(),
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
  start_test <- suppressWarnings(parse_date_time(start_date, orders ="%Y-%m-%d"))
  if (is.na(start_test)) {
    stop(function_name, " requires a valid start date formatted as 'YYYY-MM-DD'")
  }
  end_test <- suppressWarnings(parse_date_time(end_date, orders ="%Y-%m-%d"))
  if (is.na(end_test)) {
    stop(function_name, " requires a valid end date formatted as 'YYYY-MM-DD'")
  }
  if (start_date >= end_date) {
    stop(function_name, " requires end date to be later than start date")
  }
  
  # Make sure request source is specified
  if (is.na(request_source)) {
    stop(function_name, " requires name of request source")
  }
  
  # Make sure max_yes_gap is numeric
  if (!is.numeric(max_yes_gap)) {
    stop(function_name, " requires numeric value for max_yes_gap")
  }
  
  # Extract data from period extending >= 1 year before and after desired periood 
  # to get proper calculation of prior and following nos ##### KEY STEP ######
  start_date <- ymd(start_date)
  start_yr <- year(start_date)
  start_extract <- start_yr - 1
  end_date <- ymd(end_date)
  end_yr <- year(end_date)
  end_extract <- end_yr + 1
  
  # Download status-intensity data
  si_orig <- npn_download_status_data(
    request_source = request_source,
    years = start_extract:end_extract,
    station_ids = site_ids,
    species_ids = species_ids,
    phenophase_ids = phenophase_ids
  ) %>% data.frame()

  # Remove unnecssary columns for now and rename to make things easier to view
  si <- si_orig %>%
    select(-c(update_datetime, latitude, longitude, elevation_in_meters,
              state, genus, species, kingdom, phenophase_description, 
              abundance_value)) %>%
    rename(obsid = observation_id,
           site = site_id, 
           id = individual_id, 
           php = phenophase_id, 
           obsdate = observation_date,
           doy = day_of_year,
           status = phenophase_status,
           intensity_cat = intensity_category_id,
           intensity = intensity_value)

  # Remove intensity data
  ser <- si %>%
    select(-c(intensity_cat, intensity))
  
  # Arrange by individual, date
  ser <- ser %>%
    arrange(common_name, species_id, site, id, obsdate, doy, php, status)
  
  # Discard observations with unknown status (-1)
  ser <- ser %>%
    filter(status >= 0)
  
  # For now, will try to create multiple observer flag and conflict flag
  # If both yes and no reported on same day, assume yes but flag 
  # that this occurred
  # Flag if more than one observer on same day, regardless of status
  ser <- ser %>%
    group_by(site, species_id, common_name, id, php, obsdate, doy) %>%
    summarize(status = max(status),   
              flag_conflict = ifelse(n_distinct(status) > 1, 1, 0),
              flag_multobs = ifelse(n() > 1, 1, 0), 
              .groups = "keep") %>%
    data.frame()
  
  # Add some fake data to test max gap thing...
  # fake <- data.frame(site = 24702,
  #                    species_id = 210,
  #                    common_name = "saguaro",
  #                    id = 1111111,
  #                    php = 500,
  #                    obsdate = ymd(c("2021-02-01", "2021-02-02", "2021-02-05",
  #                                    "2021-02-06", "2021-08-01", "2021-08-03",
  #                                    "2021-08-05", "2021-10-01")),
  #                    doy = yday(c("2021-02-01", "2021-02-02", "2021-02-05",
  #                                 "2021-02-06", "2021-08-01", "2021-08-03",
  #                                 "2021-08-05", "2021-10-01")),
  #                    status = c(0, rep(1, 5), 0, 1),
  #                    flag_multobs = 0,
  #                    flag_conflict = 0)
  # ser <- rbind(fake, ser)
  
  # Identify unique combinations of individual plant and phenophase
  combos <- ser %>%
    distinct(site, species_id, common_name, id, php)
  
  # Loop through each individual-phenophase combination
  for (i in 1:nrow(combos)) {
    ser1 <- ser %>%
      filter(id == combos$id[i] & php == combos$php[i])
    
    rles <- rle(ser1$status)
    
    # If there are no series of yeses, skip to next combination
    nseries <- sum(rles$value == 1)
    if (nseries == 0) {next}
    
    # Get row number (index) for the start of each 0 or 1 series
    starts <- cumsum(c(1, rles$lengths))
    starts <- starts[-length(starts)]
    
    # Extract table with information about each run of 0s and 1s
    series01 <- data.frame(value = rles$values,
                           length = rles$lengths,
                           startrow = starts) %>%
      mutate(endrow = length + startrow - 1)
    # Extract table with information about each run of 1s (series)
    yesseries <- series01 %>% filter(value == 1) %>%
      rename(n_yeses = length) %>%
      select(-value) %>%
      mutate(lastnorow = ifelse(startrow == 1, NA, startrow - 1)) %>%
      mutate(nextnorow = ifelse(endrow == nrow(ser1), NA, endrow + 1))
    
    # Identify any yesseries that have problems with max_yes_gap
    ys <- yesseries
    for (j in 1:nrow(ys)) {
      
      # If there's only one yes in a series, then there is no gap (and no
      # problem). But if there's more than one yes, extract the gaps between
      # consecutive yeses to see if they exceed the user-defined limit
      if (ys$n_yeses[j] == 1) {
        gaps <- 0
      } else {
        gaps <- as.numeric(ser1$obsdate[(ys$startrow[j] + 1):ys$endrow[j]] -
                             ser1$obsdate[ys$startrow[j]:(ys$endrow[j] - 1)]) 
      }
      # Identify number of new series that need to be created
      newseries <- which(gaps > max_yes_gap)
      
      # If there were consecutive yeses separated by more than the user-defined
      # maximum (max_yes_gap), then create a new series starting at the 2nd yes
      if (length(newseries) == 0) {
        yesseries_new <- ys[j,]
      } else {
        yesseries_new <- data.frame(
          startrow = c(ys$startrow[j], ys$startrow[j] + newseries),
          endrow = c(ys$startrow[j] + newseries - 1, ys$endrow[j]),
          lastnorow = c(ys$lastnorow[j], NA),
          nextnorow = c(NA, ys$nextnorow[j])
        ) %>%
          mutate(n_yeses = endrow - startrow + 1, .before = startrow)
      }
      
      if (j == 1) {
        yesseries <- yesseries_new
      } else {
        yesseries <- rbind(yesseries, yesseries_new)
      }
    }
    
    yesseries$first_yes_date <- ser1$obsdate[yesseries$startrow]
    yesseries$last_yes_date <- ser1$obsdate[yesseries$endrow] 
    
    # If there are no prior nos for any yes series, make all prior no dates NA.
    # Otherwise, add in date of prior no
    if (sum(is.na(yesseries$lastnorow)) == nrow(yesseries)) {
      yesseries$prior_no_date <- NA
    } else {
      yesseries$prior_no_date <- ser1$obsdate[yesseries$lastnorow]
    }
    # If there are no next nos for any yes series, make all next no dates NA.
    # Otherwise, add in date of next No
    if (sum(is.na(yesseries$nextnorow)) == nrow(yesseries)) {
      yesseries$next_no_date <- NA
    } else {
      yesseries$next_no_date <- ser1$obsdate[yesseries$nextnorow]
    }
    
    # Append information about the series
    yesseries$site <- combos$site[i]
    yesseries$species_id <- combos$species_id[i]
    yesseries$common_name <- combos$common_name[i]
    yesseries$id <- combos$id[i]
    yesseries$php <- combos$php[i]
    yesseries <- yesseries %>%
      select(site, species_id, common_name, id, php, first_yes_date, 
             prior_no_date, n_yeses, last_yes_date, next_no_date)
    
    if (i == 1) {
      series <- yesseries
    } else {
      series <- rbind(series, yesseries)
    }
  }
  
  # Only proceed if there are one or more series....
  if (!exists("series")) {
    
    stop("No series for the selected species and phenophases")
    
  } else {
    
    # Only keep series that have first or last yes within user-defined 
    # start/end dates
    series <- series %>%
      filter((first_yes_date >= start_date & first_yes_date <= end_date) |
               (last_yes_date >= start_date & last_yes_date <= end_date))
    
    if (nrow(series) == 0) {
      
      stop("No series for the selected species and phenophases have a first ",
           "or last yes within user-selected start and end dates")
      
    } else {
      
      # Append flags to the series dataset....
      # If one or more yeses in a series had multiple observers, then flag
      # If one or more yeses in a series had a conflict, then flag
      flag1 <- ser %>%
        filter(status == 1)
      
      series <- series %>%
        mutate(flag_multobs = NA,
               flag_conflict = NA)
      
      for (i in 1:nrow(series)) {
        # Find observation dates in flag1 that fall in series
        flag1sub <- flag1 %>%
          filter(id == series$id[i] & php == series$php[i]) %>%
          filter(obsdate >= series$first_yes_date[i] &
                   obsdate <= series$last_yes_date[i])
        series$flag_multobs[i] <- sum(flag1sub$flag_multobs)
        series$flag_conflict[i] <- sum(flag1sub$flag_conflict)
      }
      series <- series %>%
        mutate(flag_multobs = ifelse(flag_multobs > 0, 1, 0),
               flag_conflict = ifelse(flag_conflict > 0, 1, 0))
      
      # Calculate days since prior no, until next no
      series <- series %>%
        mutate(prior_no = as.numeric(first_yes_date - prior_no_date),
               next_no = as.numeric(next_no_date - last_yes_date))
      
      # Add in Julian dates
      series <- series %>%
        mutate(first_yes_julian = julian_date(first_yes_date),
               prior_no_julian = julian_date(prior_no_date),
               last_yes_julian = julian_date(last_yes_date),
               next_no_julian = julian_date(next_no_date))
      
    }
  }
  
  return(series)
}
