#' Calculate mean value, but return NA if all values = NA
#' 
#' @param vector vector of values we'd like to average
#'
#' @details 
#' 
#' @return numeric mean or NA

mean_na <- function(x) {
  
  if (sum(is.na(x)) == length(x)) {
    mn <- NA
  } else {
    mn <- mean(x, na.rm = TRUE)
  }
  
  return(mn)
}
  
  