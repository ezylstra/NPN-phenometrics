#' Calculate standard error (SE), but return NA if all values = NA
#' 
#' @param vector vector of values we'd like to summarize oover
#'
#' @details 
#' 
#' @return numeric SE or NA

se_na <- function(x) {
  
  if (sum(is.na(x)) == length(x)) {
    se <- NA
  } else {
    x <- x[!is.na(x)]
    se <- sd(x) / sqrt(length(x))
  }
  
  return(se)
}
  
  