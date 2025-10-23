# Smart Data Setup for Survival Analysis
# Functions to prepare survival data for Bayesian imputation

#' Prepare Survival Data for Bayesian Imputation
#'
#' This function prepares survival data for use with bayesian_impute().
#' It can auto-detect time and status columns, validate the data format and create
#' a properly formatted survival dataset object.
#'
#' @param data A data.frame containing survival data
#' @param time Column name or index for survival times. If NULL, attempts auto-detection
#'   using common survival analysis column names.
#' @param status Column name or index for event status (1=event, 0=censored). 
#'   If NULL, attempts auto-detection.
#' @param auto_detect Logical; whether to attempt automatic column detection (default: TRUE)
#' @param validate Logical; whether to run data validation (default: TRUE)
#' @param verbose Logical; print information about detected columns (default: TRUE)
#' @param time_unit Character string describing the unit used for time values (default "days")
#'
#' @return A survival_data object (enhanced data.frame) ready for bayesian_impute()
#'
#' @details
#' The function attempts to automatically detect survival data columns using common
#' naming conventions:
#' 
#' **Time columns**: "time", "duration", "followup", "survival_time", "days", 
#' "months", "years", "event_time", "tte"
#' 
#' **Status columns**: "status", "event", "died", "death", "censored", "outcome", 
#' "vital_status", "event_indicator"
#' 
#' The resulting object retains all original columns but adds metadata about
#' which columns represent time and status for use by bayesian_impute().
#'
#' @examples
#' \dontrun{
#' # Auto-detection (recommended for beginners)
#' library(survival)
#' lung_data <- prepare_survival_data(lung)
#' 
#' # Manual specification  
#' my_data <- prepare_survival_data(mydata, time = "duration", status = "died")
#' 
#' # Then use directly with bayesian_impute
#' result <- bayesian_impute(lung_data, n_imputations = 20)
#' }
#'
#' @seealso \code{\link{bayesian_impute}}, \code{\link{validate_survival_data}}
#' @export
prepare_survival_data <- function(data, 
                                 time = NULL, 
                                 status = NULL,
                                 auto_detect = TRUE,
                                 validate = TRUE,
                                 verbose = TRUE,
                                 time_unit = "days") {
  
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame")
  }
  
  if (nrow(data) == 0) {
    stop("Data cannot be empty")
  }
  
  # Auto-detect time column if not specified
  if (is.null(time) && auto_detect) {
    time_candidates <- c("time", "duration", "followup", "survival_time", 
                        "days", "months", "years", "event_time", "tte",
                        "futime", "surv_time", "survival")
    
    time <- detect_survival_column(data, time_candidates, "time")
    
    if (verbose && !is.null(time)) {
      message("Auto-detected time column: '", time, "'")
    }
  }
  
  # Auto-detect status column if not specified  
  if (is.null(status) && auto_detect) {
    status_candidates <- c("status", "event", "died", "death", "censored", 
                          "outcome", "vital_status", "event_indicator",
                          "cens", "delta", "observed")
    
    status <- detect_survival_column(data, status_candidates, "status")
    
    if (verbose && !is.null(status)) {
      message("Auto-detected status column: '", status, "'")
    }
  }
  
  # Check that we have both time and status
  if (is.null(time) || is.null(status)) {
    stop("Could not detect survival data columns. Please specify 'time' and 'status' manually.\n",
         "Available columns: ", paste(names(data), collapse = ", "))
  }
  
  # Convert column names to character if they were provided as indices
  if (is.numeric(time)) {
    time <- names(data)[time]
  }
  if (is.numeric(status)) {
    status <- names(data)[status]
  }
  
  # Validate columns exist
  if (!time %in% names(data)) {
    stop("Time column '", time, "' not found in data")
  }
  if (!status %in% names(data)) {
    stop("Status column '", status, "' not found in data")
  }
  
  # Run validation if requested
  if (validate) {
    validation_result <- validate_survival_data(data, time, status)
    data <- validation_result$data  # Use the potentially converted data
  }
  
  # Create enhanced data object
  # Normalise time_unit to a simple label
  if (is.null(time_unit) || !nzchar(time_unit)) time_unit <- "days"
  time_unit <- tolower(as.character(time_unit))

  result <- structure(
    data,
    time_col = time,
    status_col = status,
    n_total = nrow(data),
    n_censored = sum(data[[status]] == 0),
    censoring_rate = mean(data[[status]] == 0),
    time_unit = time_unit,
    class = c("survival_data", class(data))
  )
  
  if (verbose) {
    message("Successfully prepared survival data with ", nrow(data), " observations")
    message("Censoring rate: ", round(attr(result, "censoring_rate") * 100, 1), "%")
    message("Time unit: ", attr(result, "time_unit"))
  }
  
  return(result)
}

#' Detect Survival Analysis Columns
#'
#' Helper function to automatically detect time or status columns in survival data
#' based on common naming conventions.
#'
#' @param data Data.frame to search
#' @param candidates Character vector of candidate column names to search for
#' @param col_type Type of column being detected ("time" or "status") for error messages
#'
#' @return Character name of detected column, or NULL if none found
#' @keywords internal
detect_survival_column <- function(data, candidates, col_type) {
  
  # Check for exact matches first (case-insensitive)
  col_names_lower <- tolower(names(data))
  candidates_lower <- tolower(candidates)
  
  for (candidate in candidates_lower) {
    matches <- which(col_names_lower == candidate)
    if (length(matches) > 0) {
      return(names(data)[matches[1]])  # Return first match
    }
  }
  
  # Check for partial matches (contains the candidate)
  for (candidate in candidates_lower) {
    matches <- which(grepl(candidate, col_names_lower, fixed = TRUE))
    if (length(matches) > 0) {
      return(names(data)[matches[1]])  # Return first match
    }
  }
  
  # If this is a status column, also check for numeric patterns with smart prioritisation
  if (col_type == "status") {
    # Get all binary columns (0/1 and 1/2 coding)
    binary_01_cols <- sapply(data, function(x) is.numeric(x) && all(x %in% c(0, 1, NA)))
    binary_12_cols <- sapply(data, function(x) is.numeric(x) && all(x %in% c(1, 2, NA)))
    
    binary_cols <- names(data)[binary_01_cols | binary_12_cols]
    
    if (length(binary_cols) > 0) {
      # Score each binary column based on how likely it is to be a status variable
      status_scores <- sapply(binary_cols, function(col) {
        score <- 0
        col_lower <- tolower(col)
        
        # High priority: clear status/event indicators
        if (grepl("status|event|outcome|death|died|censored", col_lower)) score <- score + 10
        if (grepl("response|failure|survival|recurrence", col_lower)) score <- score + 8
        if (grepl("progression|relapse|remission", col_lower)) score <- score + 8
        
        # Medium priority: could be status
        if (grepl("indicator|flag|binary|dichotomous", col_lower)) score <- score + 5
        if (grepl("yes|no|present|absent", col_lower)) score <- score + 3
        
        # Low priority: likely not status variables
        if (grepl("sex|gender|male|female", col_lower)) score <- score - 15
        if (grepl("group|treatment|arm|cohort|study", col_lower)) score <- score - 12
        if (grepl("smoking|alcohol|drug|medication", col_lower)) score <- score - 10
        if (grepl("age|bmi|weight|height", col_lower)) score <- score - 8
        if (grepl("race|ethnicity|education", col_lower)) score <- score - 8
        if (grepl("stage|grade|score|index", col_lower)) score <- score - 5
        
        # Very low priority: definitely not status
        if (grepl("id|patient|subject|sample", col_lower)) score <- score - 20
        if (grepl("date|time|duration|follow", col_lower)) score <- score - 20
        
        return(score)
      })
      
      # Return the highest scoring column
      best_col <- binary_cols[which.max(status_scores)]
      best_score <- max(status_scores)
      
      # Only return if the score is reasonable (not too negative)
      if (best_score >= -5) {
        return(best_col)
      } else {
        # If all scores are too low, return the first binary column but warn
        warning("Multiple binary columns found but none clearly identified as status. ",
                "Using '", binary_cols[1], "' as status column. ",
                "Consider specifying status column manually.")
        return(binary_cols[1])
      }
    }
  }
  
  return(NULL)
}


#' Print method for survival_data objects, primarily an older debugging/validation tool
#' @param x A survival_data object
#' @param ... Additional arguments (ignored)
#' @export
print.survival_data <- function(x, ...) {
  cat("Survival Data\n")
  cat("=============\n")
  cat("Observations:", attr(x, "n_total"), "\n")
  cat("Time column: '", attr(x, "time_col"), "'\n", sep = "")
  cat("Status column: '", attr(x, "status_col"), "'\n", sep = "")
  cat("Censored:", attr(x, "n_censored"), 
      "(", round(attr(x, "censoring_rate") * 100, 1), "%)\n", sep = "")
  cat("\nFirst few rows:\n")
  print(head(as.data.frame(x)))
  
  invisible(x)
}

#' Quick Survival Analysis Setup
#'
#' One-liner function for users who want to go from raw data to imputation
#' with sensible defaults. Useful for beginners or quick analyses.
#'
#' @param data Raw survival data
#' @param groups Name of the group variable for group comparison (optional)
#' @param n_imputations Number of imputed datasets (default: 10)
#' @param distribution Survival distribution to use (default: "weibull")
#' @param time_unit Label describing the time unit to display (default "days")
#' @param ... Additional arguments passed to `bayesian_impute()`
#'
#' @return bayesian_imputation object
#' @export
quick_impute <- function(data, groups = NULL, n_imputations = 10, distribution = "weibull", time_unit = "days", ...) {
  
  # Delegate to unified wrapper
  result <- impute(
    data = data,
    groups = groups,
    n_imputations = n_imputations,
    distribution = distribution,
    time_unit = time_unit,
    ...
  )
  
  return(result)
}