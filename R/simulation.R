#' Simulation Framework for Bayesian Survival Imputation
#'
#' This file contains functions for generating simulated survival data,
#' applying controlled censoring mechanisms.

#' Null-Coalescing Operator
#'
#' Return the left-hand side unless it is `NULL`, otherwise return the
#' right-hand side. Used internally for defaulting optional parameters.
#'
#' @name grapes-or-or-grapes
#' @aliases %||%
#' @usage x \%||\% y
#' @param x Primary value.
#' @param y Fallback value when `x` is `NULL`.
#'
#' @return `x` if non-`NULL`, otherwise `y`.
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Generate Survival Data from Known Distribution
#'
#' Generates survival times from a specified parametric distribution with known parameters.
#' This allows validation against ground truth.
#'
#' @param n Number of observations to generate
#' @param distribution Distribution family ("weibull", "exponential", "lognormal")
#' @param params Named list of distribution parameters
#' @param seed Random seed for reproducibility
#'
#' @return Vector of survival times
#'
#' @details
#' Distribution parameters:
#' - Weibull: list(shape = α, scale = λ)
#' - Exponential: list(rate = λ)
#' - Log-normal: list(meanlog = μ, sdlog = σ)
#'
#' @examples
#' \dontrun{
#' # Generate Weibull data 
#' times <- simulate_survival_data(100, "weibull", list(shape = 2, scale = 4))
#' 
#' # Generate exponential data
#' times <- simulate_survival_data(100, "exponential", list(rate = 0.5))
#' }
#' @export
simulate_survival_data <- function(n, distribution, params, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  checkmate::assert_int(n, lower = 1)
  checkmate::assert_choice(distribution, c("weibull", "exponential", "lognormal"))
  checkmate::assert_list(params, min.len = 1)
  
  # Generate survival times based on distribution
  times <- switch(distribution,
    "weibull" = {
      checkmate::assert_names(names(params), must.include = c("shape", "scale"))
      rweibull(n, shape = params$shape, scale = params$scale)
    },
    "exponential" = {
      checkmate::assert_names(names(params), must.include = "rate")
      rexp(n, rate = params$rate)
    },
    "lognormal" = {
      checkmate::assert_names(names(params), must.include = c("meanlog", "sdlog"))
      rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog)
    }
  )
  
  # Ensure all times are positive
  if (any(times <= 0)) {
    warning("Some generated times were <= 0, replacing with small positive values")
    times[times <= 0] <- 1e-6
  }
  
  return(times)
}

#' Apply Halabi-Singh Censoring Method
#'
#' Applies the Halabi-Singh censoring mechanism. This method generates realistic 
#' censoring patterns with control over the censoring rate.
#'
#' @param survival_times Vector of true survival times
#' @param censoring_rate Desired proportion of censored observations (0 to 1)
#' @param censoring_distribution Distribution for censoring times ("exponential" or "uniform")
#' @param seed Random seed for reproducibility
#'
#' @return Data frame with columns:
#'   - time: Observed time (min of survival time and censoring time)
#'   - status: Event indicator (1 = event, 0 = censored)
#'   - true_time: Original survival time (for validation)
#'   - censor_time: Censoring time used
#'
#' @details
#' The Halabi-Singh method supports two censoring mechanisms:
#' 
#' **Exponential Censoring (default):**
#' - Censoring times: C ~ Exponential(β)
#' 
#' **Uniform Censoring:**
#' - Censoring times: C ~ Uniform(0, C)  
#' - May have issues at high censoring rates (>50%)
#'
#'
#' @references
#' Halabi, S., & Singh, B. (2004)
#'
#' @examples
#' \dontrun{
#' # Generate Weibull survival times
#' true_times <- simulate_survival_data(100, "weibull", list(shape = 2, scale = 4))
#' 
#' # Apply 30% exponential censoring 
#' censored_data <- apply_halabi_singh_censoring(true_times, censoring_rate = 0.3)
#' 
#' # Apply uniform censoring (for comparison)
#' censored_uniform <- apply_halabi_singh_censoring(
#'   true_times, censoring_rate = 0.3, censoring_distribution = "uniform"
#' )
#' 
#' # Check actual censoring rates
#' mean(censored_data$status == 0)
#' }
#' @export
apply_halabi_singh_censoring <- function(survival_times, censoring_rate, 
                                     censoring_distribution = "exponential", seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  checkmate::assert_numeric(survival_times, lower = 0, finite = TRUE, min.len = 1)
  checkmate::assert_number(censoring_rate, lower = 0, upper = 1)
  checkmate::assert_choice(censoring_distribution, c("exponential", "uniform"))
  
  n <- length(survival_times)
  
  if (censoring_distribution == "exponential") {
    # Exponential censoring: C ~ Exponential(β)
    find_exponential_parameter <- function(times, target_rate) {
      # Define objective function: compute expected censoring rate for given β
      objective_function <- function(beta) {
        if (beta <= 0) return(1)  # Invalid parameter
        
        # For each survival time, compute P(C < T_i) where C ~ Exp(β)
        # P(C < t) = 1 - exp(-β*t) for exponential 
        censoring_probs <- 1 - exp(-beta * times)
        
        # Expected censoring rate 
        expected_censoring_rate <- mean(censoring_probs)
        
        return(abs(expected_censoring_rate - target_rate))
      }
      
      # Use uniroot to solve for β 
      tryCatch({
        # Find a reasonable range where the function changes sign
        test_betas <- seq(0.001, 5, length.out = 100)
        test_errors <- sapply(test_betas, function(b) {
          probs <- 1 - exp(-b * times)
          mean(probs) - target_rate
        })
        
        # Find where the function crosses zero
        sign_changes <- which(diff(sign(test_errors)) != 0)
        
        if (length(sign_changes) > 0) {
          # Use the first sign change
          idx <- sign_changes[1]
          lower <- test_betas[idx]
          upper <- test_betas[idx + 1]
          
          result <- uniroot(
            function(beta) {
              probs <- 1 - exp(-beta * times)
              mean(probs) - target_rate
            },
            interval = c(lower, upper),
            tol = 1e-6
          )
          
          return(result$root)
        } else {
          # Fallback to optimisation
          result <- optimize(
            f = objective_function,
            interval = c(0.001, 5),
            tol = 1e-6
          )
          return(result$minimum)
        }
      }, error = function(e) {
        # Fallback: use simple approximation
        # For exponential censoring, β ≈ -log(1-p) / mean(times)
        mean_time <- mean(times)
        return(-log(1 - target_rate) / mean_time)
      })
    }
    
    # Find optimal exponential parameter
    beta <- find_exponential_parameter(survival_times, censoring_rate)
    
    # Generate censoring times from exponential distribution
    censor_times <- rexp(n, rate = beta)
    
    # Store parameter for metadata
    censoring_param <- beta
    param_name <- "exponential_rate"
    
  } else {
    # Uniform censoring: C ~ Uniform(0, C)
    # Need to solve: p = E[min(T/C, 1)] = mean(pmin(times/C, 1))
    
    find_uniform_constant <- function(times, target_rate) {
      # Use direct search since we have the exact formula
      search_function <- function(C) {
        if (C <= 0) return(1)  # All censored if C = 0
        expected_censoring_rate <- mean(pmin(times/C, 1))
        return(abs(expected_censoring_rate - target_rate))
      }
      
      # Use optimise to find the best C
      result <- optimize(
        f = search_function,
        interval = c(0.001, max(times) * 3),
        tol = 0.001
      )
      
      return(result$minimum)
    }
    
    # Find optimal uniform constant
    C <- find_uniform_constant(survival_times, censoring_rate)
    
    # Generate censoring times from uniform distribution
    censor_times <- runif(n, 0, C)
    
    # Store parameter for metadata
    censoring_param <- C
    param_name <- "uniform_upper"
  }
  
  # Apply censoring logic 
  # Observed time is minimum of survival time and censoring time
  observed_times <- pmin(survival_times, censor_times)
  
  # Status: 1 = event observed (survival_time <= censor_time)
  #         0 = censored (censor_time < survival_time)
  status <- as.numeric(survival_times <= censor_times)
  
  # Create result data frame
  result <- data.frame(
    time = observed_times,
    status = status,
    true_time = survival_times,
    censor_time = censor_times
  )
  
  # Add metadata
  actual_censoring_rate <- mean(status == 0)
  attr(result, "target_censoring_rate") <- censoring_rate
  attr(result, "actual_censoring_rate") <- actual_censoring_rate
  attr(result, "censoring_distribution") <- censoring_distribution
  attr(result, param_name) <- censoring_param
  
  # Warn if censoring rate is far from target
  if (abs(actual_censoring_rate - censoring_rate) > 0.05) {
    warning(sprintf("Actual censoring rate (%.3f) differs from target (%.3f) by more than 5%% using %s censoring",
                   actual_censoring_rate, censoring_rate, censoring_distribution))
  }
  
  return(result)
}

# Note: Simulation study functions have been removed as they were
# incompatible with the current package API. The functions
# run_simulation_study(), calculate_simulation_metrics(), and 
# calculate_simulation_summary() would need to be completely rewritten
# to work with the unified impute() function and current data structures. 