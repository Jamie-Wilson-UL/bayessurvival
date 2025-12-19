test_that("get_default_priors returns correct structure", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("writexl")
  priors <- get_default_priors("weibull")
  
  expect_true(is.list(priors))
  expect_true(all(c("mu_log_shape", "sd_log_shape", 
                   "mu_log_scale", "sd_log_scale") %in% names(priors)))
  
  # Check reasonable default values
  expect_equal(priors$mu_log_shape, 0)
  expect_true(priors$sd_log_shape > 0)
  expect_equal(priors$mu_log_scale, 0)
  expect_true(priors$sd_log_scale > 0)
})

test_that("get_default_priors handles unknown distribution", {
  expect_error(
    get_default_priors("unknown"),
    "Unknown distribution: unknown"
  )
})

test_that("stan_models_ready function works", {
  # Should return FALSE initially since models aren't compiled yet
  # (This may change as we develop)
  result <- stan_models_ready()
  expect_true(is.logical(result))
})

# Note: We skip actual Stan compilation tests in automated testing
# as they require a working C++ toolchain and can be slow.
test_that("get_stan_model handles missing model gracefully", {
  skip_if_not_installed("rstan")
  skip("Skipping Stan compilation test - requires rstan toolchain")
  
  # This should trigger compilation attempt
  expect_error(
    get_stan_model("nonexistent"),
    "Stan model file not found"
  )
}) 
