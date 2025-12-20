test_that("normalize_mcmc_options merges missing defaults", {
  opts <- normalize_mcmc_options(list(iter_warmup = 50, iter_sampling = 50, chains = 2))

  expect_equal(opts$iter_warmup, 50L)
  expect_equal(opts$iter_sampling, 50L)
  expect_equal(opts$chains, 2L)

  expect_true(is.numeric(opts$adapt_delta))
  expect_length(opts$adapt_delta, 1)
  expect_true(is.finite(opts$adapt_delta))
  expect_true(opts$adapt_delta > 0 && opts$adapt_delta < 1)

  expect_true(is.integer(opts$max_treedepth))
  expect_length(opts$max_treedepth, 1)
  expect_true(is.finite(opts$max_treedepth))
  expect_true(opts$max_treedepth > 0)
})

test_that("build_rstan_control always returns scalars", {
  control <- build_rstan_control(list(iter_warmup = 50, iter_sampling = 50, chains = 2))

  expect_true(is.list(control))
  expect_true(all(c("adapt_delta", "max_treedepth") %in% names(control)))

  expect_true(is.numeric(control$adapt_delta))
  expect_length(control$adapt_delta, 1)
  expect_true(is.finite(control$adapt_delta))

  expect_true(is.integer(control$max_treedepth))
  expect_length(control$max_treedepth, 1)
  expect_true(is.finite(control$max_treedepth))
})

