# Tests for the complete() function

test_that("complete() extracts single dataset correctly", {
  # Create mock bayesian_imputation object for testing
  mock_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )
  
  # Create mock imputed datasets
  dataset1 <- mock_data
  dataset1$time[c(2,4)] <- c(18, 28)  # Replace censored times
  dataset1$status[c(2,4)] <- 1
  dataset1$dataset_id <- 1
  
  dataset2 <- mock_data  
  dataset2$time[c(2,4)] <- c(19, 29)
  dataset2$status[c(2,4)] <- 1
  dataset2$dataset_id <- 2
  
  mock_result <- structure(list(
    original_data = mock_data,
    time_col = "time",
    status_col = "status",
    imputed_datasets = list(dataset1, dataset2)
  ), class = "bayesian_imputation")
  
  # Test extracting single dataset
  result <- complete(mock_result, dataset = 1)
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true(all(result$status == 1))  # All should be events now
  expect_true(".imp" %in% names(result))  # Should have .imp column
  expect_equal(result$.imp[1], 1)  # Should be dataset 1
})

test_that("complete() extracts all datasets correctly", {
  # Create mock data (same as above)
  mock_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )
  
  dataset1 <- mock_data
  dataset1$time[c(2,4)] <- c(18, 28)
  dataset1$status[c(2,4)] <- 1
  dataset1$dataset_id <- 1
  
  dataset2 <- mock_data
  dataset2$time[c(2,4)] <- c(19, 29)  
  dataset2$status[c(2,4)] <- 1
  dataset2$dataset_id <- 2
  
  mock_result <- structure(list(
    original_data = mock_data,
    time_col = "time", 
    status_col = "status",
    imputed_datasets = list(dataset1, dataset2)
  ), class = "bayesian_imputation")
  
  # Test extracting all datasets
  result <- complete(mock_result)
  
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(names(result), c("dataset_1", "dataset_2"))
  expect_true(all(sapply(result, is.data.frame)))
})

test_that("complete() produces long format correctly", {
  # Create simple mock data
  mock_data <- data.frame(
    time = c(10, 15),
    status = c(1, 0)
  )
  
  dataset1 <- data.frame(time = c(10, 18), status = c(1, 1), dataset_id = 1)
  dataset2 <- data.frame(time = c(10, 19), status = c(1, 1), dataset_id = 2)
  
  mock_result <- structure(list(
    original_data = mock_data,
    time_col = "time",
    status_col = "status", 
    imputed_datasets = list(dataset1, dataset2)
  ), class = "bayesian_imputation")
  
  # Test long format
  result <- complete(mock_result, format = "long")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 4)  # 2 datasets × 2 observations
  expect_true(all(c(".imp", ".id") %in% names(result)))
  expect_equal(unique(result$.imp), c(1, 2))
  expect_equal(result$.id, c(1, 2, 1, 2))  # Row IDs within each dataset
})

test_that("complete() validates inputs correctly", {
  
  # Test with non-bayesian_imputation object
  expect_error(
    complete(list(a = 1)),
    "no applicable method for 'complete' applied"
  )
  
  # Create mock object without imputed datasets
  mock_empty <- structure(list(
    original_data = data.frame(time = 1, status = 1),
    imputed_datasets = NULL
  ), class = "bayesian_imputation")
  
  expect_error(
    complete(mock_empty),
    "No imputed datasets found"
  )
  
  # Create mock object for testing dataset number validation
  mock_data <- data.frame(time = c(10, 15), status = c(1, 0))
  mock_result <- structure(list(
    original_data = mock_data,
    imputed_datasets = list(mock_data, mock_data)  # 2 datasets
  ), class = "bayesian_imputation")
  
  # Test invalid dataset number
  expect_error(
    complete(mock_result, dataset = 5),
    "Dataset numbers must be between 1 and 2"
  )
  
  expect_error(
    complete(mock_result, dataset = 0),
    "Dataset numbers must be between 1 and 2"
  )
})

test_that("complete() handles include.original correctly", {
  mock_data <- data.frame(
    time = c(10, 15),
    status = c(1, 0)  # One censored observation
  )
  
  dataset1 <- data.frame(time = c(10, 18), status = c(1, 1), dataset_id = 1)
  
  mock_result <- structure(list(
    original_data = mock_data,
    time_col = "time",
    status_col = "status",
    imputed_datasets = list(dataset1)
  ), class = "bayesian_imputation")
  
  # Test including original data
  result <- complete(mock_result, include.original = TRUE)
  
  expect_true(is.list(result))
  expect_equal(length(result), 2)  # Original + 1 imputed
  expect_equal(names(result), c("dataset_0", "dataset_1"))
  
  # Check that dataset_0 is the original data
  original_extracted <- result$dataset_0
  expect_equal(original_extracted$.imp[1], 0)
  expect_true(any(original_extracted$status == 0))  # Should have censoring
})

test_that("complete() works with different format arguments", {
  mock_data <- data.frame(time = c(10, 15), status = c(1, 0))
  dataset1 <- data.frame(time = c(10, 18), status = c(1, 1), dataset_id = 1)
  
  mock_result <- structure(list(
    original_data = mock_data,
    time_col = "time",
    status_col = "status", 
    imputed_datasets = list(dataset1)
  ), class = "bayesian_imputation")
  
  # Wide format returns a list for dataset = "all" (default)
  result_wide_all <- complete(mock_result, format = "wide")
  expect_true(is.list(result_wide_all))
  expect_equal(length(result_wide_all), 1)
  expect_true(is.data.frame(result_wide_all$dataset_1))

  # A single numeric dataset returns a data.frame
  result_wide_one <- complete(mock_result, dataset = 1, format = "wide")
  expect_true(is.data.frame(result_wide_one))
  
  # Test list format  
  result_list <- complete(mock_result, format = "list")
  expect_true(is.list(result_list))
  expect_true(is.data.frame(result_list$dataset_1))
  
  # Test long format
  result_long <- complete(mock_result, format = "long")
  expect_true(is.data.frame(result_long))
  expect_true(".imp" %in% names(result_long))
}) 
