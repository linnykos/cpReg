context("Test Buhlmann estimator")

## .compute_regression_buhlmann is correct

test_that(".compute_regression_buhlmann works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .compute_regression_buhlmann(dat, start = 1,
                                   end = 100, lambda = 1, breakpoint = 50)

  expect_true(is.numeric(res))
})

############

## high_dim_buhlmann_estimate is correct

test_that("high_dim_buhlmann_estimate works", {
  set.seed(10)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))
  res <- high_dim_buhlmann_estimate(dat$X, dat$y, lambda = 0.1, K = 2)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
})
