context("Test Buhlmann estimator")

## .compute_regression_buhlmann is correct

test_that(".compute_regression_buhlmann works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .compute_regression_buhlmann(dat, start = 1,
                                   end = 100, lambda = 1, breakpoint = 50,
                                   gamma = 10)

  expect_true(is.numeric(res))
})

############

## .buhlmann_threshold_closure is correct

test_that(".buhlmann_threshold works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .buhlmann_threshold(dat, c(1,100), 10, 10)

  expect_true(is.numeric(res))
})

#############

## high_dim_buhlmann_estimate is correct

test_that("high_dim_buhlmann_estimate works", {
  set.seed(10)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))
  res <- high_dim_buhlmann_estimate(dat$X, dat$y, 0.1, 1)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
})

###############

## oracle_tune_gamma_range is correct

test_that("oracle_tune_gamma_range works", {
  set.seed(10)
  n <- 100
  partition <- c(0, 0.5, 1)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), round(partition*n))
  lambda <- oracle_tune_lambda(dat$X, dat$y, partition)
  k <- length(partition)-1

  res <- oracle_tune_gamma_range(dat$X, dat$y, lambda, k, verbose = F)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("gamma", "min_gamma", "max_gamma")))
})
