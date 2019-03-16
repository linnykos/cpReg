context("Test wild binary segmentation")

## .find_breakpoint is correct

test_that(".find_breakpoint works", {
  set.seed(10)
  n <- 100
  true_partition <- c(0, 0.5, 1)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), true_partition*n)
  lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
  tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                                factor = 1/2)
  compute_cusum_func <- .compute_regression_cusum
  tau_function <- function(data, interval, ...){
    tau
  }
  M <- 100
  delta <- 10
  max_candidates <- NA

  res <- .find_breakpoint(dat, c(0, 100), delta = delta, max_candidates = max_candidates,
                          data_length_func = function(x){nrow(x$X)},
                          compute_cusum_func = compute_cusum_func, verbose = F,
                          lambda = lambda)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("val", "b")))
})

test_that(".find_breakpoint exits gracefully when delta is too large", {
  set.seed(10)
  n <- 1000
  true_partition <- c(0, 0.5, 1)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), true_partition*n)
  lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
  tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                                factor = 1/2)
  compute_cusum_func <- .compute_regression_cusum
  tau_function <- function(data, interval, ...){
    tau
  }
  delta <- 100
  max_candidates <- 10

  res <- .find_breakpoint(dat, c(723, 1000), delta = delta, max_candidates = max_candidates,
                          data_length_func = function(x){nrow(x$X)},
                          compute_cusum_func = compute_cusum_func, verbose = F,
                          lambda = lambda)

  expect_true(is.list(res))
})

