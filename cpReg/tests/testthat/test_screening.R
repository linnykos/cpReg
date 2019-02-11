context("Test screening")

## .compute_cusum is correct

test_that(".compute_cusum works", {
  fit <- matrix(0, nrow = 200, ncol = 3)
  fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))
  res <- .compute_cusum(fit, 0, 200, 100)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
})

test_that(".compute_cusum is maximized at the right location", {
  fit <- matrix(0, nrow = 200, ncol = 3)
  fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))
  res <- sapply(1:199, .compute_cusum, fit = fit, s = 0, e = 200)

  expect_true(which.max(res) == 100)
})

###############

## .find_breakpoint is correct

test_that(".find_breakpoint works", {
  fit <- matrix(0, nrow = 200, ncol = 3)
  fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))
  interval <- c(0, 200)

  res <- .find_breakpoint(fit, interval)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("val", "b")))
})

################

## .generate_intervals is correct

test_that(".generate_intervals works", {
  res <- .generate_intervals(200, 100)

  expect_true(length(res) == 101)
  expect_true(is.list(res))
  expect_true(all(sapply(res, function(x){x[1] < x[2]})))
})

################

## .truncate is correct

test_that(".truncate works", {
  set.seed(10)
  random_intervals <- .generate_intervals(200, 100)
  res <- .truncate(c(100,120), random_intervals)

  expect_true(is.list(res))
  expect_true(all(sapply(res, function(x){x[1] < x[2]})))
})

test_that(".truncate can return the same intervals", {
  set.seed(10)
  random_intervals <- .generate_intervals(200, 100)
  res <- .truncate(c(1,200), random_intervals)

  expect_true(all(sapply(1:101, function(x){all(res[[x]] == random_intervals[[x]])})))
})

test_that(".truncate all lie in the range", {
  set.seed(10)
  random_intervals <- .generate_intervals(200, 100)
  res <- .truncate(c(100,120), random_intervals)

  expect_true(all(sapply(res, function(x){x[1]>=100})))
  expect_true(all(sapply(res, function(x){x[2]<=120})))
})
