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


  res <- .find_breakpoint(fit, interval, compute_cusum_func = .compute_cusum,
                          data_length_func = nrow)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("val", "b")))
})

test_that(".find_breakpoint works for really small intervals", {
  fit <- matrix(0, nrow = 200, ncol = 3)
  fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))
  interval <- c(100,102)

  res <- .find_breakpoint(fit, interval)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("val", "b")))
  expect_true(all(is.na(unlist(res))))
})

################

## .remove_duplicate is correct

test_that(".remove_duplicate works", {
  res <- .remove_duplicate(list(c(0,1), c(0,5), c(0,1)))

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(unlist(res) == c(0,1,0,5)))
})

###################

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

################

## screening is correct

test_that("screening works", {
  fit <- matrix(0, nrow = 200, ncol = 3)
  fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))

  res <- screening(fit, tau = 10, verbose = F)

  expect_true(is.numeric(res))
  expect_true(all(res == c(0, 100, 200)))
})

test_that("screening works for a more complex setting", {
  set.seed(10)
  fit <- matrix(0, nrow = 300, ncol = 3)
  fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))
  fit[201:300,] <- c(rep(1,100), rep(0, 100), rep(1,100))
  fit <- fit+rnorm(900, sd = 0.1)

  res <- screening(fit, tau = 10, verbose = F)

  expect_true(is.numeric(res))
  expect_true(all(res == c(0, 100, 200, 300)))
})
