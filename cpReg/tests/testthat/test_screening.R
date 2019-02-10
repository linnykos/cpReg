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

