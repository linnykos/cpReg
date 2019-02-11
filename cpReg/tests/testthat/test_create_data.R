context("Test create data")

## create_data is correct

test_that("create_data works", {
  set.seed(10)
  res <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))

  expect_true(is.list(res))
  expect_true(all(names(res) == c("X", "y")))
  expect_true(all(dim(res$X) == c(100, 3)))
  expect_true(length(res$y) == 100)
})

test_that("create_data can create difference covariance", {
  set.seed(10)
  res1 <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 500, 1000), cov_type = "identity")
  res2 <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 500, 1000), cov_type = "toeplitz")
  res3 <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 500, 1000), cov_type = "equicorrelation")

  cov1 <- diag(3)
  cov2 <- stats::toeplitz(0.8^(0:2))
  cov3 <- matrix(0.2, 3, 3); diag(cov3) <- 1

  expect_true(sum(abs(stats::cov(res1$X) - cov1)) < sum(abs(stats::cov(res1$X) - cov2)))
  expect_true(sum(abs(stats::cov(res1$X) - cov1)) < sum(abs(stats::cov(res1$X) - cov3)))

  expect_true(sum(abs(stats::cov(res2$X) - cov2)) < sum(abs(stats::cov(res2$X) - cov1)))
  expect_true(sum(abs(stats::cov(res2$X) - cov2)) < sum(abs(stats::cov(res2$X) - cov3)))

  expect_true(sum(abs(stats::cov(res3$X) - cov3)) < sum(abs(stats::cov(res3$X) - cov1)))
  expect_true(sum(abs(stats::cov(res3$X) - cov3)) < sum(abs(stats::cov(res3$X) - cov2)))
})

