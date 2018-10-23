context("Test basis")

test_that("construct_AR_basis works", {
  dat <- matrix(1:32, ncol = 4, nrow = 8)

  res <- construct_AR_basis(dat)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(8,4)))
})

test_that("construct_AR_basis works forms the correct basis for lag = 1", {
  dat <- matrix(1:32, ncol = 4, nrow = 8)

  res <- construct_AR_basis(dat)

  expect_true(all(res[1,] == 0))
  expect_true(all(res[-1,] == dat[-8,]))
})

test_that("construct_AR_basis works forms the correct basis for lag = 2", {
  dat <- matrix(1:32, ncol = 4, nrow = 8)

  res <- construct_AR_basis(dat, lag = 2)

  expect_true(all(dim(res) == c(8,8)))
  expect_true(all(res[1,] == 0))
  expect_true(all(res[2,5:8] == 0))
  expect_true(all(res[-1,1:4] == dat[-8,]))
  expect_true(all(res[-c(1:2),5:8] == dat[-c(7:8),]))
})

test_that("construct_AR_basis works forms the correct basis for lag = 2", {
  dat <- matrix(1:32, ncol = 4, nrow = 8)

  res <- construct_AR_basis(dat, lag = 2)

  expect_true(all(dim(res) == c(8,8)))
  expect_true(all(res[1,] == 0))
  expect_true(all(res[2,5:8] == 0))
  expect_true(all(res[-1,1:4] == dat[-8,]))
  expect_true(all(res[-c(1:2),5:8] == dat[-c(7:8),]))
})

test_that("construct_AR_basis can set a threshold", {
  dat <- matrix(1:32, ncol = 4, nrow = 8)

  res <- construct_AR_basis(dat, lag = 2, thres_u = 2)

  expect_true(all(res <= 2))
})

