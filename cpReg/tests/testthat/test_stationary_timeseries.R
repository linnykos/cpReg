context("Test stationary timeseries")

## generative_model is correct

test_that("generative_model works", {
  set.seed(1)
  M <- 5; TT <- 10
  nu <- 0.1*runif(5)
  A <- 0.1*matrix(runif(M*M), M, M)

  res <- generative_model(nu, A, TT, lag = 1)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(TT,M)))
})

test_that("generative_model works for lag = 2", {
  set.seed(1)
  M <- 5; TT <- 10
  nu <- 0.1*runif(5)
  A <- 0.1*matrix(runif(M*2*M), nrow = M, ncol = 2*M)

  res <- generative_model(nu, A, TT, lag = 2)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(TT,M)))
})

test_that("generative_model prevents random data from blowing up", {
  set.seed(1)
  M <- 5; TT <- 10
  nu <- runif(5)
  A <- matrix(runif(M*M), nrow = M, ncol = M)

  res <- generative_model(nu, A, TT, lag = 1, thres_u = 1)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(TT,M)))
})

test_that("generative_model stops if parameters blow up", {
  set.seed(1)
  M <- 5; TT <- 10
  nu <- runif(5)
  A <- matrix(runif(M*2*M), nrow = M, ncol = 2*M)

  expect_error(generative_model(nu, A, TT, lag = 1))
})


###################


## .nll is correct

test_that(".nll works", {
  set.seed(10)
  M <- 5; TT <- 10
  nu <- runif(5)
  A <- matrix(runif(M*M), M, M)
  dat <- matrix(stats::rpois(M*TT, lambda = 1), nrow = TT, ncol = M)
  transform_dat <- construct_AR_basis(dat)

  res <- .nll(nu, A, dat, transform_dat)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})

test_that(".nll is roughly minimized under the ground truth",{
  set.seed(1)
  M <- 5; TT <- 100
  nu <- 0.1*runif(5)
  A <- 0.1*matrix(runif(M*M), M, M)

  dat <- generative_model(nu, A, TT, lag = 1)

  transform_dat <- construct_AR_basis(dat, lag = 1)
  res <- .nll(nu, A, dat, transform_dat)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    nu_new <- runif(5)
    A_new <- matrix(runif(M*M), M, M)
    res2 <- .nll(nu_new, A_new, dat, transform_dat)

    res < res2
  })

  expect_true(all(bool_vec))
})
