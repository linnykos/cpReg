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


test_that(".nll is computed correctly",{
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    M <- 5; TT <- 10
    nu <- 0.1*runif(M)
    A <- 0.1*matrix(runif(M*2*M), M, 2*M)

    dat <- generative_model(nu, A, TT, lag = 2)
    transform_dat <- construct_AR_basis(dat, lag = 2)

    res <- .nll(nu, A, dat, transform_dat)

    res2 <- sum(sapply(1:(TT-1), function(time){
      sapply(1:M, function(m){
        exp(nu[m] + A[m,]%*%transform_dat[time,]) - dat[time+1,m] * (nu[m] + A[m,]%*%transform_dat[time,])
      })
    }))

    abs(res - res2) <= 1e-6
  })

  expect_true(all(bool_vec))
})

#######################

## .nll_row is correct

test_that(".nll_row works", {
  set.seed(10)
  M <- 5; TT <- 100
  nu <- 0.1*runif(5)
  A <- 0.1*matrix(runif(M*M), M, M)

  dat <- generative_model(nu, A, TT, lag = 1)
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res <- .nll_row(nu[1], A[1,], dat[,1], transform_dat)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})

test_that(".nll_row sums to .nll", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    M <- 5; TT <- 10
    nu <- 0.1*runif(5)
    A <- 0.1*matrix(runif(M*M), M, M)

    dat <- generative_model(nu, A, TT, lag = 1)
    transform_dat <- construct_AR_basis(dat, lag = 1)

    res <- sum(sapply(1:M, function(x){
      .nll_row(nu[x], A[x,], dat[,x], transform_dat)
    }))
    res2 <- .nll(nu, A, dat, transform_dat)

    abs(res - res2) <= 1e-6
  })

  expect_true(all(bool_vec))
})

#######################

## .objective_func_row is correct

test_that(".objective_func_row works", {
  set.seed(10)
  M <- 5; TT <- 100
  nu <- 0.1*runif(5)
  A <- 0.1*matrix(runif(M*M), M, M)

  dat <- generative_model(nu, A, TT, lag = 1)
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res <- .objective_func_row(nu[1], A[1,], dat[,1], transform_dat, lambda = 1)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})

test_that(".objective_func_row sums to .objective_func", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    M <- 5; TT <- 10
    nu <- 0.1*runif(5)
    A <- 0.1*matrix(runif(M*M), M, M)

    dat <- generative_model(nu, A, TT, lag = 1)
    transform_dat <- construct_AR_basis(dat, lag = 1)

    res <- sum(sapply(1:M, function(x){
      .objective_func_row(nu[x], A[x,], dat[,x], transform_dat, lambda = 1)
    }))
    res2 <- .objective_func(nu, A, dat, transform_dat, lambda = 1)

    abs(res - res2) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".objective_func is roughly minimized under the ground truth",{
  set.seed(1)
  M <- 5; TT <- 100
  nu <- 0.2*runif(5)
  A <- 0.2*matrix(runif(M*M), M, M)
  A[sample(1:prod(dim(A)), round(0.7*prod(dim(A))))] <- 0

  dat <- generative_model(nu, A, TT, lag = 1)
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res <- .objective_func(nu, A, dat, transform_dat, lambda = 1)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    nu_new <- runif(5)
    A_new <- matrix(runif(M*M), M, M)
    res2 <- .objective_func(nu_new, A_new, dat, transform_dat, lambda = 1)

    res < res2
  })

  expect_true(all(bool_vec))
})


test_that(".objective_func is computed correctly",{
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    M <- 5; TT <- 10
    nu <- 0.1*runif(M)
    A <- 0.1*matrix(runif(M*M), M, M)

    dat <- generative_model(nu, A, TT, lag = 1)
    transform_dat <- construct_AR_basis(dat, lag = 1)
    lambda <- 1

    res <- .objective_func(nu, A, dat, transform_dat, lambda = lambda)

    res2 <- sum(sapply(1:(TT-1), function(time){
      sapply(1:M, function(m){
        exp(nu[m] + A[m,]%*%transform_dat[time,]) - dat[time+1,m] * (nu[m] + A[m,]%*%transform_dat[time,])
      })
    }))
    res2 <- res2 + lambda*sum(abs(A))

    abs(res - res2) <= 1e-6
  })

  expect_true(all(bool_vec))
})

######################

## .rowwise_glmnet is correct

test_that(".rowwise_glmnet works", {
  set.seed(10)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  dat <- generative_model(nu, A, TT, lag = 1)
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res <- .rowwise_glmnet(dat[,1], transform_dat)
  expect_true(nrow(res$glmnet.fit$beta) == 5)

  expect_true(class(res) == "cv.glmnet")

  res2 <- .rowwise_glmnet(dat[,1], transform_dat, lambda = 0.25)
  expect_true(is.numeric(res2))
  expect_true(length(res2) == 5+1)
  expect_true(!is.matrix(res2))
})

#########################
