context("Test high dimensional infeasible")

## .reformat_covariates is correct

test_that(".reformat_covariates works", {
  X <- matrix(1:30, 6)
  partition <- c(0,2,5,6)
  res <- .reformat_covariates(X, partition)

  expect_true(is.matrix(res))
  expect_true(!is.vector(res))
  expect_true(all(dim(res) == c(nrow(X), ncol(X)*(length(partition)-1))))
})

test_that(".reformat_covariates gives the correct matrix", {
  set.seed(10)
  X <- matrix(rnorm(30), 6)
  partition <- c(0,2,5,6)
  res <- .reformat_covariates(X, partition)

  # construct the right matrix at each stage
  n <- nrow(X); d <- ncol(X); k <- length(partition)-1
  for(i in 1:k){
    res2 <- matrix(0, n, d)
    res2[(partition[i]+1):partition[i+1],] <- X[(partition[i]+1):partition[i+1],]/sqrt(partition[i+1]-partition[i])
    expect_true(sum(abs(res2 - res[,(d*(i-1)+1):(d*i)])) < 1e-6)
  }
})

###############

## .high_dim_infeasible_subroutine is correct

test_that(".high_dim_infeasible_subroutine works", {
  set.seed(10)
  true_partition <- c(0,0.5,1)
  n <- 100
  dat <- create_data(list(rep(10, 3), rep(-10, 3)), true_partition * n)
  lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
  maxl2 <- .l2norm(rep(10, 3))

  res <- .high_dim_infeasible_subroutine(dat$X, dat$y, lambda, maxl2, true_partition * n)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("obj_val", "coef_list")))
})

#############

## high_dim_infeasible_estimate is correct

test_that("high_dim_infeasible_estimate works", {
  set.seed(10)
  true_partition <- c(0,0.5,1)
  n <- 100
  dat <- create_data(list(rep(10, 3), rep(-10, 3)), true_partition * n)
  lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
  maxl2 <- .l2norm(rep(10, 3))
  K <- 1
  delta <- 10
})
