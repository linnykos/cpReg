context("Test high dimensional feasible estimator")

## .lasso_regression is correct

test_that(".lasso_regression works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .lasso_regression(dat$X, dat$y, lambda = 1)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 3)
})

##########

## .compute_regression_cusum is correct

test_that(".compute_regression_cusum works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .compute_regression_cusum(dat, start = 1,
                           end = 100, lambda = 1, breakpoint = 50)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that(".compute_regression_cusum is maximized appropriately", {
  set.seed(10)
  dat1 <- create_data(list(c(-10,-10,-10), c(20,20,20)), c(0, 50, 100))
  dat2 <- create_data(list(c(-10,-10,-10), c(20,20,20)), c(0, 75, 100))
  res1 <- sapply(5:95, function(x){
    .l2norm(.compute_regression_cusum(dat1, start = 1,
                      end = 100, lambda = 1, breakpoint = x))
  })
  res2 <- sapply(5:95, function(x){
    .l2norm(.compute_regression_cusum(dat2, start = 1,
                              end = 100, lambda = 1, breakpoint = x))
  })

  expect_true(abs(which.max(res1) - 45) < abs(which.max(res2) - 45))
  expect_true(abs(which.max(res2) - 70) < abs(which.max(res1) - 70))
})

############

## high_dim_feasible_estimate is correct

test_that("high_dim_feasible_estimate works" ,{
  set.seed(10)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))
  res <- high_dim_feasible_estimate(dat$X, dat$y, 0.1, 40, M = 1)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
})

test_that("high_dim_feasible_estimate exits gracefully when delta is large" ,{
  vec <- c("n" = 100, "X_type" = 1, "d/n" = 0.5)
  X_type_vec <- c("identity", "toeplitz", "equicorrelation")
  true_partition <- c(0,0.3,0.7,1)

  create_coef <- function(vec, full = F){
    d <- 50
    beta1 <- c(rep(1, 10), rep(0, d-10))
    beta2 <- c(rep(0, d-10), rep(1, 10))
    lis <- list(beta1 = beta1, beta2 = beta2)

    if(!full){
      lis
    } else {
      mat <- matrix(0, nrow = vec["n"], ncol = vec["d/n"]*vec["n"])
      idx <- round(true_partition*vec["n"])
      for(i in 1:(length(idx)-1)){
        zz <- i %% 2; if(zz == 0) zz <- 2
        mat[(idx[i]+1):idx[i+1],] <- rep(lis[[zz]], each = idx[i+1]-idx[i])
      }
      mat
    }
  }

  rule <- function(vec){
    lis <- create_coef(vec, full = F)

    cpReg::create_data(list(lis$beta1, lis$beta2, lis$beta1), round(true_partition*vec["n"]),
                       cov_type = X_type_vec[vec["X_type"]])
  }

  set.seed(1)
  dat <- rule(vec)
  lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
  tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                                factor = 3/4)
  delta <- max(round(vec["n"]/10), 10)

  res <- cpReg::high_dim_feasible_estimate(dat$X, dat$y, lambda = lambda, tau = tau,
                                            verbose = F, max_candidates = NA, delta = delta, M = 0)

  expect_true(is.list(res))
})

################

## oracle_tune_lambda is correct

test_that("oracle_tune_lambda works" ,{
  set.seed(10)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))
  lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, c(0,0.5,1))

  expect_true(is.numeric(lambda))
  expect_true(lambda > 0)
})

#################

## oracle_tune_tau is correct

test_that("oracle_tune_tau works",{
  set.seed(20)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))
  lambda <- oracle_tune_lambda(dat$X, dat$y, c(0,0.5,1))
  tau <- oracle_tune_tau(dat$X, dat$y, lambda, c(0,0.5,1),
                         factor = 7/8)

  expect_true(is.numeric(tau))
  expect_true(tau > 0)
})

