context("Test Buhlmann estimator")

## .compute_regression_buhlmann is correct

test_that(".compute_regression_buhlmann works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .compute_regression_buhlmann(dat, start = 1,
                                   end = 100, lambda = 1, breakpoint = 50)

  expect_true(is.numeric(res))
})

############

## high_dim_buhlmann_estimate is correct

test_that("high_dim_buhlmann_estimate works", {
  set.seed(10)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))
  res <- high_dim_buhlmann_estimate(dat$X, dat$y, lambda = 0.1, K = 2)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
})

test_that("high_dim_buhlmann_estimate works with given gamma", {
  set.seed(10)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))

  true_partition <-  c(0, .5, 1)
  lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
  delta <- 10
  gamma_range <- oracle_tune_gamma_range(dat$X, dat$y, lambda = lambda,
                                         min_gamma = 1000, max_gamma= 25000,
                                         K = length(true_partition)-1, delta = delta,
                                         verbose = F)

  if(any(is.na(gamma_range$gamma))) gamma <- gamma_range$min_gamma else gamma <- mean(gamma_range$gamma)

  res <- high_dim_buhlmann_estimate(dat$X, dat$y, lambda = lambda, gamma = gamma)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
})

############

## .buhlmann_threshold is correct

test_that(".buhlmann_threshold works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .buhlmann_threshold(dat, c(1,100), 10, 10)

  expect_true(is.numeric(res))
})


###############

## oracle_tune_gamma_range is correct

test_that("oracle_tune_gamma_range works", {
  set.seed(10)
  n <- 100
  partition <- c(0, 0.5, 1)
  dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), round(partition*n))
  lambda <- oracle_tune_lambda(dat$X, dat$y, partition)
  K <- length(partition)-1

  res <- oracle_tune_gamma_range(dat$X, dat$y, lambda, K = K, verbose = F)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("gamma", "min_gamma", "max_gamma")))
})


test_that("oracle_tune_gamma_range can select appropriate gamma", {
  paramMat <- as.matrix(expand.grid(round(exp(seq(log(100), log(1000), length.out = 10))), c(1,2),
                                    1/2))
  colnames(paramMat) <- c("n", "X_type", "d/n")

  X_type_vec <- c("identity", "toeplitz", "equicorrelation")
  true_partition <- c(0,0.3,0.7,1)

  create_coef <- function(vec, full = F){
    d <- 50 # d <- vec["d/n"]*vec["n"]
    beta1 <- c(rep(1, 10), rep(0, d-10))
    beta2 <- c(rep(0, d-10), rep(1, 10))
    lis <- list(beta1 = beta1, beta2 = beta2)

    if(!full){
      lis
    } else {
      mat <- matrix(0, nrow = vec["n"], ncol = d)
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

    create_data(list(lis$beta1, lis$beta2, lis$beta1), round(true_partition*vec["n"]),
                cov_type = X_type_vec[vec["X_type"]])
  }

  #####

  set.seed(49)
  vec <- paramMat[1,]
  dat <- rule(vec)

  true_beta <- create_coef(vec, full = T)

  delta <- max(round(vec["n"]/10), 10)
  lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
  res <- oracle_tune_gamma_range(dat$X, dat$y, lambda = lambda, K = length(true_partition)-1, delta = delta,
                                         verbose = F)

  expect_true(length(res$gamma) == 2)

})
