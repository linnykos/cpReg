context("Test changepoint screening")

## .unravel is correct
test_that(".unravel works", {
  s <- 3
  A_list <- lapply(1:s, function(x){
    matrix((1:25)+(x-1)*25, 5, 5)
  })
  partition <- c(0,5,10,15)

  res <- .unravel(partition, A_list)

  expect_true(is.list(res))
  expect_true(length(res) == max(partition))
})

test_that(".unravel gets the right enumeration", {
  s <- 3
  A_list <- lapply(1:s, function(x){
    matrix((1:25)+(x-1)*25, 5, 5)
  })
  partition <- c(0,5,10,15)

  res <- .unravel(partition, A_list)

  diff_vec <- sapply(2:max(partition), function(x){
    sum(abs(Reduce("-", res[(x-1):x])))
  })

  expect_true(length(which(diff_vec != 0)) == length(partition)-2)

  expect_true(all(diff_vec[partition[2:(length(partition)-1)]] != 0))
  expect_true(all(diff_vec[-partition[2:(length(partition)-1)]] == 0))
})

#####################

## .multidimen_cusum is works

test_that(".multidimen_cusum works", {
  s <- 3
  A_list <- lapply(1:s, function(x){
    matrix((1:25)+(x-1)*25, 5, 5)
  })
  partition <- c(0,5,10,15)

  A_long <- .unravel(partition, A_list)
  res <- .multidimen_cusum(A_long[1:5], A_long[6:10])

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(A_list[[1]])))
})

test_that(".multidimen_cusum can output varied numbers", {
  set.seed(10)
  s <- 3
  A_list <- lapply(1:s, function(x){
    matrix(s+rnorm(25), 5, 5)
  })
  partition <- c(0,5,10,15)

  A_long <- .unravel(partition, A_list)
  res <- .multidimen_cusum(A_long[1:5], A_long[6:10])

  expect_true(length(unique(as.numeric(res))) == prod(dim(res)))
})

test_that(".multidimen_cusum is maximized appropriately at a changepoint", {
  set.seed(10)
  s <- 2
  A_list <- lapply(1:s, function(x){
    matrix(10*s+rnorm(25), 5, 5)
  })
  partition <- c(0,5,10)

  A_long <- .unravel(partition, A_list)

  res_vec <- sapply(2:8, function(x){
    .l2norm(.multidimen_cusum(A_long[1:x], A_long[(x+1):10]))
  })

  idx <- which(2:8 == partition[2])
  expect_true(which.max(res_vec) == idx)
})

###################

## changepoint_screening is correct

test_that("changepoint_screening works", {
  set.seed(10)
  M <- 5; TT <- 25
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  fit <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                         min_spacing = 10, verbose = F)

  res <- changepoint_screening(fit, gamma = 2)

  expect_true(is.list(res))
})

test_that("changepoint_screening can remove changepoints", {
  set.seed(10)
  M <- 5; TT <- 25
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  fit <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                        min_spacing = 10, verbose = F)

  res <- changepoint_screening(fit, gamma = 5)

  expect_true(length(res$partition) < length(fit$partition))
})


test_that("changepoint_screening can handle harder cases", {
  set.seed(10)
  M <- 5; TT <- 50
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  fit <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                         min_spacing = 10, verbose = F)

  res <- changepoint_screening(fit, gamma = 5)

  expect_true(length(res$partition) == 2)
})
