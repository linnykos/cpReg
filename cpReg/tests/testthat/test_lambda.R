context("Test lambda")

## .extract_lambdas is correct

test_that(".extract_lambdas works", {
  set.seed(1)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  obj <- generative_model(nu, A, TT, lag = 1)
  dat <- obj$dat
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, cv = F)
  })

  res <- .extract_lambdas(res_list, lambda = 0.1)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("lambda", "nu", "A")))
  expect_true(length(res$lambda) == 1)
  expect_true(length(res$nu) == M)
  expect_true(all(dim(res$A) == c(M,M)))
})

test_that(".extract_lambdas implies that penalization doesn't affect intercepts", {
  set.seed(10)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  obj <- generative_model(nu, A, TT, lag = 1)
  dat <- obj$dat
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, cv = F)
  })

  res1 <- .extract_lambdas(res_list, lambda = 5)
  res2 <- .extract_lambdas(res_list, lambda = 0.05)

  expect_true(length(which(res1$nu != 0)) == length(which(res2$nu != 0)))
  expect_true(length(which(res1$A != 0)) < length(which(res2$A != 0)))
})

########################

## .combine_lambdas is correct

test_that(".combine_lambdas works", {
  set.seed(10)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  obj <- generative_model(nu, A, TT, lag = 1)
  dat <- obj$dat
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, cv = T)
  })

  res <- .combine_lambdas(res_list, dat, transform_dat, verbose = F)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("lambda", "nu", "A")))
  expect_true(length(res$lambda) == 1)
  expect_true(length(res$nu) == M)
  expect_true(all(dim(res$A) == c(M,M)))
})

test_that(".combine_lambdas can fit properly, based on .extract_lambdas", {
  set.seed(5)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  obj <- generative_model(nu, A, TT, lag = 1)
  dat <- obj$dat
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, cv = T)
  })

  res <- .combine_lambdas(res_list, dat, transform_dat, verbose = F)

  res2_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, cv = F)
  })
  res2 <- .extract_lambdas(res2_list, lambda = res$lambda)

  expect_true(sum(abs(res$nu - res2$nu)) <= 1e-6)
  expect_true(sum(abs(res$A - res2$A)) <= 1e-6)
})

###################

## .lambda_oracle is correct

test_that(".lambda_oracle works", {
  set.seed(10)
  M <- 5; TT <- 300
  nu <- rep(0, M)
  A_list <- list(0.75*diag(M),
                 matrix(0, M, M),
                 0.5*matrix(runif(M*M), M, M))
  changepoint_perc <- c(0.4, 0.6)

  obj <- generative_model_cp(nu, A_list, changepoint_perc, timesteps = TT)

  res <- .lambda_oracle(obj, intercept = F)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
})
