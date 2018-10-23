context("Test lambda")

## .extract_lambdas is correct

test_that(".extract_lambdas works", {
  set.seed(1)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  dat <- generative_model(nu, A, TT, lag = 1)
  transform_dat <- construct_AR_basis(dat, lag = 1)

  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, cv = F)
  })

  res <- .extract_lambdas(res_list, lambda = 0.1)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("nu", "A")))
  expect_true(length(res$nu) == M)
  expect_true(all(dim(res$A) == c(M,M)))
})

test_that(".extract_lambdas implies that penalization doesn't affect intercepts", {
  set.seed(1)
  M <- 5; TT <- 100
  nu <- 0.1*runif(M)
  A <- 0.1*matrix(runif(M*M), M, M)

  dat <- generative_model(nu, A, TT, lag = 1)
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
