context("Test gamma")

## .gamma_oracle is correct

test_that(".gamma_oracle works", {
  set.seed(10)
  M <- 5; TT <- 5000
  nu <- rep(0, M)
  A_list <- list(0.75*diag(M),
                 matrix(0, M, M),
                 0.5*matrix(runif(M*M), M, M))
  changepoint_perc <- c(0.3, 0.6)

  obj <- generative_model_cp(nu, A_list, changepoint_perc, timesteps = TT)

  lambda <- .lambda_oracle(obj, intercept = F)
  res <- .gamma_oracle(obj, lambda, intercept = F, min_spacing = 100)

  zz <-  changepoint_dp(obj$dat, lambda = lambda, gamma = res,
                        min_spacing = 50, skip_interval = 500, verbose = T)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
})
