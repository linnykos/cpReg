context("Test changepoint dp")

test_that("changepoint_dp works", {
  set.seed(10)
  M <- 5; TT <- 25
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  res <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 15,
                         min_spacing = 10, verbose = F)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("obj_val", "partition", "A_list")))
})

test_that("changepoint_dp reacts properly to gamma", {
  set.seed(10)
  M <- 5; TT <- 50
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  res1 <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                         min_spacing = 10, verbose = F)
  res2 <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 11,
                           min_spacing = 10, verbose = F)
  res3 <- changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 15,
                          min_spacing = 10, verbose = F)

  expect_true(length(res1$partition) >= length(res2$partition))
  expect_true(length(res2$partition) >= length(res3$partition))
})

test_that("changepoint_dp can complete faster with spacing", {
  set.seed(10)
  M <- 5; TT <- 50
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  res1 <- system.time(changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                         min_spacing = 10, verbose = F))
  res2 <- system.time(changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                                     min_spacing = 10, skip_interval = 5, verbose = F))

  expect_true(res1["elapsed"] > res2["elapsed"])
})

################

## generative_model_cp is correct

test_that("generative_model_cp works", {
  set.seed(10)
  M <- 5; TT <- 100
  nu <- rep(0, M)
  A_list <- list(0.75*diag(M),
                 matrix(0, M, M),
                 0.5*matrix(runif(M*M), M, M))
  changepoint_perc <- c(0.4, 0.6)

  res <- generative_model_cp(nu, A_list, changepoint_perc, timesteps = TT)

  expect_true(is.list(res))
  expect_true(class(res) == "SEPP_cp")
})
