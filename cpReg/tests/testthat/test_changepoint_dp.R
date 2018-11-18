context("Test changepoint dp")

test_that(".changepoint_dp works", {
  set.seed(10)
  M <- 5; TT <- 25
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  res <- .changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 15,
                         min_spacing = 10, verbose = F)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("obj_val", "partition", "A_list")))
})

test_that(".changepoint_dp reacts properly to gamma", {
  set.seed(10)
  M <- 5; TT <- 50
  nu <- rep(0, M)
  A <- 0.3*diag(M)
  thres_u <- 5
  dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = thres_u)

  res1 <- .changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 1,
                         min_spacing = 10, verbose = F)
  res2 <- .changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 11,
                           min_spacing = 10, verbose = F)
  res3 <- .changepoint_dp(dat, thres_u = thres_u, lambda = 0.05, gamma = 15,
                          min_spacing = 10, verbose = F)

  expect_true(length(res1$partition) >= length(res2$partition))
  expect_true(length(res2$partition) >= length(res3$partition))
})
