rm(list=ls())
set.seed(10)
d <- 50; n <- 100; snr <- 10
beta1 <- c(rep(snr, 2), rep(0, d-2))
beta2 <- c(rep(0, d-2), rep(snr, 2))
true_partition <- c(0,0.25,0.75,1)
dat <- create_data(list(beta1, beta2, beta1),
                   round(true_partition*n),
                   cov_type = "identity")

lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
tau <- oracle_tune_tau(dat$X, dat$y, lambda, true_partition)

res <- high_dim_feasible_estimate(dat$X, dat$y, lambda = lambda, tau = tau,
                                  verbose = T)
