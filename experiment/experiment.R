rm(list=ls())
set.seed(10)
dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), c(0, 50, 100))

true_partition <-  c(0, .5, 1)
lambda <- oracle_tune_lambda(dat$X, dat$y, true_partition)
delta <- 10
# gamma_range <- oracle_tune_gamma_range(dat$X, dat$y, lambda = lambda, k = length(true_partition)-1, delta = delta,
#                                        verbose = T)
X = dat$X
y = dat$y
K = length(true_partition)-1
verbose = T
min_gamma = 0.01
max_gamma = 1000
max_iter = 10

high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = 1700.634,
                           delta = delta)

# min_gamma <- .initial_gamma_overshoot(X, y, lambda, K, min_gamma, delta = delta, smaller = T,
#                                       verbose = verbose)
# high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = min_gamma,
#                            delta = delta)
# max_gamma <- .initial_gamma_overshoot(X, y, lambda, K, max_gamma, delta = delta, smaller = F,
#                                       verbose = verbose)
# high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = max_gamma,
#                            delta = delta)
