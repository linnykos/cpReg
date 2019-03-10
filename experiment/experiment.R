rm(list=ls())
set.seed(10)
n <- 100
partition <- c(0, 0.5, 1)
dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), round(partition*n))
lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, partition)
min_gamma = 0.01
max_gamma = 1000
max_iter = 10
delta = 10
verbose = T
X = dat$X
y = dat$y

##########

zz <- oracle_tune_gamma_range(X, y, lambda, partition, delta, min_gamma,
                        max_gamma, max_iter, F)
