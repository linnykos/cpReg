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
#tau <- oracle_tune_tau(dat$X, dat$y, lambda, true_partition)
gamma <- oracle_tune_gamma(dat$X, dat$y, lambda, true_partition)

#res <- high_dim_feasible_estimate(dat$X, dat$y, lambda = lambda, tau = tau,
#                                  verbose = T)

res <- high_dim_buhlmann_estimate(dat$X, dat$y, lambda = lambda, gamma = gamma,
                                  verbose = F)

##############

X <- dat$X; y <- dat$y; partition = true_partition; factor = 3/4

coef_list <- .refit_high_dim(X, y, lambda, partition)
n <- nrow(X)
partition_idx <- round(partition*n)
k <- length(coef_list)

print(partition_idx)

res <- sapply(1:(k-1), function(x){
  idx1 <- (partition_idx[x]+1):partition_idx[x+1]
  idx2 <- (partition_idx[x+1]+1):partition_idx[x+2]

  fit <- glmnet::glmnet(X[c(idx1, idx2),,drop = F], y[c(idx1, idx2)],
                        intercept = F)
  alternative_coef <- glmnet::coef.glmnet(fit, s = lambda*sqrt(length(c(idx1,idx2))))[-1]
  loss <- as.numeric(.l2norm(X[c(idx1, idx2),,drop=F]%*%alternative_coef - y[c(idx1, idx2)])^2)/n

  loss1 <- as.numeric(.l2norm(X[idx1,,drop=F]%*%coef_list[[x]] - y[idx1])^2)/n
  loss2 <- as.numeric(.l2norm(X[idx2,,drop=F]%*%coef_list[[x+1]] - y[idx2])^2)/n

  loss - (loss1+loss2)
})

print(res)
