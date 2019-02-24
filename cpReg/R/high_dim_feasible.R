high_dim_feasible_estimate <- function(X, y, lambda, tau, M = 100,
                                       delta = 10, verbose = F){
  data <- list(X = X, y = y)
  wbs(data, data_length_func = function(x){nrow(x$X)},
      compute_cusum_func = .compute_regression_cusum,
      tau = tau, M = M, delta = delta, verbose = verbose,
      lambda = lambda)
}

#########

.compute_regression_cusum <- function(data, start, end, breakpoint, lambda){
  X1 <- data$X[(start+1):breakpoint,,drop = F]
  y1 <- data$y[(start+1):breakpoint]
  X2 <- data$X[(breakpoint+1):end,,drop = F]
  y2 <- data$y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, lambda)
  beta2 <- .lasso_regression(X2, y2, lambda)

  .l2norm(sqrt((breakpoint - start)*(end-breakpoint)/(end-start))*(beta1 - beta2))
}

.lasso_regression <- function(X, y, lambda){
  glmnet_res <- glmnet::glmnet(X, y, intercept = F)
  as.numeric(glmnet::coef.glmnet(glmnet_res, s = lambda))[-1]
}

