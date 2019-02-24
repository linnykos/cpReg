high_dim_feasible_estimate <- function(X, y, lambda, tau, M = 100,
                                       delta = 10, verbose = F){
  data <- list(X = X, y = y)
  wbs(data, data_length_func = function(x){nrow(x$X)},
      compute_cusum_func = .compute_regression_cusum,
      tau = tau, M = M, delta = delta, verbose = verbose,
      lambda = lambda)
}

oracle_tune_lambda <- function(X, y, partition){
  stopifnot(partition[1] == 0, partition[length(partition)] == 1)

  n <- nrow(X)
  k <- length(partition)-1
  partition_idx <- round(partition*n)

  median(sapply(1:k, function(x){
    fit <- glmnet::cv.glmnet(X[(partition_idx[x]+1):partition_idx[x+1],,drop = F],
                      y[(partition_idx[x]+1):partition_idx[x+1]],
                      intercept = F, grouped = F)
    fit$lambda.1se
  }))
}

oracle_tune_tau <- function(X, y, partition, lambda, factor = 3/4){
  stopifnot(partition[1] == 0, partition[length(partition)] == 1)

  n <- nrow(X)
  k <- length(partition)-1
  partition_idx <- round(partition*n)

  coef_list <- lapply(1:k, function(x){
    fit <- glmnet::glmnet(X[(partition_idx[x]+1):partition_idx[x+1],,drop = F],
                             y[(partition_idx[x]+1):partition_idx[x+1]],
                             intercept = F)
    glmnet::coef.glmnet(fit, s = lambda)[-1]
  })

  res <- min(sapply(1:(k-1), function(x){
    s <- partition_idx[x]
    b <- partition_idx[x+1]
    e <- partition_idx[x+2]
    .l2norm(sqrt((b - s)*(e-b)/(e-s))*(coef_list[[x]] - coef_list[[x+1]]))
  }))

  res*factor
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

