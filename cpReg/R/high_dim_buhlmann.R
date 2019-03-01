high_dim_buhlmann_estimate <- function(X, y, lambda, gamma,
                                       delta = 10, verbose = F){
  data <- list(X = X, y = y)
  partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
                   compute_cusum_func = .compute_regression_buhlmann,
                   tau_function = .buhlmann_threshold,
                   M = 0, delta = delta, verbose = verbose,
                   lambda = lambda, gamma = gamma)

  n <- nrow(X)
  list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
}

oracle_tune_gamma <- function(X, y, lambda, partition, factor = 3/4){
  coef_list <- .refit_high_dim(X, y, lambda, partition)
  n <- nrow(X)
  partition_idx <- round(partition*n)
  k <- length(coef_list)

  res <- min(sapply(1:(k-1), function(x){
    idx1 <- (partition_idx[x]+1):partition_idx[x+1]
    idx2 <- (partition_idx[x+1]+1):partition_idx[x+2]

    fit <- glmnet::glmnet(X[c(idx1, idx2),,drop = F], y[c(idx1, idx2)],
                          intercept = F)
    alternative_coef <- glmnet::coef.glmnet(fit, s = lambda*sqrt(length(c(idx1,idx2))))[-1]
    loss <- as.numeric(.l2norm(X[c(idx1, idx2),,drop=F]%*%alternative_coef - y[c(idx1, idx2)])^2)/n

    loss1 <- as.numeric(.l2norm(X[idx1,,drop=F]%*%coef_list[[x]] - y[idx1])^2)/n
    loss2 <- as.numeric(.l2norm(X[idx2,,drop=F]%*%coef_list[[x+1]] - y[idx2])^2)/n

    loss - (loss1+loss2)
  }))

  res*factor
}

##################

.compute_regression_buhlmann <- function(data, start, end, breakpoint, lambda, gamma){
  X1 <- data$X[(start+1):breakpoint,,drop = F]
  y1 <- data$y[(start+1):breakpoint]
  X2 <- data$X[(breakpoint+1):end,,drop = F]
  y2 <- data$y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, lambda*sqrt(breakpoint-start))
  beta2 <- .lasso_regression(X2, y2, lambda*sqrt(end-breakpoint))

  -1*(as.numeric(.l2norm(X1%*%beta1 - y1)^2)/nrow(data$X) +
    as.numeric(.l2norm(X2%*%beta2 - y2)^2)/nrow(data$X) + 2*gamma)
}

.buhlmann_threshold <- function(data, interval, lambda, gamma){
  stopifnot(interval[1] < interval[2])
  X <- data$X[(interval[1]+1):interval[2],,drop = F]
  y <- data$y[(interval[1]+1):interval[2]]

  len <- interval[2] - interval[1]
  beta <- .lasso_regression(X, y, lambda*sqrt(len))
  -1*(as.numeric(.l2norm(X%*%beta - y)^2)/nrow(data$X) + gamma)
}
