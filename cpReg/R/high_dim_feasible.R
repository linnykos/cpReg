#' High dimensional estimator (feasible)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param tau numeric
#' @param M numeric
#' @param delta numeric
#' @param max_candidates numeric
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_feasible_estimate <- function(X, y, lambda, tau, M = 100,
                                       delta = 10, max_candidates = NA,
                                       verbose = F){
  tau_function <- function(data, interval, ...){
    tau
  }

  data <- list(X = X, y = y)
  partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
      compute_cusum_func = .compute_regression_cusum,
      tau_function = tau_function, M = M, delta = delta, max_candidates = max_candidates,
      verbose = verbose, lambda = lambda)

  n <- nrow(X)
  list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
}

#' Tune lambda (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param partition vector with values between 0 and 1
#'
#' @return numeric
#' @export
oracle_tune_lambda <- function(X, y, partition){
  stopifnot(partition[1] == 0, partition[length(partition)] == 1)

  n <- nrow(X)
  k <- length(partition)-1
  partition_idx <- round(partition*n)

  stats::median(sapply(1:k, function(x){
    fit <- glmnet::cv.glmnet(X[(partition_idx[x]+1):partition_idx[x+1],,drop = F],
                      y[(partition_idx[x]+1):partition_idx[x+1]],
                      intercept = F, grouped = F)
    fit$lambda.min/sqrt(partition_idx[x+1]-partition_idx[x])
  }))
}

#' Tune tau (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param partition vector with values between 0 and 1
#' @param factor numeric
#'
#' @return numeric
#' @export
oracle_tune_tau <- function(X, y, lambda, partition, factor = 3/4){
  coef_list <- .refit_high_dim(X, y, lambda, partition)
  n <- nrow(X)
  partition_idx <- round(partition*n)
  k <- length(coef_list)

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

  beta1 <- .lasso_regression(X1, y1, lambda*sqrt(breakpoint-start))
  beta2 <- .lasso_regression(X2, y2, lambda*sqrt(end-breakpoint))

  .l2norm(sqrt((breakpoint - start)*(end-breakpoint)/(end-start))*(beta1 - beta2))
}

.lasso_regression <- function(X, y, lambda){
  glmnet_res <- glmnet::glmnet(X, y, intercept = F)
  as.numeric(glmnet::coef.glmnet(glmnet_res, s = lambda))[-1]
}

.refit_high_dim <- function(X, y, lambda, partition){
  stopifnot(partition[1] == 0, partition[length(partition)] == 1)

  n <- nrow(X)
  k <- length(partition)-1
  partition_idx <- round(partition*n)

  lapply(1:k, function(x){
    fit <- glmnet::glmnet(X[(partition_idx[x]+1):partition_idx[x+1],,drop = F],
                          y[(partition_idx[x]+1):partition_idx[x+1]],
                          intercept = F)
    glmnet::coef.glmnet(fit, s = lambda*sqrt(partition_idx[x+1]-partition_idx[x]))[-1]
  })
}

