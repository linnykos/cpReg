#' High dimensional estimator (feasible)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param K numeric
#' @param tau numeric
#' @param M numeric
#' @param delta numeric
#' @param max_candidates numeric
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_feasible_estimate <- function(X, y, lambda, K = NA, tau = NA,
                                       M = 0,
                                       delta = 10, max_candidates = NA,
                                       verbose = F){
  stopifnot(!is.na(K) | !is.na(tau))
  stopifnot(is.na(K) | is.na(tau))
  data <- list(X = X, y = y)

  if(!is.na(K)){
    partition <- cp_fixedstep(data, data_length_func = function(x){nrow(x$X)},
                              compute_cusum_func = .compute_regression_cusum,
                              K = K, delta = delta, max_candidates = max_candidates,
                              verbose = verbose, lambda = lambda)

  } else {
    stopifnot(!is.na(tau))
    tau_function <- function(data, interval, ...){
      tau
    }

    partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
                     compute_cusum_func = .compute_regression_cusum,
                     tau_function = tau_function, M = M, delta = delta,
                     max_candidates = max_candidates,
                     verbose = verbose, lambda = lambda)
  }

  n <- nrow(X)
  list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
}

#########

.compute_regression_cusum <- function(data, start, end, breakpoint, lambda){
  n <- nrow(data$X)
  X1 <- data$X[(start+1):breakpoint,,drop = F]
  y1 <- data$y[(start+1):breakpoint]
  X2 <- data$X[(breakpoint+1):end,,drop = F]
  y2 <- data$y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, .cp_to_glmnet(lambda, breakpoint-start))
  beta2 <- .lasso_regression(X2, y2, .cp_to_glmnet(lambda, end-breakpoint))

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
    glmnet::coef.glmnet(fit, s = .cp_to_glmnet(lambda, partition_idx[x+1]-partition_idx[x]))[-1]
  })
}

