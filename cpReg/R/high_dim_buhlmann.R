#' High dimensional estimator - Buhlmann
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param gamma numeric
#' @param delta numeric
#' @param max_candidates numeric
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_buhlmann_estimate <- function(X, y, lambda, gamma,
                                       delta = 10, max_candidates = NA,
                                       verbose = F){
  data <- list(X = X, y = y)
  partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
                   compute_cusum_func = .compute_regression_buhlmann,
                   tau_function = .buhlmann_threshold,
                   M = 0, delta = delta, max_candidates = 10,
                   verbose = verbose,
                   lambda = lambda, gamma = gamma)

  n <- nrow(X)
  list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
}

#' Tune gamma (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param partition vector with values between 0 and 1
#' @param factor numeric
#'
#' @return numeric
#' @export
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

#' Tune gamma range (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param k numeric
#' @param delta numeric
#' @param min_gamma numeric
#' @param max_gamma numeric
#' @param max_iter numeric
#' @param verbose boolean
#'
#' @return list
#' @export
oracle_tune_gamma_range <- function(X, y, lambda, k, delta = 10, min_gamma = 0.01,
                                    max_gamma = 1000, max_iter = 10, verbose = T){

  # find a suitable lower and upper startpoint
  min_gamma <- .initial_gamma_overshoot(X, y, lambda, k, min_gamma, delta = delta, smaller = T,
                                        verbose = verbose)
  max_gamma <- .initial_gamma_overshoot(X, y, lambda, k, max_gamma, delta = delta, smaller = F,
                                        verbose = verbose)

  # find a suitable starting point
  res <- .initialize_gamma_binarysearch(X, y, lambda, k = k, delta = delta,
                                        min_gamma = min_gamma, max_gamma = max_gamma, verbose = verbose)

  min_gamma_vec <- res$min_gamma_vec; max_gamma_vec <- res$max_gamma_vec; gamma_vec <- res$gamma
  if(length(gamma_vec) == 0) return(list(gamma = NA, min_gamma = max(min_gamma_vec), max_gamma = min(max_gamma_vec)))
  iter <- 1

  # throttle on both sides
  while(iter <= max_iter){
    try_min <- mean(c(max(min_gamma_vec), min(gamma_vec)))
    res <-  high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = try_min,
                                           delta = delta)
    if(length(res$partition) - 1 > k){
      min_gamma_vec <- c(min_gamma_vec, try_min)
    } else if(length(res$partition) - 1 == k){
      gamma_vec <- c(gamma_vec, try_min)
    } else stop()

    try_max <- mean(c(min(max_gamma_vec), max(gamma_vec)))
    res <-  high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = try_max,
                                           delta = delta)

    if(length(res$partition) - 1 < k){
      max_gamma_vec <- c(max_gamma_vec, try_max)
    } else if(length(res$partition) - 1 == k){
      gamma_vec <- c(gamma_vec, try_max)
    } else stop()

    iter <- iter + 1
    if(verbose) print(list(min_gamma_vec = min_gamma_vec, max_gamma_vec = max_gamma_vec, gamma_vec = gamma_vec))
  }

  stopifnot(max(min_gamma_vec) < min(gamma_vec), min(max_gamma_vec) > max(gamma_vec))

  list(gamma = range(gamma_vec), min_gamma = max(min_gamma_vec), max_gamma = min(max_gamma_vec))
}

#' Tune screening tau for high dimensional infeasible (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param partition vector with values between 0 and 1
#' @param factor numeric
#'
#' @return numeric
#' @export
oracle_tune_screeningtau <- function(X, y, lambda, partition, factor = 1/4){
  stopifnot(all(partition >= 0))

  n <- nrow(X)
  coef_list <- .refit_high_dim(X, y, lambda, partition)
  partition_idx <- round(partition*n)
  fit <- list(partition = partition_idx, coef_list = coef_list)
  mat <- unravel(fit)

  vec <- sapply(2:(length(partition_idx)-1), function(x){
    .compute_cusum(mat, partition_idx[x-1], partition_idx[x+1], partition_idx[x])
  })

  min(vec)*factor
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

.initial_gamma_overshoot <- function(X, y, lambda, k, gamma, delta = 10, smaller = T,
                                     verbose = T, max_iter = 30){
  iter <- 1

  while(iter <= max_iter){
    if(verbose) print(paste0("Still trying to initialize ", ifelse(smaller, "minimum", "maximum"), " gamma"))
    res_low <- high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = gamma,
                                          delta = delta)
    if(smaller){
      if(length(res_low$partition)-1 <= k) gamma <- gamma/2 else break()
    }

    if(!smaller){
      if(length(res_low$partition)-1 >= k) gamma <- gamma*2 else break()
    }

    iter <- iter+1
  }

  gamma
}

.initialize_gamma_binarysearch <- function(X, y, lambda, k, delta = 10, min_gamma = 0.01,
                                           max_gamma = 1000, tol = 1e-3, verbose = T,
                                           max_iter = 30){
  min_gamma_vec <- min_gamma; max_gamma_vec <- max_gamma; gamma <- numeric(0)
  iter <- 1

  # phase 1: find ONE suitable gamma
  while(abs(max(min_gamma_vec) - min(max_gamma_vec)) > tol & iter < max_iter){
    try_gamma <- mean(c(max(min_gamma_vec), min(max_gamma_vec)))
    res <- high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = try_gamma,
                                          delta = delta)

    if(verbose) print(paste0("Trying ", try_gamma, " resulting in ", length(res$partition)-1))

    if(length(res$partition)-1 == k) {
      gamma <- try_gamma; break()
    } else if(length(res$partition)-1 > k){
      min_gamma_vec <- c(try_gamma, min_gamma_vec)
    } else {
      max_gamma_vec <- c(try_gamma, max_gamma_vec)
    }

    iter <- iter + 1
  }

  list(gamma = gamma, min_gamma_vec = min_gamma_vec, max_gamma_vec = max_gamma_vec)
}
