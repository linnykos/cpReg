#' Tune lambda (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param partition vector with values between 0 and 1
#' @param lambda.min boolean on whether or not to use the \code{lambda.min} or \code{lambda.1se} in the cross validation
#'
#' @return numeric
#' @export
oracle_tune_lambda <- function(X, y, partition, lambda.min = F){
  stopifnot(partition[1] == 0, partition[length(partition)] == 1)

  n <- nrow(X)
  k <- length(partition)-1
  partition_idx <- round(partition*n)

  stats::median(sapply(1:k, function(x){
    fit <- glmnet::cv.glmnet(X[(partition_idx[x]+1):partition_idx[x+1],,drop = F],
                             y[(partition_idx[x]+1):partition_idx[x+1]],
                             intercept = F, grouped = F)
    val <- ifelse(lambda.min, fit$lambda.min, fit$lambda.1se)
     .glmnet_to_cp(val, partition_idx[x+1]-partition_idx[x])
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


#' Tune group lambda (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param partition vector with values between 0 and 1
#' @param return_refit boolean
#'
#' @return numeric
#' @export
oracle_tune_grouplambda <- function(X, y, partition, return_refit = F){
  n <- nrow(X); d <- ncol(X); k <- length(partition)-1
  partition_idx <- round(partition*n)

  X_new <- .reformat_covariates(X, partition_idx)
  group_vec <- rep(1:d, times = k)

  fit <- grpreg::cv.grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                           group.multiplier = rep(1, d))

  if(!return_refit){
    .compute_lambda1se(fit$lambda, fit$cve, fit$cvse)
  } else {
    lambda <- .compute_lambda1se(fit$lambda, fit$cve, fit$cvse)
    fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                          group.multiplier = rep(1, d))
    coef_vec <- as.numeric(stats::coef(fit, lambda = lambda)[-1])

    .convert_grouplasso_to_cp(coef_vec, partition_idx, n, d)
  }
}

#' Tune screening tau for high dimensional infeasible (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param partition vector with values between 0 and 1
#' @param factor numeric
#'
#' @return numeric
#' @export
oracle_tune_group_screeningtau <- function(X, y, partition, factor = 1/4){
  stopifnot(all(partition >= 0))

  n <- nrow(X)
  mat <- oracle_tune_grouplambda(X, y, partition, return_refit = T)

  partition_idx <- round(partition*n)
  vec <- sapply(2:(length(partition_idx)-1), function(x){
    .compute_cusum(mat, partition_idx[x-1], partition_idx[x+1], partition_idx[x])
  })

  min(vec)*factor
}


#' Tune gamma range (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param K the number of desired segments (number of changepoints + 1)
#' @param delta numeric
#' @param min_gamma numeric
#' @param max_gamma numeric
#' @param max_iter numeric
#' @param verbose boolean
#'
#' @return list
#' @export
oracle_tune_gamma_range <- function(X, y, lambda, K, delta = 10, min_gamma = 0.01,
                                    max_gamma = 16000, max_iter = 10, verbose = T){

  # find a suitable lower and upper startpoint
  min_gamma <- .initial_gamma_overshoot(X, y, lambda, K, min_gamma, delta = delta, smaller = T,
                                        verbose = verbose)
  max_gamma <- .initial_gamma_overshoot(X, y, lambda, K, max_gamma, delta = delta, smaller = F,
                                        verbose = verbose)

  # find a suitable starting point
  res <- .initialize_gamma_binarysearch(X, y, lambda, K = K, delta = delta,
                                        min_gamma = min_gamma, max_gamma = max_gamma, verbose = verbose)

  min_gamma_vec <- res$min_gamma_vec; max_gamma_vec <- res$max_gamma_vec; gamma_vec <- res$gamma
  if(length(gamma_vec) == 0) return(list(gamma = NA, min_gamma = max(min_gamma_vec), max_gamma = min(max_gamma_vec)))
  iter <- 1

  # throttle on both sides
  while(iter <= max_iter){
    try_min <- mean(c(max(min_gamma_vec), min(gamma_vec)))
    res <-  high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = try_min,
                                       delta = delta)
    if(length(res$partition) - 1 > K){
      min_gamma_vec <- c(min_gamma_vec, try_min)
    } else if(length(res$partition) - 1 == K){
      gamma_vec <- c(gamma_vec, try_min)
    } else stop()

    try_max <- mean(c(min(max_gamma_vec), max(gamma_vec)))
    res <-  high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = try_max,
                                       delta = delta)

    if(length(res$partition) - 1 < K){
      max_gamma_vec <- c(max_gamma_vec, try_max)
    } else if(length(res$partition) - 1 == K){
      gamma_vec <- c(gamma_vec, try_max)
    } else stop()

    iter <- iter + 1
    if(verbose) print(list(min_gamma_vec = min_gamma_vec, max_gamma_vec = max_gamma_vec, gamma_vec = gamma_vec))
  }

  stopifnot(max(min_gamma_vec) < min(gamma_vec), min(max_gamma_vec) > max(gamma_vec))

  list(gamma = range(gamma_vec), min_gamma = max(min_gamma_vec), max_gamma = min(max_gamma_vec))
}

#' Tune gamma to minimize hausdorff distance (oracle)
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param partition vector with values between 0 and 1
#' @param k_vec numeric vector of the number of desired segments (number of changepoints + 1)
#' @param delta numeric
#' @param min_gamma numeric
#' @param max_gamma numeric
#' @param verbose boolean
#'
#' @return numeric
#' @export
oracle_tune_gamma_hausdorff <- function(X, y, lambda, partition,
                                        k_vec,
                                        delta = 10, min_gamma = 0.01,
                                        max_gamma = 1000,
                                        verbose = F){
  n <- nrow(X)
  res <- lapply(k_vec, function(i){
    if(verbose) print(i)
    oracle_tune_gamma_range(X, y, lambda, i, delta = delta, min_gamma = min_gamma,
                            max_gamma = max_gamma, verbose = F)
  })

  quality_vec <- sapply(1:length(res), function(i){
    if(verbose) print(i)
    if(any(is.na(res[[i]]$gamma))) return(NA)

    gamma <- mean(res[[i]]$gamma)
    tmp <- high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = gamma,
                                      delta = delta)$partition
    hausdorff(tmp, round(partition*n))
  })

  if(all(is.na(quality_vec))){
    return(max(sapply(res, function(i){i$min_gamma}), na.rm = T))
  }

  idx <- which.min(quality_vec)
  mean(res[[idx]]$gamma)
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

#############################


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

.initialize_gamma_binarysearch <- function(X, y, lambda, K, delta = 10, min_gamma = 0.01,
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

    if(length(res$partition)-1 == K) {
      gamma <- try_gamma; break()
    } else if(length(res$partition)-1 > K){
      min_gamma_vec <- c(try_gamma, min_gamma_vec)
    } else {
      max_gamma_vec <- c(try_gamma, max_gamma_vec)
    }

    iter <- iter + 1
  }

  list(gamma = gamma, min_gamma_vec = min_gamma_vec, max_gamma_vec = max_gamma_vec)
}



.compute_lambda1se <- function(lambda, cv_error, cv_sd){
  idx <- which.min(cv_error)
  max_error <- cv_error[idx] + cv_sd[idx]
  max(lambda[cv_error < max_error])
}



