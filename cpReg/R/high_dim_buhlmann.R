#' High dimensional estimator - Buhlmann
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param K numeric
#' @param gamma numeric
#' @param delta numeric
#' @param max_candidates numeric
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_buhlmann_estimate <- function(X, y, lambda, K = NA, gamma = NA,
                                       delta = 10, max_candidates = NA,
                                       verbose = F){
  stopifnot(nrow(X) == length(y))
  stopifnot(is.na(gamma) | is.na(K))
  stopifnot(!is.na(gamma) | !is.na(K))
  data <- list(X = X, y = y)

  if(!is.na(gamma)){
    partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
                     compute_cusum_func = .compute_regression_buhlmann,
                     tau_function = .buhlmann_threshold,
                     M = 0, delta = delta,
                     verbose = verbose, lambda = lambda, gamma = gamma,
                     with_penalty = T)
  } else {
    partition <- cp_fixedstep(data, data_length_func = function(x){nrow(x$X)},
                              compute_cusum_func = .compute_regression_buhlmann,
                              K = K, delta = delta, max_candidates = max_candidates,
                              verbose = verbose, lambda = lambda)
  }

  n <- nrow(X)
  list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
}

##################

.compute_regression_buhlmann <- function(data, start, end, breakpoint, lambda,
                                         with_penalty = F, gamma = 0, ...){
  n <- nrow(data$X)
  X1 <- data$X[(start+1):breakpoint,,drop = F]
  y1 <- data$y[(start+1):breakpoint]
  X2 <- data$X[(breakpoint+1):end,,drop = F]
  y2 <- data$y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, .cp_to_glmnet(lambda, breakpoint-start))
  beta2 <- .lasso_regression(X2, y2, .cp_to_glmnet(lambda, end-breakpoint))

  if(with_penalty) penalty <- 2*gamma*n else penalty <- 0

  -1*(as.numeric(.l2norm(X1%*%beta1 - y1)^2) +
        as.numeric(.l2norm(X2%*%beta2 - y2)^2) + penalty)
}


.buhlmann_threshold <- function(data, interval, lambda, gamma, ...){
  n <- nrow(data$X)
  stopifnot(interval[1] < interval[2])
  X <- data$X[(interval[1]+1):interval[2],,drop = F]
  y <- data$y[(interval[1]+1):interval[2]]

  len <- interval[2] - interval[1]
  beta <- .lasso_regression(X, y, .cp_to_glmnet(lambda, len))
  -1*(as.numeric(.l2norm(X%*%beta - y)^2) + lambda*sqrt(len)*sum(abs(beta)) + gamma*n)
}
