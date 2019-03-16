#' High dimensional estimator - Buhlmann
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param K numeric
#' @param delta numeric
#' @param max_candidates numeric
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_buhlmann_estimate <- function(X, y, lambda, K,
                                       delta = 10, max_candidates = NA,
                                       verbose = F){
  stopifnot(nrow(X) == length(y))
  data <- list(X = X, y = y)

  partition <- cp_fixedstep(data, data_length_func = function(x){nrow(x$X)},
                        compute_cusum_func = .compute_regression_buhlmann,
                        K = K, delta = delta, max_candidates = max_candidates,
                        verbose = verbose, lambda = lambda)

  n <- nrow(X)
  list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
}

##################

.compute_regression_buhlmann <- function(data, start, end, breakpoint, lambda){
  n <- nrow(data$X)
  X1 <- data$X[(start+1):breakpoint,,drop = F]
  y1 <- data$y[(start+1):breakpoint]
  X2 <- data$X[(breakpoint+1):end,,drop = F]
  y2 <- data$y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, lambda/sqrt(breakpoint-start))
  beta2 <- .lasso_regression(X2, y2, lambda/sqrt(end-breakpoint))

  -1*(as.numeric(.l2norm(X1%*%beta1 - y1)^2) +
        as.numeric(.l2norm(X2%*%beta2 - y2)^2))
}
