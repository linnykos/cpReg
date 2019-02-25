#' Create data
#'
#' @param coef_list list of numerics
#' @param partition vector
#' @param cov_type character
#'
#' @return list containing \code{X} and code{y}
#' @export
create_data <- function(coef_list, partition, cov_type = "identity"){
  stopifnot(length(partition) == length(coef_list)+1, partition[1] == 0)
  stopifnot(length(unique(sapply(coef_list, length))) == 1)

  n <- partition[length(partition)]
  d <- length(coef_list[[1]])
  X <- matrix(0, nrow = n, ncol = d)
  class(X) <- c(cov_type, class(X)[length(class(X))])
  X <- .cov_generation(X)
  y <- numeric(n)

  for(i in 1:(length(partition)-1)){
    y[(partition[i]+1):partition[i+1]] <- X[(partition[i]+1):partition[i+1],] %*% coef_list[[i]] + stats::rnorm(partition[i+1] - partition[i])
  }

  list(X = X, y = y)
}

.cov_generation <- function(X, ...) {
  UseMethod(".cov_generation")
}

.cov_generation.identity <- function(X){
  n <- nrow(X); d <- ncol(X)
  matrix(stats::rnorm(n*d), nrow = n, ncol = d)
}

.cov_generation.toeplitz <- function(X){
  n <- nrow(X); d <- ncol(X)
  cov_mat <- stats::toeplitz(0.8^(0:(d-1)))
  MASS::mvrnorm(n, rep(0, d), Sigma = cov_mat)
}

.cov_generation.equicorrelation <- function(X){
  n <- nrow(X); d <- ncol(X)
  cov_mat <- matrix(0.2, d, d)
  diag(cov_mat) <- 1
  MASS::mvrnorm(n, rep(0, d), Sigma = cov_mat)
}
