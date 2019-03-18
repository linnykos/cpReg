#' Screening
#'
#' @param fit object from estimator
#' @param tau numeric
#' @param M numeric
#' @param delta numeric
#' @param verbose boolean
#'
#' @return numeric
#' @export
screening <- function(fit, tau, M = 100, delta = 1, verbose = F){
  tau_function <- function(data, interval, ...){
    tau
  }

  wbs(fit, data_length_func = nrow,
      compute_cusum_func = .compute_cusum,
      tau_function = tau_function, M = M, delta = delta, verbose = verbose)
}

#' Hausdorff distance
#'
#' @param set1 numeric
#' @param set2 numeric
#' @param one.sided boolean
#'
#' @return numeric
#' @export
hausdorff <- function(set1, set2, one.sided = FALSE){
  #handle corner cases
  if(length(set1) == 0 | length(set2) == 0) return(NA)
  if(length(set2) == 1) set2 = c(set2, set2)
  if(length(set1) == 1) set1 = c(set1, set1)

  dist.mat = sapply(set1, function(i){abs(i-set2)})
  if(class(dist.mat) != "matrix") dist.mat = as.matrix(dist.mat, nrow = 1)
  dist.vecx = apply(dist.mat, 2, min)

  if(!one.sided) dist.vecy = apply(dist.mat, 1, min) else dist.vecy = 0

  max(dist.vecx, dist.vecy)
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

################

.compute_cusum <- function(fit, s, e, b){
  stopifnot(s < b, b < e, s >= 0, e <= nrow(fit))

  factor1 <- sqrt((e-b)/((e-s)*(b-s)))
  factor2 <- sqrt((b-s)/((e-s)*(e-b)))
  .l2norm(factor1*colSums(fit[(s+1):b,,drop = F]) - factor2*colSums(fit[(b+1):e,,drop = F]))
}

.l2norm <- function(x){
  sqrt(sum(x^2))
}

.create_full_matrix <- function(coef_list, partition){
  n <- partition[length(partition)]
  mat <- matrix(0, nrow = n, ncol = length(coef_list[[1]]))
  len <- length(partition)

  for(i in 1:(len-1)){
    mat[(partition[i]+1):partition[i+1],] <- rep(coef_list[[i]], each = partition[i+1]-partition[i])
  }

  mat
}
