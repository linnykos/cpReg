#' High dimesional infeasible estimator
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param maxl2 numeric
#' @param K numeric
#' @param delta numeric (minimal spacing between the points to be brute searched over)
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_infeasible_estimate <- function(X, y, lambda, maxl2, K,
                                         delta = 10, verbose = F){
  n <- nrow(X)
  combn_mat <- .enumerate_possibilites(n, K, delta)
  res_list <- lapply(1:ncol(combn_mat), function(x){
    .high_dim_infeasible_subroutine(X, y, lambda, maxl2, combn_mat[,x])
  })

  idx <- which.min(sapply(res_list, function(x){x$obj_val}))

  list(partition = as.numeric(combn_mat[,idx]), coef_list = res_list[[idx]]$coef_list,
       obj_val = res_list[[idx]]$obj_val)
}

##################

.enumerate_possibilites <- function(n, K, delta = 10){
  stopifnot(delta < n)
  seq_vec <- seq(delta, n-delta, by = delta)
  rbind(0, utils::combn(seq_vec, K), n)
}

.high_dim_infeasible_subroutine <- function(X, y, lambda, maxl2, partition){
  stopifnot(partition[1] == 0, partition[length(partition)] == nrow(X))

  n <- nrow(X); d <- ncol(X); k <- length(partition)-1
  X_new <- .reformat_covariates(X, partition)
  group_vec <- rep(1:d, times = k)
  # warning: grpreg is not able to turn off the intercept
  fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                        group.multiplier = rep(1, d))

  # do something if the lambda asked is too small
  if(min(fit$lambda) > lambda){
    lambda_seq <- seq(log(max(fit$lambda)), log(lambda/2), length.out = 100)
    fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                          lambda = lambda_seq,
                          group.multiplier = rep(1, d))
  }

  coef_vec <- as.numeric(stats::coef(fit, lambda = lambda)[-1])

  beta_list <- .convert_grouplasso_to_cp(coef_vec, partition, n, d, as_list = T)
  beta_mat <- .convert_grouplasso_to_cp(coef_vec, partition, n, d)

  # check that the l2norm constraint is met
  for(i in 1:d){
    val <- .l2norm(beta_mat[,i])
    if(val > maxl2){
      beta_mat[,i] <- beta_mat[,i]/val*maxl2
    }
  }

  # evaluate objective function
  y_pred <- sapply(1:n, function(i){X[i,]%*%beta_mat[i,]})
  obj_val <- .l2norm(y - y_pred)^2/n + lambda*sum(apply(beta_mat, 2, .l2norm))

  list(obj_val = obj_val, coef_list = beta_list)
}

.convert_grouplasso_to_cp <- function(coef_vec, partition, n, d, as_list = F){

  k <- length(partition)-1

  if(as_list){
    beta_list <- vector("list", k)
  } else {
    beta_mat <- matrix(0, nrow = n, ncol = d)
  }

  for(i in 1:k){
    idx <- (partition[i]+1):partition[i+1]
    tmp <- coef_vec[(d*(i-1)+1):(d*i)]/sqrt(length(idx))
    if(as_list){
      beta_list[[i]] <- tmp
    } else {
      beta_mat[idx,] <- t(sapply(1:length(idx), function(x){tmp}))
    }
  }

  if(as_list){beta_list} else{beta_mat}
}

.reformat_covariates <- function(X, partition){
  stopifnot(is.matrix(X))
  stopifnot(partition[1] == 0, partition[length(partition)] == nrow(X),
            all(sort(partition) == partition), length(unique(partition)) == length(partition))

  n <- nrow(X); d <- ncol(X); k <- length(partition)-1
  X_new <- matrix(0, n, d*k)
  for(i in 1:k){
    len <- partition[i+1]-partition[i]
    X_new[(partition[i]+1):partition[i+1], (d*(i-1)+1):(d*i)] <- X[(partition[i]+1):partition[i+1],,drop = F]/sqrt(len)
  }

  X_new
}
