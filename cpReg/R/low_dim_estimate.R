#' Low dimensional estimator
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param gamma numeric
#' @param delta numeric
#' @param verbose boolean
#'
#' @return list containing \code{obj_val}, \code{partition}, \code{coef_list}
#' @export
low_dim_estimate <- function(X, y, gamma, delta = 5, verbose = F){
  # create list of results to store DP
  n <- nrow(X); d <- ncol(X)
  gamma2 <- gamma*n
  n_seq <- .format_sequence(n, delta)
  h <- vector("list", length(n_seq))

  # special case for the first split
  partition <- NA; coef_list <- NA
  h[[1]] <- .construct_list_obj(obj_val = 0, partition = 0, coef_list = list(numeric(0)))

  # search through all possible split points
  for(i in n_seq[-1]){
    idx_i <- which(n_seq == i)
    if(verbose & idx_i %% max(floor(length(n_seq)/10),1) == 0) cat('*')
    n_subseq <- n_seq[n_seq < i]
    obj_vec <- rep(NA, length(n_subseq))
    lis <- vector("list", length = length(n_subseq))

    for(j in n_subseq){
      idx_j <- which(n_seq == j)
      lis[[idx_j]] <- .local_regression(X = X[(j+1):i, , drop = F], y = y[(j+1):i])
      obj_vec[idx_j] <- h[[idx_j]]$obj_val + lis[[idx_j]]$obj_val + gamma2
    }

    # select
    idx <- which.min(obj_vec)

    # finalize
    h_prev <- h[[idx]]
    partition <- .special_concatenate(h_prev$partition, i, F)
    coef_list <- .special_concatenate(h_prev$coef_list, lis[[idx]]$coef, T)

    h[[idx_i]]  <- .construct_list_obj(min(obj_vec), partition = partition,
                                 coef_list = coef_list)
  }

  h[[length(h)]]
}

#' Unravel estimate into matrix
#'
#' @param fit object from one of the estimation functions
#'
#' @return \code{n} by \code{d} matrix
#' @export
unravel <- function(fit){
  mat <- matrix(0, nrow = max(fit$partition), ncol = length(fit$coef_list[[1]]))

  idx <- fit$partition
  for(i in 1:(length(idx)-1)){
    mat[(idx[i]+1):idx[i+1],] <- rep(fit$coef_list[[i]], each = idx[i+1]-idx[i])
  }

  mat
}

################

.format_sequence <- function(n, delta){
  seq_vec <- seq(0, n, by = delta)
  seq_vec[length(seq_vec)] <- n

  stopifnot(length(seq_vec) >= 2)
  seq_vec
}

.construct_list_obj <- function(obj_val, partition, coef_list){
  stopifnot(length(coef_list) == length(partition)-1 | partition == 0)
  stopifnot(all(partition == sort(partition)),
            length(partition) == length(unique(partition)))

  list(obj_val = obj_val, partition = partition, coef_list = coef_list)
}

.local_regression <- function(X, y){
  res <- stats::lm(y ~ X - 1)
  list(obj_val = sum(res$residuals^2), coef = stats::coef(res))
}

.special_concatenate <- function(obj1, obj2, is_list = F){
  if(is_list){
    stopifnot(is.list(obj1))
    obj1[[length(obj1)+1]] <- obj2

    tmp <- sapply(obj1, length)
    obj1 <- obj1[tmp != 0]
    tmp <- tmp[tmp != 0]
    stopifnot(length(unique(tmp)) == 1)

  } else {
    stopifnot(!is.list(obj1))
    obj1 <- c(obj1, obj2)
  }

  obj1
}
