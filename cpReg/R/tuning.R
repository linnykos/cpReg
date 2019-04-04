#' Tuning via cross validation
#'
#' @param method function to fit
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param K_range vector of number of changepoints to consider
#' @param cv_verbose boolean
#' @param max_iter numeric
#' @param ... additional parameters for \code{method}
#'
#' @return list
#' @export
tuning_cross_validation <- function(method, X, y, K_range = c(1:5),
                                    cv_verbose = F, max_iter = 5,
                                    ...){
  dat <- list(X = X, y = y)
  paramMat <- .tuning_lambda_K_pairing(method, dat, K_range = K_range,
                                       cv_verbose = cv_verbose, max_iter = max_iter,
                                       ...)

  cv_val <- sapply(1:nrow(paramMat), function(x){
    if(cv_verbose) print(paste0("On K=", paramMat[x,2]))
    .cross_validate(method = method, dat, K = paramMat[x,2],
                    lambda = paramMat[x,1], cv_verbose = cv_verbose, ...)
  })

  list(K = paramMat[which.min(cv_val), "K"],
       lambda = paramMat[which.min(cv_val), "lambda"],
       paramMat = cbind(paramMat, cv_val))
}

########################

.tuning_lambda_K_pairing <- function(method, dat, K_range = c(1:5),
                                    cv_verbose = F, max_iter = 5,
                                    ...){
  n <- nrow(dat$X)

  # first fit the entire dataset via lasso
  fit <- glmnet::cv.glmnet(dat$X, dat$y, intercept = F, grouped = F)
  new_lambda_vec <- rep(.glmnet_to_cp(fit$lambda.min, n), length(K_range))
  lambda_vec <- rep(Inf, length(K_range))
  iter <- 1

  # iterate between lambda and finding partition
  while(iter < max_iter & sum(abs(new_lambda_vec - lambda_vec)) > 1e-3){
    if(cv_verbose) print(paste0("Lambda values: ",
                                paste0(round(new_lambda_vec, 2), collapse = ", ")))

    lambda_vec <- new_lambda_vec

    # try fitting the method now for different values of K
    res_list <- lapply(1:length(K_range), function(i){
      method(dat$X, dat$y, K= K_range[i], lambda = lambda_vec[i], ...)
    })

    # based on the estimated partitions, refit lambda
    new_lambda_vec <- sapply(1:length(res_list), function(i){
      oracle_tune_lambda(dat$X, dat$y, res_list[[i]]$partition/n, lambda.min = T)
    })

    iter <- iter + 1
  }

  mat <- cbind(lambda_vec, K_range)
  colnames(mat) <- c("lambda", "K")
  mat
}

.cross_validate <- function(method, dat, K, lambda,
                            nfolds = 5, cv_verbose = F, ...){
  n <- nrow(dat$X)
  fold_id <- rep(1:nfolds, times = ceiling(n/nfolds))[1:n]

  if(cv_verbose) print("Fitting each model")
  res_list <- lapply(1:nfolds, function(i){
    idx <- which(fold_id != i)
    method(dat$X[idx,,drop = F], dat$y[idx], K = K, lambda = lambda, ...)
  })

  if(cv_verbose) print("Evaluating the error")
  stats::median(sapply(1:nfolds, function(i){
    .out_of_sample_prediction(dat, fold_id, fold = i, res_list[[i]])
  }))
}

.out_of_sample_prediction <- function(dat, fold_id, fold, res){
  stopifnot(length(res$partition) > 2)

  n <- nrow(dat$X)
  idx <- which(fold_id == fold)
  cp_idx <- .convert_cp_idx(res$partition, fold_id, fold)
  stopifnot(cp_idx[1] == 0 & cp_idx[length(cp_idx)] == n)
  cp_idx_wide <- sort(unique(unlist(lapply(cp_idx[2:(length(cp_idx)-1)], function(x){
    pmax(0, pmin(n, x+c(-1,0,1)))
  }))))

  # handle cases away from the margin
  idx_away <- idx[which(!idx %in% cp_idx_wide)]
  if(length(idx_away) > 0){
    res_away <- sapply(idx_away, function(x){
      i <- max(which(cp_idx < x))
      .l2norm(dat$y[x] - dat$X[x,,drop=F]%*%res$coef_list[[i]])^2
    })
  } else {
    res_away <- numeric(0)
  }

  idx_close <- idx[which(idx %in% cp_idx_wide)]
  if(length(idx_close) > 0){
    res_close <- sapply(idx_close, function(x){
      i <- which.min(abs(x - cp_idx))
      stopifnot(i != 1 & i != length(cp_idx))

      coef_avg <- apply(cbind(res$coef_list[[i-1]], res$coef_list[[i]]), 1, mean)
      .l2norm(dat$y[x] - dat$X[x,,drop=F]%*%res$coef_list[[i]])^2
    })
  } else {
    res_close <- numeric(0)
  }

  stats::median(c(res_away, res_close))
}

.convert_cp_idx <- function(partition, fold_id, fold){
  stopifnot(length(fold_id) >= max(partition))
  tmp <- rep(NA, length(fold_id))
  tmp[which(fold_id != fold)] <- 1:length(fold_id[fold_id != fold])

  c(0, sapply(2:length(partition), function(x){
    which(tmp == partition[x])
  }))
}
