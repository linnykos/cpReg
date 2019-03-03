high_dim_feasible_estimate <- function(X, y, lambda, tau, verbose = F){

}

##################

.high_dim_feasible_subroutine <- function(X, y, lambda, maxl2, partition){
  stopifnot(partition[1] == 0, partition[length(partition)] == 1)

  n <- nrow(X); d <- ncol(X); k <- length(partition)-1
  partition_idx <- round(partition*n)
  X_new <- .reformat_covariates(X, partition_idx)
  group_vec <- rep(1:d, times = k)
  fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                        group.multiplier = rep(1, d))
  coef_vec <- grpreg:::coef.grpreg(fit, lambda = lambda)[-1]

  # form beta matrix
  beta_mat <- matrix(0, nrow = n, ncol = d)
  beta_list <- vector("list", k)
  for(i in 1:k){
    idx <- (partition[i]+1):partition[i+1]
    beta_mat[idx,] <- coef_vec[(d*(i-1)+1):(d*i)]/sqrt(length(idx))
    beta_list[[i]] <- beta_mat[idx,]
  }

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

.compute_lambda1se <- function(lambda, cv_error, cv_sd){
  idx <- which.min(cv_error)
  max_error <- cv_error[idx] + cv_sd[idx]
  max(lambda[cv_error < max_error])
}


