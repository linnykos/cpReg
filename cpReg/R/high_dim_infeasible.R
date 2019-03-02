# high_dim_feasible_estimate <- function(X, y, lambda, tau, verbose = F){
#
# }
#
# ##################
#
# .high_dim_feasible_subroutine <- function(X, y, lambda, maxl2, partition){
#   stopifnot(partition[1] == 0, partition[length(partition)] == 1)
#
#   n <- nrow(X); d <- ncol(X); k <- length(partition)-1
#   partition_idx <- round(partition*n)
#   X_new <- .reformat_covariates(X, partition_idx)
#   group_vec <- rep(1:d, times = k)
#   fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso")
#   coef_vec <- grpreg:::coef.grpreg(fit, lambda = lambda)[-1]
#
#   # check that the l2norm constraint is met
#   for(i in 1:k){
#     idx <- which(group_vec == i)
#     val <- .l2norm(coef_vec[idx])
#     if(val > maxl2){
#       coef_vec[idx] <- coef_vec[idx]/.l2norm(coef_vec[idx])*maxl2
#     }
#   }
#
#   # evaluate objective function
#   .l2norm(y - X_new %*% coef_vec)^2/n + lambda*
#
#   list(partition = partition, coef_list = .refit_high_dim(X, y, lambda, partition/n))
# }
#
# .reformat_covariates <- function(X, partition){
#   stopifnot(is.matrix(X))
#   stopifnot(partition[1] == 0, partition[length(partition)] == nrow(X),
#             all(sort(partition) == partition), length(unique(partition)) == length(partition))
#
#   n <- nrow(X); d <- ncol(X); k <- length(partition)-1
#   X_new <- matrix(0, n, d*k)
#   for(i in 1:k){
#     X_new[(partition[i]+1):n, (d*(i-1)+1):(d*i)] <- X[(partition[i]+1):n,,drop = F]
#   }
#
#   X_new
# }
#
# .compute_lambda1se <- function(lambda, cv_error, cv_sd){
#   idx <- which.min(cv_error)
#   max_error <- cv_error[idx] + cv_sd[idx]
#   max(lambda[cv_error < max_error])
# }
#
#
