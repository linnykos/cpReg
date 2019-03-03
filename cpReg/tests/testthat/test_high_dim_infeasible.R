# context("Test high dimensional infeasible")
#
# ## .reformat_covariates is correct
#
# test_that(".reformat_covariates works", {
#   X <- matrix(1:30, 6)
#   partition <- c(0,2,5,6)
#   res <- .reformat_covariates(X, partition)
#
#   expect_true(is.matrix(res))
#   expect_true(!is.vector(res))
#   expect_true(all(dim(res) == c(nrow(X), ncol(X)*(length(partition)-1))))
# })
#
# test_that(".reformat_covariates gives the correct matrix", {
#   set.seed(10)
#   X <- matrix(rnorm(30), 6)
#   partition <- c(0,2,5,6)
#   res <- .reformat_covariates(X, partition)
#
#   # construct the right matrix at each stage
#   n <- nrow(X); d <- ncol(X); k <- length(partition)-1
#   for(i in 1:k){
#     res2 <- matrix(0, n, d)
#     res2[(partition[i]+1):n,] <- X[(partition[i]+1):n,]
#     expect_true(sum(abs(res2 - res[,(d*(i-1)+1):(d*i)])) < 1e-6)
#   }
# })
