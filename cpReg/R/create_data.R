create_data <- function(coef_list, partition){
  stopifnot(length(partition) == length(coef_list)+1, partition[1] == 0)
  stopifnot(length(unique(sapply(coef_list, length))) == 1)

  n <- partition[length(partition)]
  d <- length(coef_list[[1]])
  X <- matrix(stats::rnorm(n*d), nrow = n, ncol = d)
  y <- numeric(n)

  for(i in 1:(length(partition)-1)){
    y[(partition[i]+1):partition[i+1]] <- X[(partition[i]+1):partition[i+1],] %*% coef_list[[i]] + stats::rnorm(partition[i+1] - partition[i])
  }

  list(X = X, y = y)
}
