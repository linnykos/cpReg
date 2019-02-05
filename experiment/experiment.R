rm(list=ls())
set.seed(10)
dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
X = dat$X
y = dat$y
gamma = 1
delta = 10
verbose = T

#############

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
  if(verbose & idx_i %% floor(length(n_seq)/10) == 0) cat('*')
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
