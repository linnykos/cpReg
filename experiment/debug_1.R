# set up the data
rm(list=ls())
set.seed(10)
M <- 5; TT <- 50
nu <- rep(0.1, M)
A <- 0.3*diag(M)
dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = 5)

thres_u <- round(stats::quantile(dat[dat > 0], probs = 0.75))
lambda <- 0.25
gamma <- 10
verbose = T

TT <- nrow(dat); M <- ncol(dat)
h <- vector("list", TT)
min_spacing <- 10

for(i in 2:(TT-1)){
  print(i)
  if(verbose & i %% floor(TT/10) == 0) cat('*')

  obj_vec <- rep(NA, i-1)
  lis <- vector("list", length = i-1)

  # search
  for(j in 1:(i-1)){
    if(i - j <= min_spacing){
      lis[[j]] <- list(A = NA)
      obj_vec[j] <- Inf
    } else {
      lis[[j]] <- stationary_ar(dat = dat[j:i, , drop = F], thres_u = thres_u,
                                lambda = lambda * sqrt(i-j), verbose = F,
                                intercept = F)
      tmp <- h[[j]]$obj_val
      if(is.null(tmp)) tmp <- 0
      obj_vec[j] <- lis[[j]]$obj_val + tmp + gamma
    }
  }

  # select
  idx <- which.min(obj_vec)

  # finalize
  h_prev <-  h[[idx]]
  if(all(is.null(h_prev))){
    partition <- c(1, i); A_list <- list(lis[[idx]]$A)
  } else {
    partition <- c(h_prev$partition, i)
    A_list <- h_prev$A_list; A_list[[length(A_list)+1]] <- lis[[idx]]$A
  }

  h_obj <- .construct_list_obj(min(obj_vec), partition = partition,
                               A_list = A_list)
  h[[i]] <- h_obj
}
