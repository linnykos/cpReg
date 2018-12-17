rm(list=ls())
set.seed(10)
M <- 5; TT <- 5000
nu <- rep(0, M)
A_list <- list(0.75*diag(M),
               matrix(0, M, M),
               0.5*matrix(runif(M*M), M, M))
changepoint_perc <- c(0.3, 0.6)

obj <- generative_model_cp(nu, A_list, changepoint_perc, timesteps = TT)
lambda <- .lambda_oracle(obj, intercept = F)

###########

# res <- changepoint_dp(obj$dat, lambda = lambda, gamma = 5,
#                       min_spacing = 500, skip_interval = 500, verbose = T)

dat <- obj$dat
thres_u = round(stats::quantile(dat[dat > 0], probs = 0.75))
gamma = 5
min_spacing = 500
skip_interval = 500
verbose = T

TT <- nrow(dat); M <- ncol(dat)
TT_seq <- seq(1, TT, by = skip_interval)
TT_seq[length(TT_seq)] <- TT
h <- vector("list", length(TT_seq)-1)

for(i in TT_seq){
  if(i == 1) next()
  idxi <- which(TT_seq == i)
  if(verbose & idxi %% floor(length(TT_seq)/10) == 0) cat('*')

  obj_vec <- rep(NA, idxi-1)
  lis <- vector("list", length = length(TT_seq[TT_seq < i]))

  # search
  for(j in TT_seq[TT_seq < i]){
    idxj <- which(TT_seq == j)
    if(i - j < min_spacing){
      lis[[idxj]] <- list(A = NA)
      obj_vec[idxj] <- Inf
    } else {
      lis[[idxj]] <- stationary_ar(dat = dat[j:i, , drop = F], thres_u = thres_u,
                                   lambda = lambda * sqrt(i-j), verbose = F,
                                   intercept = F)
      tmp <- h[[idxj]]$obj_val
      if(is.null(tmp)) tmp <- 0
      obj_vec[idxj] <- lis[[idxj]]$obj_val + tmp + gamma*TT
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
  h[[idxi]] <- h_obj
}

TT_seq2 <- TT_seq; TT_seq2[1] <- 0
obj_vec_gamma <- sapply(1:(length(TT_seq2)-1), function(i){
  len <- TT_seq2[i+1]-TT_seq2[i]

  stationary_ar(obj$dat[(TT_seq2[i]+1):TT_seq2[i+1],], thres_u = thres_u,
                lambda = lambda*sqrt(len), intercept = F)$obj_val
})
sum(obj_vec_gamma)
sum(obj_vec_gamma)+(length(TT_seq2)-1)*5*TT

h[[10]]$obj_val-(length(h[[10]]$partition)-1)*5*TT

#################

# do some inspection
TT_seq
length(h)
# let's start with h[[1]], should be innocuous
h[[1]]$obj_val



stopifnot(abs(h[[1]]$obj_val - 5*TT - obj_vec_gamma[1]) <= 1e-6)

