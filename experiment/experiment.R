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

#under estimate
res1 <- changepoint_dp(obj$dat, lambda = lambda, gamma = 2,
                      skip_interval = 500, verbose = T)

#just right
res2 <- changepoint_dp(obj$dat, lambda = 10, gamma = 40,
                      skip_interval = 500, verbose = T)

#too many estimates
res3 <- changepoint_dp(obj$dat, lambda = lambda, gamma = 0.0001,
                      skip_interval = 500, verbose = T)

############################
lambda <- 10
changepoint_idx <- c(0,obj$partition,TT)
k <- length(changepoint_idx)
obj_middle <- sum(sapply(1:(k-1), function(i){
  len <- changepoint_idx[i+1]-changepoint_idx[i]

  stationary_ar(obj$dat[(changepoint_idx[i]+1):changepoint_idx[i+1],],
                lambda = lambda*sqrt(len), intercept = F)$obj_val
}))

changepoint_idx <- c(0,750,1500,2250,3000,4250,5000)
k <- length(changepoint_idx)
obj_many <- sum(sapply(1:(k-1), function(i){
  len <- changepoint_idx[i+1]-changepoint_idx[i]

  stationary_ar(obj$dat[(changepoint_idx[i]+1):changepoint_idx[i+1],],
                lambda = lambda*sqrt(len), intercept = F)$obj_val
}))

changepoint_idx <- c(0,1500,5000)
k <- length(changepoint_idx)
obj_few <- sum(sapply(1:(k-1), function(i){
  len <- changepoint_idx[i+1]-changepoint_idx[i]

  stationary_ar(obj$dat[(changepoint_idx[i]+1):changepoint_idx[i+1],],
                lambda = lambda*sqrt(len), intercept = F)$obj_val
}))

plot(c(obj_few, obj_middle, obj_many))
