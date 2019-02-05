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

##########

changepoint_idx <- c(0,obj$partition,TT)
k <- length(changepoint_idx)

i = 3
len <- changepoint_idx[i+1]-changepoint_idx[i]

stationary_ar(obj$dat[(changepoint_idx[i]+1):changepoint_idx[i+1],],
              lambda = 0.00002*sqrt(len), intercept = F)
