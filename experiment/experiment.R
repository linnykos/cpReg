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

res <- changepoint_dp(obj$dat, lambda = lambda, gamma = 5,
                      min_spacing = 500, skip_interval = 500, verbose = T)

obj_val <- res$obj_val - (length(res$partition)-1)*5*nrow(obj$dat)
# thought experiment:
# so does the objective function above (for res, w/o penalty) match
#  the one i get when i run gamma?
#  let's try this: I want to decompose the obj function into (nll + lambda) and
#   (gamma). run the above function. it should be exact, so i'll get something
#   for (nll+lambda) and (gamma). Then, use the gamma function and see
#   what the (nll+lambda) is when I add a lot of extraneous spacing


obj3 <- sum(sapply(1:(length(res$partition)-1), function(i){
  len <- res$partition[i+1]-res$partition[i]

  stationary_ar(obj$dat[(res$partition[i]+1):res$partition[i+1],],
                lambda = lambda*sqrt(len), intercept = F)$obj_val
}))

#something is deeply wrong..
