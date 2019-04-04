rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

paramMat <- as.matrix(expand.grid(round(exp(seq(log(100), log(1000), length.out = 10))),
                                  c(1), 1/2))
colnames(paramMat) <- c("n", "X_type", "d/n")

X_type_vec <- c("identity", "toeplitz", "equicorrelation", "block")
true_partition <- c(0,0.3,0.7,1)

#############

create_coef <- function(vec, full = F){
  d <- 50 # d <- vec["d/n"]*vec["n"]
  beta1 <- c(rep(1, 10), rep(0, d-10))
  beta2 <- c(rep(0, d-10), rep(1, 10))
  lis <- list(beta1 = beta1, beta2 = beta2)

  if(!full){
    lis
  } else {
    mat <- matrix(0, nrow = vec["n"], ncol = d)
    idx <- round(true_partition*vec["n"])
    for(i in 1:(length(idx)-1)){
      zz <- i %% 2; if(zz == 0) zz <- 2
      mat[(idx[i]+1):idx[i+1],] <- rep(lis[[zz]], each = idx[i+1]-idx[i])
    }
    mat
  }
}

rule <- function(vec){
  lis <- create_coef(vec, full = F)

  cpReg::create_data(list(lis$beta1, lis$beta2, lis$beta1), round(true_partition*vec["n"]),
                     cov_type = X_type_vec[vec["X_type"]])
}

set.seed(1)
vec = paramMat[10,]
dat = rule(vec)

##################

set.seed(1)
paramMat <- .tuning_lambda_K_pairing(method = high_dim_buhlmann_estimate,
                                     dat, cv_verbose = T)

# lambda = paramMat[2,1]
# K = paramMat[2,2]
# nfolds = 5
# method = high_dim_buhlmann_estimate
#
# n <- nrow(dat$X)
# margin <- ceiling(n/nfolds)
# fold_id <- rep(1:nfolds, times = margin)[1:n]
#
# res_list <- lapply(1:nfolds, function(i){
#   idx <- which(fold_id != i)
#   method(dat$X[idx,,drop = F], dat$y[idx], K = K, lambda = lambda)
# })
#
# .out_of_sample_prediction(dat, fold_id, fold = 1, res_list[[1]])
#

# .cross_validate(method = high_dim_buhlmann_estimate, dat, K = paramMat[2,2],
#                 lambda = paramMat[2,1])


cv_val <- sapply(1:nrow(paramMat), function(x){
  if(cv_verbose) print(paste0("On K=", paramMat[x,2]))
  .cross_validate(method = high_dim_buhlmann_estimate, dat, K = paramMat[x,2],
                  lambda = paramMat[x,1], cv_verbose = T)
})
