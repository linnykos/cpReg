rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

paramMat <- as.matrix(expand.grid(round(exp(seq(log(20), log(1000), length.out = 10))), c(1,2),
                                  5))
colnames(paramMat) <- c("n", "X_type", "d")

X_type_vec <- c("identity", "toeplitz", "equicorrelation")
true_partition <- c(0,0.3,0.7,1)
#############

create_coef <- function(vec, full = F){
  beta1 <- c(rep(1, 2), rep(0, vec["d"]-2))
  beta2 <- c(rep(0, vec["d"]-2), rep(1, 2))
  lis <- list(beta1 = beta1, beta2 = beta2)

  if(!full){
    lis
  } else {
    mat <- matrix(0, nrow = vec["n"], ncol = vec["d"])
    idx <- round(true_partition*vec["n"])
    for(i in 1:(length(idx)-1)){
      mat[(idx[i]+1):idx[i+1],] <- rep(lis[[(i %% 2)+1]], each = idx[i+1]-idx[i])
    }
    mat
  }
}

rule <- function(vec){
  lis <- create_coef(vec, full = F)

  create_data(list(lis$beta1, lis$beta2, lis$beta1), round(true_partition*vec["n"]),
              cov_type = X_type_vec[vec["X_type"]])
}


############

set.seed(1)
vec <- paramMat[19,]
dat <- rule(vec)
omega <- min(eigen(stats::cov(dat$X))$value)
res1 <- low_dim_estimate(dat$X, dat$y, gamma = 5*log(vec["n"])/(vec["n"]),
                                 delta = max(min(10, vec["n"]/10),3), verbose = F)
