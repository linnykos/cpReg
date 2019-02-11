rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

paramMat <- as.matrix(expand.grid(round(exp(seq(log(20), log(1000), length.out = 10))), c(1,2),
                                  5))
colnames(paramMat) <- c("n", "X_type", "d")

X_type_vec <- c("identity", "toeplitz", "equicorrelation")

#############

create_coef <- function(vec, full = F){
  beta1 <- c(rep(1, 2), rep(0, vec["d"]-2))
  beta2 <- c(rep(0, vec["d"]-2), rep(1, 2))
  lis <- list(beta1 = beta1, beta2 = beta2)

  if(!full){
    lis
  } else {
    mat <- matrix(0, nrow = vec["n"], ncol = vec["d"])
    idx <- round(c(0,0.3,0.7,1)*vec["n"])
    for(i in 1:(length(idx)-1)){
      mat[(idx[i]+1):idx[i+1],] <- rep(lis[[(2 %% i)+1]], each = idx[i+1]-idx[i])
    }
    mat
  }
}

rule <- function(vec){
  lis <- create_coef(vec, full = F)

  create_data(list(lis$beta1, lis$beta2, lis$beta1), round(c(0,0.3,0.7,1)*vec["n"]),
              cov_type = X_type_vec[vec["X_type"]])
}

set.seed(10)
vec <- paramMat[1,]
dat <- rule(vec)
X <- dat$X; y <- dat$y; gamma = 1; delta = max(3, vec["n"]/10)

n <- nrow(X); d <- ncol(X)
gamma2 <- gamma*n
n_seq <- .format_sequence(n, delta)
h <- vector("list", length(n_seq))

# special case for the first split
partition <- NA; coef_list <- NA
h[[1]] <- .construct_list_obj(obj_val = 0, partition = 0, coef_list = list(numeric(0)))
i = n_seq[2]
idx_i <- which(n_seq == i)
verbose = T
if(verbose & idx_i %% max(floor(length(n_seq)/10),1) == 0) cat('*')
