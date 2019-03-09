rm(list=ls())
library(simulation)
library(cpReg)

paramMat <- as.matrix(expand.grid(round(exp(seq(log(100), log(1000), length.out = 10))), c(1,2),
                                  1/2))
colnames(paramMat) <- c("n", "X_type", "d/n")

X_type_vec <- c("identity", "toeplitz", "equicorrelation")
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

#################

set.seed(10)
vec = paramMat[10,]
dat = rule(vec)

#################

true_beta <- create_coef(vec, full = T)

lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                              factor = 1/2)
gamma <- cpReg::oracle_tune_gamma(dat$X, dat$y, lambda, true_partition,
                                  factor = 3/8)
grouplambda <- cpReg::oracle_tune_grouplambda(dat$X, dat$y, true_partition)

screeningtau <- cpReg::oracle_tune_screeningtau(dat$X, dat$y, lambda, true_partition)
group_screeningtau <-  cpReg::oracle_tune_group_screeningtau(dat$X, dat$y, true_partition)

maxl2 <- max(apply(true_beta, 1, cpReg:::.l2norm))
delta <- max(round(vec["n"]/10), 10)

res4 <- cpReg::high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, 4, delta = round(vec["n"]/10))
res3 <- cpReg::high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, 3, delta = round(vec["n"]/12))
res2 <- cpReg::high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, 2, delta = round(vec["n"]/20))

c(res2$obj_val, res3$obj_val, res4$obj_val)
