rm(list=ls())
library(simulation)
library(cpReg)

paramMat <- as.matrix(expand.grid(round(exp(seq(log(100), log(300), length.out = 10))), c(1,2),
                                  1/2))
colnames(paramMat) <- c("n", "X_type", "d/n")

X_type_vec <- c("identity", "toeplitz", "equicorrelation")
true_partition <- c(0,0.3,0.7,1)

#############

create_coef <- function(vec, full = F){
  d <- vec["d/n"]*vec["n"]
  beta1 <- c(rep(1, 2), rep(0, d-2))
  beta2 <- c(rep(0, d-2), rep(1, 2))
  lis <- list(beta1 = beta1, beta2 = beta2)

  if(!full){
    lis
  } else {
    mat <- matrix(0, nrow = vec["n"], ncol = vec["d/n"]*vec["n"])
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

criterion <- function(dat, vec, y){
  grouplambda <- cpReg::oracle_tune_grouplambda(dat$X, dat$y, true_partition)
  true_beta <- create_coef(vec, full = T)
  maxl2 <- max(apply(true_beta, 1, cpReg:::.l2norm))
  K <- 2
  delta <- vec["n"]/10

  res3 <- high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, K, delta)

  beta_mat3 <- cpReg::unravel(res3)

  beta_error3 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat3[x,] - true_beta[x,])^2}))/vec["n"]

  haus3 <- cpReg::hausdorff(res3$partition, round(true_partition*vec["n"]))

  list(beta_error3 = beta_error3, haus3 = haus3,
       partition3 = res3$partition,
       grouplambda = grouplambda)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(2); criterion(rule(paramMat[4,]), paramMat[4,], 2)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = 20,
                                        cores = 15, as_list = T,
                                        filepath = "../results/high_dim_simulation_tmp.RData",
                                        verbose = T)
save.image("../results/high_dim_simulation.RData")
