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

criterion <- function(dat, vec, y){
  true_beta <- create_coef(vec, full = T)

  delta <- max(round(vec["n"]/10), 10)
  lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
  tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                                factor = 1/2)
  grouplambda <- cpReg::oracle_tune_grouplambda(dat$X, dat$y, true_partition)
  group_screeningtau <-  cpReg::oracle_tune_group_screeningtau(dat$X, dat$y, true_partition)

  maxl2 <- max(apply(true_beta, 1, cpReg:::.l2norm))
  K <- 4

  res1 <- cpReg::high_dim_feasible_estimate(dat$X, dat$y, lambda = lambda, tau = tau,
                                    verbose = F, max_candidates = NA, delta = delta, M = 50)
  res2 <- cpReg::high_dim_buhlmann_estimate(dat$X, dat$y, lambda = lambda, k = 2,
                                            verbose = F, max_candidates = NA, delta = delta)
  res3 <- cpReg::high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, K, delta = delta)

  beta_mat1 <- cpReg::unravel(res1)
  beta_mat2 <- cpReg::unravel(res2)
  beta_mat3 <- cpReg::unravel(res3)

  beta_error1 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat1[x,] - true_beta[x,])^2}))/vec["n"]
  beta_error2 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat2[x,] - true_beta[x,])^2}))/vec["n"]
  beta_error3 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat3[x,] - true_beta[x,])^2}))/vec["n"]

  haus1 <- cpReg::hausdorff(res1$partition, round(true_partition*vec["n"]))
  haus2 <- cpReg::hausdorff(res2$partition, round(true_partition*vec["n"]))
  haus3 <- cpReg::hausdorff(res3$partition, round(true_partition*vec["n"]))

  #### now see what changes after regression

  partition3b <- cpReg::screening(beta_mat3, group_screeningtau, M = 0)

  beta_mat3b <- cpReg::unravel(list(partition = partition3b,
                                    coef_list = cpReg:::.refit_high_dim(dat$X, dat$y, lambda, partition3b/vec["n"])))

  beta_error3b <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat3b[x,] - true_beta[x,])^2}))/vec["n"]

  haus3b <- cpReg::hausdorff(partition3b, round(true_partition*vec["n"]))

  list(beta_error = list(beta_error1, beta_error2, beta_error3,
                         beta_error3b),
       haus = list(haus1, haus2, haus3, haus3b),
       partition = list(res1$partition, res2$partition, res3$partition,
                        partition3b),
       parameters = list(lambda = lambda, tau = tau,
                         grouplambda = grouplambda,
                         group_screeningtau = group_screeningtau))
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[10,]), paramMat[10,], 1)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = 50,
                                        cores = 15, as_list = T,
                                        filepath = "../results/high_dim_simulation_tmp.RData",
                                        verbose = T)
save.image("../results/high_dim_simulation.RData")
