rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

paramMat <- as.matrix(expand.grid(round(exp(seq(log(100), log(1000), length.out = 10))), c(4),
                                  1/2))
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

criterion <- function(dat, vec, y){
  true_beta <- create_coef(vec, full = T)

  n <- nrow(dat$X); p <- ncol(dat$X)

  delta <- max(round(n/10), 10)

  # parameter for feasible (and others)
  feasible_paramMat <- cpReg::tuning_cross_validation(high_dim_feasible_estimate, X = dat$X,
                                                      y = dat$y, K_range = c(1:5),
                                                      max_iter = 10, cv_verbose = F)

  # parameter for buhlmann
  buhlmann_paramMat <- cpReg::tuning_cross_validation(high_dim_buhlmann_estimate, X = dat$X,
                                                      y = dat$y, K_range = c(1:5),
                                                      max_iter = 10, cv_verbose = F)


  res1 <- cpReg::high_dim_feasible_estimate(dat$X, dat$y, lambda = feasible_paramMat$lambda,
                                            K = feasible_paramMat$K, tau = NA,
                                            verbose = F, delta = delta, M = 0)
  res2 <- cpReg::high_dim_buhlmann_estimate(dat$X, dat$y, lambda = buhlmann_paramMat$lambda,
                                            K = buhlmann_paramMat$K, gamma = NA,
                                            verbose = F, delta = delta)

  beta_mat1 <- cpReg::unravel(res1)
  beta_mat2 <- cpReg::unravel(res2)

  beta_error1 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat1[x,] - true_beta[x,])^2}))/n
  beta_error2 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat2[x,] - true_beta[x,])^2}))/n

  haus1 <- cpReg::hausdorff(res1$partition, round(true_partition*n))
  haus2 <- cpReg::hausdorff(res2$partition, round(true_partition*n))

  list(beta_error = list(beta_error1, beta_error2),
       haus = list(haus1, haus2),
       partition = list(res1$partition, res2$partition),
       parameters = list(feasible_paramMat = feasible_paramMat,
                         buhlmann_paramMat = buhlmann_paramMat))
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[10,]), paramMat[10,], 1)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = 20,
                                        cores = 15, as_list = T,
                                        filepath = "../results/high_dim_simulation_tmp.RData",
                                        verbose = T)
save.image("../results/high_dim_simulation_cv.RData")
