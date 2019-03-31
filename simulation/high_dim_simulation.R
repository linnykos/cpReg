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
  lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
  tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                                factor = 1/2)

  # parameter for buhlmann
  gamma_range <- cpReg::oracle_tune_gamma_range(dat$X, dat$y, lambda, K = 3, verbose = F)
  if(any(is.na(gamma_range$gamma))) gamma <- gamma_range$min_gamma else gamma <- mean(gamma_range$gamma)

  # parameter for infeasible
  grouplambda <- cpReg::oracle_tune_grouplambda(dat$X, dat$y, true_partition)
  maxl2 <- max(apply(true_beta, 1, cpReg:::.l2norm))
  K <- 2

  # parameter for GLL1
  screening_tau <- oracle_tune_screeningtau(dat$X, dat$y, lambda, true_partition)
  S <- 20
  lambda2 <- 5*sqrt(K*log(p*n))
  gamma2 <- 5*sqrt(n*S*log(p))

  res1 <- cpReg::high_dim_feasible_estimate(dat$X, dat$y, lambda = lambda, K = 2, tau = NA,
                                    verbose = F, max_candidates = NA, delta = delta, M = 50)
  res2 <- cpReg::high_dim_buhlmann_estimate(dat$X, dat$y, lambda = lambda, gamma = gamma, K = NA,
                                            verbose = F, max_candidates = NA, delta = delta)
  res3 <- cpReg::high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, K = K, delta = delta)
  res4 <- GLL1_solver(dat$X, dat$y, lambda = lambda2, gamma = gamma2)

  beta_mat1 <- cpReg::unravel(res1)
  beta_mat2 <- cpReg::unravel(res2)
  beta_mat3 <- cpReg::unravel(res3)

  beta_error1 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat1[x,] - true_beta[x,])^2}))/n
  beta_error2 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat2[x,] - true_beta[x,])^2}))/n
  beta_error3 <- sum(sapply(1:vec["n"], function(x){cpReg:::.l2norm(beta_mat3[x,] - true_beta[x,])^2}))/n
  beta_error4 <- sum(sapply(1:n, function(x){cpReg:::.l2norm(res4[x,] - true_beta[x,])^2}))/n

  haus1 <- cpReg::hausdorff(res1$partition, round(true_partition*n))
  haus2 <- cpReg::hausdorff(res2$partition, round(true_partition*n))
  haus3 <- cpReg::hausdorff(res3$partition, round(true_partition*n))
  partition4b <- cpReg::screening(res4, tau = screening_tau, M = 0)
  haus4b <- cpReg::hausdorff(partition4b, round(true_partition*n))

  beta_mat4b <- cpReg::unravel(list(partition = partition4b,
                                    coef_list = cpReg:::.refit_high_dim(dat$X, dat$y, lambda, partition4b/n)))
  beta_error4b <- sum(sapply(1:n, function(x){cpReg:::.l2norm(beta_mat4b[x,] - true_beta[x,])^2}))/n


  list(beta_error = list(beta_error1, beta_error2, beta_error3, beta_error4,
                         beta_error4b),
       haus = list(haus1, haus2, haus3, haus4b),
       partition = list(res1$partition, res2$partition, res3$partition,
                        partition4b),
       parameters = list(lambda = lambda, tau = tau,
                         gamma_range = gamma_range, gamma = gamma,
                         grouplambda = grouplambda,
                         screening_tau = screening_tau))
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[10,]), paramMat[10,], 1)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = 50,
                                        cores = 15, as_list = T,
                                        filepath = "../results/high_dim_simulation_tmp.RData",
                                        verbose = T)
save.image("../results/high_dim_simulation_block.RData")
