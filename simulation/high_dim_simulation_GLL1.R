rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

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
  lambda2 <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
  tau2 <- oracle_tune_screeningtau(dat$X, dat$y, lambda2, true_partition)

  K <- 2
  S <- 20

  n <- nrow(dat$X); p <- ncol(dat$X)

  lambda <- 5*sqrt(K*log(p*n))
  gamma <- 5*sqrt(n*S*log(p))

  res4 <- GLL1_solver(dat$X, dat$y, lambda = lambda, gamma = gamma)

  beta_error4 <- sum(sapply(1:n, function(x){cpReg:::.l2norm(res4[x,] - true_beta[x,])^2}))/n

  partition4b <- cpReg::screening(res4, tau = tau2, M = 0)
  haus4b <- cpReg::hausdorff(partition4b, round(true_partition*n))

  beta_mat4b <- cpReg::unravel(list(partition = partition4b,
                             coef_list = .refit_high_dim(dat$X, dat$y, lambda2, partition4b/n)))
  beta_error4b <- sum(sapply(1:n, function(x){cpReg:::.l2norm(beta_mat4b[x,] - true_beta[x,])^2}))/n

  list(beta_error = list(beta_error4, beta_error4b),
       haus = list(haus4b),
       partition = list(partition4b))
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
