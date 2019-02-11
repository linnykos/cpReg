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

criterion <- function(dat, vec, y){
  omega <- min(eigen(stats::cov(dat$X))$value)
  res1 <- unravel(low_dim_estimate(dat$X, dat$y, gamma = 5*log(vec["n"])/(vec["n"]*omega),
                                   delta = max(min(10, vec["n"]/10),3), verbose = F))
  res2 <- SGL_solver(dat$X, dat$y)

  true_beta <- create_coef(vec, full = T)

  beta_error1 <- sum(sapply(1:vec["n"], function(x){.l2norm(res1[x,] - true_beta[x,])^2}))/vec["n"]
  beta_error2 <- sum(sapply(1:vec["n"], function(x){.l2norm(res2[x,] - true_beta[x,])^2}))/vec["n"]

  partition1 <- screening(res1, 1, verbose = F)
  partition2 <- screening(res2, 1, verbose = F)

  haus1 <- hausdorff(partition1, round(true_partition*vec["n"]))
  haus2 <- hausdorff(partition2, round(true_partition*vec["n"]))

  list(beta_error1 = beta_error1, beta_error2 = beta_error2,
       haus1 = haus1, haus2 = haus2)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = paramMat[,"trials"],
                                        cores = 15, as_list = F,
                                        filepath = "main_powercurve_onejump_fl_tmp.RData",
                                        verbose = T)
save.image("main_powercurve_onejump_fl.RData")
