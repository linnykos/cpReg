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

set.seed(1)
vec <- paramMat[13,]
dat <- rule(vec)

true_beta <- create_coef(vec, full = T)

lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
tau <- cpReg::oracle_tune_tau(dat$X, dat$y, lambda, true_partition,
                              factor = 3/4)
gamma <- cpReg::oracle_tune_gamma(dat$X, dat$y, lambda, true_partition,
                                  factor = 1/2)
grouplambda <- cpReg::oracle_tune_grouplambda(dat$X, dat$y, true_partition)

screeningtau <- cpReg::oracle_tune_screeningtau(dat$X, dat$y, lambda, true_partition)
group_screeningtau <-  cpReg::oracle_tune_group_screeningtau(dat$X, dat$y, true_partition)

maxl2 <- max(apply(true_beta, 1, cpReg:::.l2norm))
K <- 4
delta <- max(round(vec["n"]/10), 10)

# res3 <- cpReg::high_dim_infeasible_estimate(dat$X, dat$y, grouplambda, maxl2, K, delta = delta)
###############
X = dat$X
y = dat$y
n <- nrow(X)
combn_mat <- .enumerate_possibilites(n, K, delta)
# res_list <- lapply(1:ncol(combn_mat), function(x){
#   print(x)
#   .high_dim_infeasible_subroutine(X, y, lambda, maxl2, combn_mat[,x])
# })

#.high_dim_infeasible_subroutine(X, y, lambda, maxl2, combn_mat[,1])

################
partition = combn_mat[,1]
stopifnot(partition[1] == 0, partition[length(partition)] == nrow(X))

n <- nrow(X); d <- ncol(X); k <- length(partition)-1
X_new <- .reformat_covariates(X, partition)
group_vec <- rep(1:d, times = k)
fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                      group.multiplier = rep(1, d), lamba.min = lambda)


# do something if the lambda asked is too small
if(min(fit$lambda) > lambda){
  lambda_seq <- seq(log(max(fit$lambda)), log(lambda/2), length.out = 100)
  fit <- grpreg::grpreg(X_new, y, group = group_vec, penalty = "grLasso",
                        lambda = lambda_seq,
                        group.multiplier = rep(1, d))
}

coef_vec <- as.numeric(grpreg:::coef.grpreg(fit, lambda = lambda)[-1])

beta_list <- .convert_grouplasso_to_cp(coef_vec, partition, n, d, as_list = T)
beta_mat <- .convert_grouplasso_to_cp(coef_vec, partition, n, d)

# check that the l2norm constraint is met
for(i in 1:d){
  val <- .l2norm(beta_mat[,i])
  if(val > maxl2){
    beta_mat[,i] <- beta_mat[,i]/val*maxl2
  }
}
