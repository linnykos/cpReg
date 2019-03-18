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

###########

set.seed(10)
vec = paramMat[1,]
dat = rule(vec)
X = dat$X
y = dat$y
K = 2
S = 20

lambda2 <- oracle_tune_lambda(dat$X, dat$y, true_partition)
tau2 <- oracle_tune_screeningtau(X, y, lambda2, true_partition)

p <- ncol(X); n <- nrow(X); X2 <- t(X)
lambda <- 5*sqrt(K*log(p*n))
gamma <- 5*sqrt(n*S*log(p))

beta <- CVXR::Variable(n,p)
D <- matrix(0, nrow=n, ncol=n)
diag(D) <- 1; diag(D[,-1]) <- -1
D <- D[-n,]

prob <- CVXR::Problem(CVXR::Minimize((CVXR::cvxr_norm(y-CVXR::diag(beta %*% X2),p=2)^2)/n   +
                                       lambda*CVXR::sum_entries( CVXR::cvxr_norm(beta, axis = 2, p=2)) +
                                       gamma*CVXR::cvxr_norm(D%*%beta,p=1)))

result <- CVXR::psolve(prob, solver="SCS")
zz <- result$getValue(beta)
plot(zz[,1])

yy <- screening(zz, tau = tau2, M = 0)

xx <- unravel(list(partition = yy, coef_list = .refit_high_dim(X, y, lambda2, yy/n)))
plot(xx[,1])
