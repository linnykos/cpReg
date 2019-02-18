rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

paramMat <- as.matrix(expand.grid(round(exp(seq(log(20), log(1000), length.out = 10))), c(1,2),
                                  5))
colnames(paramMat) <- c("n", "X_type", "d")

X_type_vec <- c("identity", "toeplitz", "equicorrelation")
true_partition <- c(0,0.3,0.7,1)

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
      zz <- i %% 2; if(zz == 0) zz <- 2
      mat[(idx[i]+1):idx[i+1],] <- rep(lis[[zz]], each = idx[i+1]-idx[i])
    }
    mat
  }
}

vec = paramMat[4,]
y = 2
set.seed(y)

rule <- function(vec){
  lis <- create_coef(vec, full = F)

  create_data(list(lis$beta1, lis$beta2, lis$beta1), round(true_partition*vec["n"]),
              cov_type = X_type_vec[vec["X_type"]])
}

dat = rule(vec)
set.seed(y)
res1 <- unravel(low_dim_estimate(dat$X, dat$y, gamma = 5*log(vec["n"])/(vec["n"]),
                                 delta = max(min(10, vec["n"]/10),3), verbose = F))
set.seed(y)
# partition1 <- screening(res1, 1, verbose = F)

#########
tau = 1
M = 100
verbose = T
fit = res1

# initialize
n <- nrow(fit)
random_intervals <- .generate_intervals(n, M)
q <- dequer::queue()
dequer::pushback(q, c(0,n))
b_vec <- numeric(0)
counter <- 0

# iterate over all intervals recursively
while(length(q) > 0){
  if(verbose) print(counter)
  interval <- dequer::pop(q)
  interval_list <- .truncate(interval, random_intervals)
  res_list <- lapply(interval_list, .find_breakpoint, fit = fit)
  res <- res_list[[which.max(sapply(res_list, function(x){x$val}))]]

  # if passes threshold, recurse
  if(res$val >= tau){
    b_vec <- c(b_vec, res$b)
    dequer::pushback(q, c(interval[1], res$b))
    dequer::pushback(q, c(res$b+1, interval[2]))
  }

  counter <- counter+1
}


