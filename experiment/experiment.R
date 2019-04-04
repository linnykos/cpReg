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

#######

y <- 2
set.seed(y)
vec <- paramMat[2,]
dat <- rule(vec)

true_beta <- create_coef(vec, full = T)

n <- nrow(dat$X); p <- ncol(dat$X)

delta <- max(round(n/20), 5)

# parameter for feasible (and others)
set.seed(10)
# feasible_paramMat <- cpReg::tuning_cross_validation(cpReg::high_dim_feasible_estimate, X = dat$X,
#                                                     y = dat$y, K_range = c(1:5), delta = delta,
#                                                     max_iter = 10, cv_verbose = F)

method = cpReg::high_dim_feasible_estimate
X = dat$X
y = dat$y
K_range = c(1:5)
max_iter = 10
cv_verbose = F

dat <- list(X = X, y = y)
paramMat <- .tuning_lambda_K_pairing(method, dat, K_range = K_range,
                                     cv_verbose = cv_verbose, max_iter = max_iter,
                                     delta = delta)
#
# cv_val <- sapply(1:nrow(paramMat), function(x){
#   print(x)
#   if(cv_verbose) print(paste0("On K=", paramMat[x,2]))
#   .cross_validate(method = method, dat, K = paramMat[x,2],
#                   lambda = paramMat[x,1], cv_verbose = cv_verbose, delta = delta)
# })

.cross_validate(method = method, dat, K = paramMat[1,2],
                lambda = paramMat[1,1], cv_verbose = cv_verbose, delta = delta)

#################

K = paramMat[1,2]
lambda = paramMat[1,1]
nfolds = 5

n <- nrow(dat$X)
fold_id <- rep(1:nfolds, times = ceiling(n/nfolds))[1:n]

if(cv_verbose) print("Fitting each model")
res_list <- lapply(1:nfolds, function(i){
  idx <- which(fold_id != i)
  method(dat$X[idx,,drop = F], dat$y[idx], K = K, lambda = lambda, delta = delta)
})

if(cv_verbose) print("Evaluating the error")
# stats::median(sapply(1:nfolds, function(i){
#   print(i)
#   .out_of_sample_prediction(dat, fold_id, fold = i, res_list[[i]])
# }))


# .out_of_sample_prediction(dat, fold_id, fold = 5, res_list[[5]])

##########
fold = 5
res = res_list[[5]]

stopifnot(length(res$partition) > 2)

n <- nrow(dat$X)
idx <- which(fold_id == fold)
cp_idx <- .convert_cp_idx(res$partition, fold_id, fold)
stopifnot(cp_idx[1] == 0 & cp_idx[length(cp_idx)] == n)
