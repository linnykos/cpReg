rm(list=ls())
library(simulation)
library(cpReg)
filename <- "../results/high_dim_tune_gamma.RData"
filename_tmp <- "../results/high_dim_tune_gamma_tmp.RData"

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

  lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
  delta <- max(round(vec["n"]/20), 10)
  k <- length(true_partition)-1

  res <- oracle_tune_gamma_range(dat$X, dat$y, lambda = lambda, k = k, delta = delta,
                                 verbose = F, max_iter = 10)
  print(res)

  list(lambda = lambda, gamma = res$gamma, min_gamma = res$min_gamma,
       max_gamma = res$max_gamma)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[10,]), paramMat[10,], 1)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = 100,
                                        cores = 20, as_list = T,
                                        filepath = filename_tmp,
                                        verbose = T)
save.image(filename)

## res <- res[which(sapply(res, length) !=0)]; sapply(res, function(zz){sum(sapply(zz, function(y){!all(is.na(y$gamma))}))/length(zz)})
