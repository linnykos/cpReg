rm(list=ls())
library(simulation)
library(cpReg)

set.seed(10)
trials <- 25
M <- 50
nu <- rep(0, M)
A <- matrix(0, M, M)
idx <- sample(M^2, 10)
A[idx] <- runif(length(idx), -.7,.3)

paramMat <- matrix(c(100, 150, 200, 250, 300, 350, 400), ncol = 1)
colnames(paramMat) <- "TT"

############

rule <- function(vec){
  cpReg::generative_model(nu, A, vec["TT"], lag = 1, thres_u = 6)
}

criterion <- function(dat, vec, y){
  cpReg::stationary_ar(dat, thres_u = 6,
                       basis_function = construct_AR_basis,
                       lambda = NA, verbose = F, lag = 1)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

###############

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 10, as_list = T,
                                        filepath = "../results/SEPP_tmp.RData",
                                        verbose = T)

save.image("../results/SEPP_reproduce_simulation.RData")
