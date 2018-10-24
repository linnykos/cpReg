# uncomment the following line to install
# devtools::install_github("linnylin92/cpReg", subdir = "cpReg")

# load library
library(cpReg)

## simple example to demonstrate how to set up the algorithm
# set up the data
set.seed(10)
M <- 5; TT <- 1000
nu <- rep(0.1, M)
A <- 0.3*diag(M)
dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = 5)

# fit the model
fit <- cpReg::stationary_ar(dat, thres_u = 5, basis_function = construct_AR_basis,
                            lambda = NA, verbose = T, lag = 1)

# compute the Forbenius difference between fitted A and true A
sum((A - fit$A)^2)

###############

## do a larger example
set.seed(10)
M <- 50; TT <- 500
nu <- rep(0, M)
# form block-wise A
A <- matrix(0, M, M)
# A[sample(1:prod(dim(A)), 0.5*prod(dim(A)))] <- 0
k <- 5; M_per_k <- round(M/k)
for(i in 1:k){
  tmp <- matrix(-.5, M_per_k, M_per_k)
  tmp[sample(1:prod(dim(tmp)), round(0.5*prod(dim(tmp))))] <- 0
  A[((i-1)*M_per_k+1):(i*M_per_k),((i-1)*M_per_k+1):(i*M_per_k)] <- tmp
}

dat <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = Inf)

# fit the model
fit <- cpReg::stationary_ar(dat, thres_u = Inf, basis_function = construct_AR_basis,
                            lambda = NA, sparsity = length(which(A == 0))/prod(dim(A)),
                            verbose = T, lag = 1, intercept = F)
length(which(fit$A == 0))/prod(dim(fit$A))

# visualize the sparsity patterns
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

par(mfrow = c(1,2))
tmp <- A; tmp[tmp != 0] <- 1
image(.rotate(tmp), breaks = c(-.5,.5,1.5), col = c("gray", "red"), asp = T)

tmp <- fit$A; tmp[tmp != 0] <- 1
image(.rotate(tmp), breaks = c(-.5,.5,1.5), col = c("gray", "red"), asp = T)
quantile(fit$A[fit$A != 0])

