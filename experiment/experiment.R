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
vec = paramMat[10,]
dat = rule(vec)

true_beta <- create_coef(vec, full = T)

lambda <- cpReg::oracle_tune_lambda(dat$X, dat$y, true_partition)
delta <- max(round(vec["n"]/10), 10)
k <- length(true_partition)-1

# cpReg::oracle_tune_gamma_range(dat$X, dat$y, lambda = lambda, k = k, delta = delta,
#                                verbose = T, max_iter = 10)

#############

X = dat$X; y = dat$y; verbose = T; max_iter = 10
min_gamma = 0.01; max_gamma = 1000; max_iter = 10
# min_gamma <- .initial_gamma_overshoot(X, y, lambda, k, min_gamma, delta = delta, smaller = T,
#                                       verbose = verbose)

##############

gamma = min_gamma; smaller = T
# res_low <- high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = gamma,
#                                       delta = delta)
data <- list(X = X, y = y)
# partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
#                  compute_cusum_func = .compute_regression_buhlmann,
#                  tau_function = .buhlmann_threshold,
#                  M = 0, delta = delta, max_candidates = 10,
#                  verbose = verbose,
#                  lambda = lambda, gamma = gamma)

##################

data_length_func = function(x){nrow(x$X)}
compute_cusum_func = .compute_regression_buhlmann
tau_function = .buhlmann_threshold
M = 0
max_candidates = 10

# initialize
n <- data_length_func(data)
random_intervals <- .generate_intervals(n, M)
q <- dequer::queue()
dequer::pushback(q, c(0,n))
b_vec <- numeric(0)
counter <- 0

# iterate over all intervals recursively
while(length(q) > 0){
  if(verbose) print(counter)
  interval <- dequer::pop(q)
  if(interval[2] - interval[1] < 2*delta+1) next()
  interval_list <- .truncate(interval, random_intervals)
  res_list <- lapply(interval_list, function(x){
    .find_breakpoint(data = data, interval = x, delta = delta,
                     max_candidates = max_candidates,
                     data_length_func = data_length_func,
                     compute_cusum_func = compute_cusum_func,
                     verbose = verbose, lambda = lambda, gamma = gamma)}) # error here
  tmp <- sapply(res_list, function(x){x$val})
  if(all(is.na(tmp))) break()
  res <- res_list[[which.max(tmp)]]

  # if passes threshold, recurse
  if(res$val >= tau_function(data = data, interval = interval, lambda = lambda, gamma = gamma)){
    b_vec <- c(b_vec, res$b)
    if(verbose) print(paste0("Pushing ", interval[1], " - ", res$b, " - ", interval[2]))
    dequer::pushback(q, c(interval[1], res$b))
    dequer::pushback(q, c(res$b+1, interval[2]))
  }

  counter <- counter+1
}

# .find_breakpoint(data = data, interval = interval_list[[1]], delta = delta,
#                  max_candidates = max_candidates,
#                  data_length_func = data_length_func,
#                  compute_cusum_func = compute_cusum_func,
#                  verbose = verbose, lambda = lambda, gamma = gamma)

##############

interval = interval_list[[1]]
stopifnot(interval[2] > interval[1])
if(interval[2] - interval[1] < 2*delta+1) return(list(val = NA, b = NA))
n <- data_length_func(data)

seq_vec <- (interval[1]+delta):(interval[2]-1-delta)

# if too many candidates, cut in half
if(!is.na(max_candidates) && length(seq_vec) > max_candidates){
  seq_vec <- sort(intersect(seq_vec, seq(delta, n - delta, by = delta)))
}

vec <- sapply(seq_vec, function(x){
  compute_cusum_func(data, interval[1], interval[2], x, lambda = lambda, gamma = gamma)
})
