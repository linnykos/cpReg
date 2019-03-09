set.seed(1)
vec <- paramMat[1,]
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

res1 <- cpReg::high_dim_feasible_estimate(dat$X, dat$y, lambda = lambda, tau = tau,
                                          verbose = F, max_candidates = NA, delta = delta, M = 0)

##########################
X = dat$X; y = dat$y
M = 100
delta = 10; max_candidates = NA; verbose = F

tau_function <- function(data, interval, ...){
  tau
}

M <- 0
data <- list(X = X, y = y)
# partition <- wbs(data, data_length_func = function(x){nrow(x$X)},
#                  compute_cusum_func = .compute_regression_cusum,
#                  tau_function = tau_function, M = M, delta = delta, max_candidates = max_candidates,
#                  verbose = verbose, lambda = lambda)
data_length_func = function(x){nrow(x$X)}
compute_cusum_func = .compute_regression_cusum

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
                     verbose = verbose, lambda = lambda)})
  res <- res_list[[which.max(sapply(res_list, function(x){x$val}))]]

  # if passes threshold, recurse
  if(res$val >= tau_function(data = data, interval = interval)){
    b_vec <- c(b_vec, res$b)
    if(verbose) print(paste0("Pushing ", interval[1], " - ", res$b, " - ", interval[2]))
    dequer::pushback(q, c(interval[1], res$b))
    dequer::pushback(q, c(res$b+1, interval[2]))
  }

  counter <- counter+1
}

#########
# investigating the .find_breakpoint...
interval = interval_list[[1]]
stopifnot(interval[2] > interval[1])
if(interval[2] - interval[1] < 2*delta+1) return(list(val = NA, b = NA))
n <- data_length_func(data)

seq_vec <- (interval[1]+delta):(interval[2]-1-delta)

# if too many candidates, cut in half
if(!is.na(max_candidates) && length(seq_vec) > max_candidates){
  seq_vec <- sort(intersect(seq_vec, seq(delta, n - delta, by = delta)))
}



