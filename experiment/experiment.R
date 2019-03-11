rm(list=ls())
set.seed(10)
n <- 100
partition <- c(0, 0.5, 1)
dat <- create_data(list(c(10,10,10), c(-10,-10,-10)), round(partition*n))
lambda <- oracle_tune_lambda(dat$X, dat$y, partition)
k <- length(partition)-1
k_vec <- c(2:5)

#######

X <- dat$X; y <- dat$y
delta = 10
min_gamma = 0.01
max_gamma = 1000
verbose = F
res <- lapply(k_vec, function(i){
  print(i)
  oracle_tune_gamma_range(X, y, lambda, i, delta = delta, min_gamma = min_gamma,
                          max_gamma = max_gamma, verbose = F)
})

quality_vec <- lapply(1:length(res), function(i){
  print(i)
  if(any(is.na(res[[i]]$gamma))) return(NA)

  gamma <- mean(res[[i]]$gamma)
  tmp <- high_dim_buhlmann_estimate(X, y, lambda = lambda, gamma = gamma,
                             delta = delta)$partition
  cpReg::hausdorff(tmp, round(partition*n))
})
