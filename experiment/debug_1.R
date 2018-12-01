set.seed(10)
M <- 5; TT <- 200
nu <- rep(0, M)
A_list <- list(0.75*diag(M),
               matrix(0, M, M),
               0.5*matrix(runif(M*M), M, M))
changepoint_perc <- c(0.4, 0.6)

obj <- generative_model_cp(nu, A_list, changepoint_perc, timesteps = TT)

lambda <- .lambda_oracle(obj, intercept = F)
gamma <- .gamma_oracle(obj, lambda, intercept = F)
res1 <- changepoint_dp(obj$dat, lambda = lambda, gamma = gamma,
                       min_spacing = 10, skip_interval = 20, verbose = T)

res2 <- changepoint_dp(obj$dat, lambda = lambda, gamma = 50000,
                       min_spacing = 10, skip_interval = 20, verbose = T)
