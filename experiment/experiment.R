y = 1
set.seed(y)
vec = paramMat[1,]
dat = rule(vec)

true_beta <- create_coef(vec, full = T)

n <- nrow(dat$X); p <- ncol(dat$X)

delta <- max(round(n/20), 5)

# parameter for feasible (and others)
set.seed(10)
# feasible_paramMat <- cpReg::tuning_cross_validation(cpReg::high_dim_feasible_estimate, X = dat$X,
#                                                     y = dat$y, K_range = c(1:5), delta = delta,
#                                                     max_iter = 10, cv_verbose = T)

##########

method = cpReg::high_dim_feasible_estimate
K_range = c(1:5)
cv_verbose = T
max_iter = 10
set.seed(10)
# paramMat <- cpReg:::.tuning_lambda_K_pairing(method, dat, K_range = K_range,
#                                      cv_verbose = cv_verbose, max_iter = max_iter,
#                                      delta = delta)

#################
set.seed(10)
n <- nrow(dat$X)

# first fit the entire dataset via lasso
fit <- glmnet::cv.glmnet(dat$X, dat$y, intercept = F, grouped = F)
new_lambda_vec <- rep(cpReg:::.glmnet_to_cp(fit$lambda.min, n), length(K_range))
lambda_vec <- rep(Inf, length(K_range))
iter <- 1

# iterate between lambda and finding partition
# while(iter < max_iter & sum(abs(new_lambda_vec - lambda_vec)) > 1e-3){
#   if(cv_verbose) print(paste0("Lambda values: ",
#                               paste0(round(new_lambda_vec, 2), collapse = ", ")))
#
#   lambda_vec <- new_lambda_vec
#
#   # try fitting the method now for different values of K
#   res_list <- lapply(1:length(K_range), function(i){
#     method(dat$X, dat$y, K= K_range[i], lambda = lambda_vec[i], delta = delta)
#   })
#
#   # based on the estimated partitions, refit lambda
#   new_lambda_vec <- sapply(1:length(res_list), function(i){
#     cpReg::oracle_tune_lambda(dat$X, dat$y, res_list[[i]]$partition/n, lambda.min = T)
#   })
#
#   iter <- iter + 1
# }

###############

lambda_vec <- new_lambda_vec

# try fitting the method now for different values of K
res_list <- lapply(1:length(K_range), function(i){
  print(i)
  method(dat$X, dat$y, K= K_range[i], lambda = lambda_vec[i], delta = delta)
})

new_lambda_vec <- sapply(1:length(res_list), function(i){
  print(i)
  cpReg::oracle_tune_lambda(dat$X, dat$y, res_list[[i]]$partition/n, lambda.min = T)
})

cpReg::oracle_tune_lambda(dat$X, dat$y, res_list[[4]]$partition/n, lambda.min = T)

##########
partition <- res_list[[4]]$partition/n
stopifnot(partition[1] == 0, partition[length(partition)] == 1)

X <- dat$X
y <- dat$y
n <- nrow(X)
lambda.min = T
k <- length(partition)-1
partition_idx <- round(partition*n)

stats::median(sapply(1:k, function(x){
  print(x)
  #if(partition_idx[x+1] - (partition_idx[x]+1) < 10) return(NA)
  fit <- glmnet::cv.glmnet(X[(partition_idx[x]+1):partition_idx[x+1],,drop = F],
                           y[(partition_idx[x]+1):partition_idx[x+1]],
                           intercept = F, grouped = F)
  val <- ifelse(lambda.min, fit$lambda.min, fit$lambda.1se)
  cpReg:::.glmnet_to_cp(val, partition_idx[x+1]-partition_idx[x])
}), na.rm = T)
