rm(list=ls())
set.seed(10)
M <- 5; TT <- 5000
nu <- rep(0, M)
A_list <- list(0.75*diag(M),
               matrix(0, M, M),
               0.5*matrix(runif(M*M), M, M))
changepoint_perc <- c(0.3, 0.6)

obj <- generative_model_cp(nu, A_list, changepoint_perc, timesteps = TT)

lambda <- .lambda_oracle(obj, intercept = F)

###########

min_spacing = 500
thres_u = round(stats::quantile(obj$dat[obj$dat > 0], probs = 0.75))
basis_function = construct_AR_basis
intercept = T

TT <- nrow(obj$dat)
changepoint_idx <- c(0, obj$partition, nrow(obj$dat))
k <- length(changepoint_idx)-1

# fit one big model
obj1 <- stationary_ar(obj$dat, thres_u = thres_u, lambda = lambda*sqrt(TT),
                      basis_function = basis_function, intercept = intercept)$obj_val

# fit the correct model
obj2 <- sum(sapply(1:k, function(i){
  len <- changepoint_idx[i+1]-changepoint_idx[i]

  stationary_ar(obj$dat[(changepoint_idx[i]+1):changepoint_idx[i+1],],
                thres_u = thres_u, lambda = lambda*sqrt(len),
                basis_function = basis_function, intercept = intercept)$obj_val
}))

if(obj1 < obj2) stop()
upper_val <- ((obj1-obj2)/(k-1))/TT

# fit the overfitted model
changepoint_idx3 <- sort(unique(unlist(lapply(1:k, function(i){
  len <- changepoint_idx[i+1]-changepoint_idx[i]
  if(len <= min_spacing) return(changepoint_idx[i:(i+1)])

  tmp <- seq(changepoint_idx[i], changepoint_idx[i+1], by = min_spacing)
  tmp[length(tmp)] <- changepoint_idx[i+1]
  tmp
}))))

obj3 <- sum(sapply(1:(length(changepoint_idx3)-1), function(i){
  len <- changepoint_idx3[i+1]-changepoint_idx3[i]

  stationary_ar(obj$dat[(changepoint_idx3[i]+1):changepoint_idx3[i+1],],
                thres_u = thres_u, lambda = lambda*sqrt(len),
                basis_function = basis_function, intercept = intercept)$obj_val
}))

lower_val <- ((obj2-obj3)/(length(changepoint_idx3)-k))/TT
