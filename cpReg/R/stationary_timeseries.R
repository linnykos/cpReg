# Poisson time series, penality is L_{1,1}
stationary_ar <- function(dat, thres_u = round(quantile(dat[dat > 0], probs = 0.75)), lambda = NA,
                  basis_function = construct_AR_basis,
                  verbose = F, ...){

  stopifnot(is.matrix(dat), nrow(dat) > 1)
  stopifnot(all(dat >= 0), all(dat %% 1 == 0))
  stopifnot(any(dat < thres_u))

  M <- ncol(dat); TT <- nrow(dat)
  transform_dat <- construct_AR_basis(dat, thres_u = thres_u, ...)
  stopifnot(nrow(dat) == nrow(transform_dat), is.matrix(transform_dat))

  # perform row-wise glmnets
  if(verbose) print("Starting to estimate coefficients row-wise: ")
  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, is.na(lambda))
  })

  # if lambda is NA, combine the results together
  if(is.na(lambda)) {
    print("Starting to combine across different lambdas")
    est <- .combine_lambdas(res_list, dat, transform_dat, verbose = verbose)
  } else {
    est <- .extract_lambdas(res_list, lambda)
  }

  obj_val <- .objective_func(est$nu, est$A, dat, transform_dat)
  stopifnot(nrow(est$A) == nrow(dat), ncol(est$A) == ncol(transform_dat),
            length(est$nu) == ncol(dat))
  structure(list(lambda = est$lambda, obj_val = obj_val, nu = est$nu, A = est$A), class = "cp_ar1")
}

generative_model <- function(nu, A, timesteps = 10, thres_u = 5,
                             basis_function = construct_AR_basis,
                             warning_val = 1000, ...){
  stopifnot(nrow(A) == length(nu))

  M <- length(nu); MK <- ncol(A); TT <- timesteps
  dat <- matrix(0, ncol = M, nrow = TT)

  for(time in 1:TT){
    obs_vec <- sapply(1:M, function(i){
      vec <- construct_AR_basis(dat[which(1:TT < time),,drop = F], thres_u = thres_u, ...)
      vec <- vec[nrow(vec),]
      if(!is.na(warning_val) & as.numeric(nu[i] + A[i,] %*% vec) > log(warning_val)){
        stop(paste0("Natural parameters are exploding, with values above log of ", warning_val))
      }
      stats::rpois(1, lambda = exp(as.numeric(nu[i] + A[i,] %*% vec)))
    })

    dat[time,] <- obs_vec
  }

  dat
}

########

.objective_func <- function(nu, A, dat, transform_dat, lambda){
  .nll(nu, A, dat, transform_dat) + lambda *sum(abs(A))
}

.objective_func_row <- function(nu_val, A_vec, dat_vec, transform_dat, lambda){
  .nll_row(nu_val, A_vec, dat_vec, transform_dat) + lambda * sum(abs(A_vec))
}

.nll <- function(nu, A, dat, transform_dat){
  stopifnot(nrow(dat) == nrow(transform_dat), nrow(A) == ncol(dat),
            ncol(A) == ncol(transform_dat))
  TT <- nrow(dat)
  intercept_mat <- t(sapply(1:(TT-1), function(x){nu}))

  sum(exp(intercept_mat+ transform_dat[1:(TT-1),]%*%t(A))) - sum(dat[2:TT,] *(intercept_mat + transform_dat[1:(TT-1),]%*%t(A)))
}

# for each dimension
.nll_row <- function(nu_val, A_vec, dat_vec, transform_dat){
  stopifnot(length(dat_vec) == nrow(transform_dat))
  TT <- length(dat_vec)

  sum(exp(nu_val + A_vec %*% t(transform_dat[1:(TT-1),]))) - sum(dat_vec[2:TT] * (nu_val + A_vec %*% t(transform_dat[1:(TT-1),])))
}

# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#poi
.rowwise_glmnet <- function(response, covariates, cv = TRUE){
  stopifnot(length(response) == nrow(covariates))
  TT <- length(response)

  if(cv){
    if(TT < 30) stop("Not enough observations to perform cross validation")
    glmnet::cv.glmnet(x = covariates[1:(TT-1),], y = response[2:TT], family = "poisson", alpha = 1,
                      nfolds = min(max(floor(TT/10), 3), 10))

  } else {
    glmnet::glmnet(x = covariates[1:(TT-1),], y = response[2:TT], family = "poisson", alpha = 1)
  }
}
