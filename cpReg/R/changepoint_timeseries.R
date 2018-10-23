# Poisson time series, penality is L_{1,1}
cp_ar <- function(dat, thres_u = round(quantile(dat[dat > 0], probs = 0.75)), lambda = NA,
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
    .rowwise_glmnet(dat[,x], transform_dat[,x])
  })

  # if lambda is NA, combine the results together
  if(is.na(lambda)) {
    print("Starting to combine across different lambdas")
    est <- .combine_lambdas(res_list, dat, transform_dat, verbose = verbose)
  } else {
    est <- .extract_lambdas(res_list, lambda)
  }

  stopifnot(nrow(A) == nrow(dat), ncol(A) == ncol(transform_dat))
  structure(list(nu = est$nu, A = est$A), class = "cp_ar1")
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
  intercept_mat <- sapply(1:TT, function(x){nu})

  sum(exp(intercept_mat+ A%*%t(transform_dat))) - sum(dat*(intercept_mat + A%*%t(transfrom_dat)))
}

# for each dimension
.nll_row <- function(nu_val, A_vec, dat_vec, transform_dat){
  sum(exp(nu_val + A_vec %*% t(transform_dat_vec))) - sum(dat_vec * (nu_val + A_vec %*% t(transform_dat_vec)))
}

# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#poi
# TODO: Check if the intercept is penalized
.rowwise_glmnet <- function(response, covariates, lambda = NA){
  if(is.na(lambda)){
    glmnet::cv.glmnet(covariates, response, family = "poisson", alpha = 1)

  } else {
    res <- glmnet::glmnet(covariates, response, family = "poisson", alpha = 1)
    as.numeric(glmnet::coef.glmnet(res, s = lambda))
  }
}
